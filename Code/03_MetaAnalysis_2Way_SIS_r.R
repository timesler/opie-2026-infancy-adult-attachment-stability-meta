# =============================================================================
# 03_MetaAnalysis_2Way_SIS_r.R
# =============================================================================
# Purpose: Meta-analysis of 2-way Secure/Insecure attachment continuity
#          from infancy to adulthood using Pearson's r (Fisher's z)
# Input:   data_2way_sis_es.csv (with effect sizes)
# Output:  Forest plot, funnel plot, summary statistics, RData object
# =============================================================================

# Load required packages
library(robumeta)
library(metafor)
library(dplyr)

# Set working directory to script location
if (interactive()) {
  script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
  setwd(script_dir)
}

# Source utilities
source("Utilities/RVE.r")
source("Utilities/forest.col.r")

# Define paths
data_dir <- "../Data"
output_dir <- "../Output/2way_SIS"
figures_dir <- "../Figures"

# Create output directory if it doesn't exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

cat("=============================================================================\n")
cat("META-ANALYSIS: 2-WAY SECURE/INSECURE (Pearson's r)\n")
cat("=============================================================================\n\n")

# =============================================================================
# STEP 1: Load and prepare data
# =============================================================================

cat("STEP 1: Loading data...\n")

# Load data with effect sizes
data <- read.csv(file.path(data_dir, "data_2way_sis_es.csv"), stringsAsFactors = FALSE)

cat(sprintf("  - Loaded %d samples with effect sizes\n", nrow(data)))
cat(sprintf("  - From %d unique studies\n", length(unique(data$study_group))))

# Prepare data for RVE analysis
data <- data %>%
  mutate(
    # Create study number for clustering (use study_group)
    studynum = as.numeric(as.factor(study_group)),
    
    # Effect size and variance (already calculated as Fisher's z)
    es = fisher_z,
    var_es = var_z,
    
    # Sample size
    n_total = n,
    
    # Create moderator variables
    measure_type_binary = ifelse(measure_type == "SSP", 0, 1),  # 0=SSP, 1=AQS
    parent_risk_binary = ifelse(parent_risk == "Yes", 1, 0),
    child_risk_binary = ifelse(child_risk == "Yes", 1, 0),
    any_risk_binary = ifelse(parent_risk == "Yes" | child_risk == "Yes", 1, 0),
    sample_type_binary = ifelse(sample_type == "Risk", 1, 0),
    year_centered = year - mean(year, na.rm = TRUE),
    
    # Time lag between assessments
    time_lag = adult_age_yr - infant_age_mo / 12,
    
    # Study label for plots
    study_label = paste0(author_year, " (", sample_subsample, ")")
  ) %>%
  # Remove rows with missing effect sizes
  filter(!is.na(es) & !is.na(var_es)) %>%
  # Shorten author names by extracting first author and year
  mutate(
    # Extract everything before the year in parentheses
    first_author = gsub("^([^,(&]+).*", "\\1", author_year),
    # Extract year in parentheses
    year_part = gsub(".*\\((\\d{4})\\)", "\\1", author_year),
    # Combine to create clean label
    study_label_short = paste0(trimws(first_author), " (", year_part, ")")
  ) %>%
  select(-first_author, -year_part) %>%
  # Only add sample description if study has multiple samples
  group_by(study_group) %>%
  mutate(
    n_samples = n(),
    study_label = ifelse(n_samples > 1, 
                        paste0(study_label_short, " (", sample_subsample, ")"),
                        study_label_short)
  ) %>%
  ungroup() %>%
  select(-study_label_short, -n_samples)

cat(sprintf("  - %d samples with complete data for analysis\n", nrow(data)))

# Print descriptive statistics
cat("\nDescriptive Statistics:\n")
cat(sprintf("  Mean r: %.3f (SD = %.3f)\n", mean(data$r, na.rm = TRUE), sd(data$r, na.rm = TRUE)))
cat(sprintf("  Range: [%.3f, %.3f]\n", min(data$r, na.rm = TRUE), max(data$r, na.rm = TRUE)))
cat(sprintf("  Mean Fisher's z: %.3f (SD = %.3f)\n", mean(data$es, na.rm = TRUE), sd(data$es, na.rm = TRUE)))
# Calculate total N using max N from each independent study
total_n <- data %>% group_by(studynum) %>% summarise(max_n = max(n_total, na.rm = TRUE)) %>% pull(max_n) %>% sum()
cat(sprintf("  Total participants: %d\n", total_n))
cat(sprintf("  Mean time lag: %.1f years\n", mean(data$time_lag, na.rm = TRUE)))

# =============================================================================
# STEP 2: Primary RVE meta-analysis (overall effect)
# =============================================================================

cat("\n\nSTEP 2: Running primary RVE meta-analysis...\n")

# Run RVE model with correlated effects
rve_model <- robu(
  formula = es ~ 1,
  data = data,
  studynum = studynum,
  var.eff.size = var_es,
  rho = 0.8,  # Assumed correlation between effect sizes from same study
  small = TRUE,  # Small sample correction
  modelweights = "CORR"
)

# Print summary
cat("\n")
print(rve_model)

# =============================================================================
# STEP 3: Sensitivity analysis for rho
# =============================================================================

cat("\n\nSTEP 3: Sensitivity analysis for rho...\n")

rho_sensitivity <- sensitivity(rve_model)
print(rho_sensitivity)

# =============================================================================
# STEP 4: Publication bias test (Egger's regression)
# =============================================================================

cat("\n\nSTEP 4: Testing for publication bias...\n")

# Egger's regression test (regress ES on SE) - PUBLISHED STUDIES ONLY
# Note: Unpublished studies should not be included in Egger's test as
# they are not subject to publication bias in the traditional sense
data_published <- data %>% filter(published_binary == 1)

if (nrow(data_published) >= 3) {
  egger_test <- robu(
    formula = es ~ 1 + sqrt(var_es),
    data = data_published,
    studynum = studynum,
    var.eff.size = var_es,
    rho = 0.8,
    small = TRUE,
    modelweights = "CORR"
  )
  
  cat("\nEgger's Regression Test (Published Studies Only):\n")
  cat(sprintf("  N published studies: %d\n", nrow(data_published)))
  cat(sprintf("  Intercept: %.3f, p = %.3f\n", 
              egger_test$reg_table$b.r[1], 
              egger_test$reg_table$prob[1]))
  cat(sprintf("  SE predictor: %.3f, p = %.3f %s\n", 
              egger_test$reg_table$b.r[2], 
              egger_test$reg_table$prob[2],
              ifelse(egger_test$reg_table$prob[2] < 0.05, "**", "")))
} else {
  cat("\nEgger's Regression Test: Not enough published studies (N < 3)\n")
  egger_test <- NULL
}

# =============================================================================
# STEP 5: Load moderator results from comprehensive analysis
# =============================================================================
# Note: Moderator analyses with 9 moderators are now handled by separate
# script 09_ModeratorAnalysis.R for consistency across all analyses.
# Load results here for display in summary file.

cat("\n\nSTEP 5: Loading moderator results from comprehensive analysis...\n")

# =============================================================================
# STEP 6: Create forest plot
# =============================================================================

cat("\n\nSTEP 6: Creating forest plot...\n")

# Prepare data for forest plot - using RVE weights for accuracy
# Match weights from RVE model by effect size values to handle potential reordering
model_weights <- rve_model$data.full$r.weights
model_effect_sizes <- rve_model$data.full$effect.size
data_effect_sizes <- data$es

# Match weights by effect size values (unique identifier for each observation)
matched_weights <- sapply(seq_len(nrow(data)), function(i) {
  matches <- which(abs(model_effect_sizes - data_effect_sizes[i]) < 1e-10)
  if (length(matches) >= 1) {
    return(model_weights[matches[1]])
  } else {
    return(1 / data$var_es[i])  # Fallback to inverse-variance
 }
})
matched_weights <- matched_weights / sum(matched_weights) * 100

data_forest <- data %>%
  mutate(
    # Convert back to r for display
    r_display = tanh(es),
    ci_lower = tanh(es - 1.96 * sqrt(var_es)),
    ci_upper = tanh(es + 1.96 * sqrt(var_es)),
    weight = matched_weights,
    study = gsub("^([^,(&]+).*?(\\(\\d{4}\\))", "\\1 \\2", author_year),
    sample = sample_subsample,
    n_total = n  # Copy n to n_total for consistency with header rows
  ) %>%
  # Calculate mean effect size per study for sorting
  group_by(study_group) %>%
  mutate(study_mean_r = mean(r_display)) %>%
  ungroup() %>%
  # Sort by study mean effect size
  arrange(study_mean_r, sample) %>%
  mutate(ID = row_number())

# Add header rows for studies with multiple samples (Opie2019 approach)
for (sg in unique(data_forest$study_group)) {
  data.grp <- data_forest[data_forest$study_group == sg, ]
  if (nrow(data.grp) > 1) {
    minID <- min(data.grp$ID)
    study_name <- data.grp$study[1]
    # Move study-level symbols (*/^) onto the sample labels for multi-sample
    # studies so the annotation is consistently shown on the sample, not the
    # study header (single-sample studies keep the symbol on the author).
    study_symbols <- regmatches(study_name, regexpr("[*^]+", study_name))
    if (length(study_symbols) == 0) study_symbols <- ""
    study_name <- trimws(gsub("[*^]+", "", study_name))
    # Indent sub-sample labels for multi-sample studies only (append symbols)
    indx <- data_forest$study_group == sg
    data_forest$study[indx] <- paste0("        ", data_forest$sample[indx], study_symbols)  # 8 spaces for sub-samples
    # Add header row just before first sub-sample
    header_row <- data.frame(
      study_group = sg,
      study = paste0("    ", study_name),  # 4 spaces for header
      ID = minID - 0.5,
      n_total = max(data.grp$n),  # Use max N to avoid double-counting subsamples
      weight = sum(data.grp$weight),
      r_display = NA,
      ci_lower = NA,
      ci_upper = NA,
      study_mean_r = data.grp$study_mean_r[1],
      stringsAsFactors = FALSE
    )
    data_forest <- bind_rows(data_forest, header_row)
  } else {
    # Single-sample studies: just add 4 spaces
    indx <- data_forest$study_group == sg
    data_forest$study[indx] <- paste0("    ", data_forest$study[indx])
  }
}

# Re-sort by ID to get headers in correct positions
data_forest <- data_forest %>% arrange(ID)

# Create forest plot
forest_file <- file.path(figures_dir, "forest_2way_sis.tiff")

# Calculate dimensions based on number of studies
n_studies <- nrow(data_forest)
plot_height <- 3 + 30 * n_studies / 75  # cm, adaptive height (increased from 25 to 30)

tiff(
  filename = forest_file,
  width = 18,
  height = plot_height,
  units = "cm",
  res = 600,
  compression = "lzw",
  pointsize = 6
)

# Split out header rows (those with NA for r_display) vs data rows
n_total_rows <- nrow(data_forest)
studyHeadingLocations <- (1:n_total_rows)[is.na(data_forest$r_display)]
sampleLocations <- (1:n_total_rows)[!is.na(data_forest$r_display)]
studyHeadings <- data_forest$study[is.na(data_forest$r_display)]
studyHeadingsN <- data_forest$n_total[is.na(data_forest$r_display)]
studyHeadingsWeight <- data_forest$weight[is.na(data_forest$r_display)]
data_for_plot <- data_forest[!is.na(data_forest$r_display), ]

# Set up plot parameters
cex <- 1.5
ylim <- c(-1, n_total_rows + 3)
xlim <- c(-5, 2.5)  # Wide range to prevent overlap
alim <- c(-1, 1)

# Prepare additional info columns (only for data rows)
ilab <- cbind(data_for_plot$n, format(round(data_for_plot$weight, 2), trim = FALSE))
ilab.xpos <- c(-1.8, -1.2)  # Position for N and Weight columns
ilab.pos <- 2
anno.pos <- 2.5  # Position for effect size and CI annotation

# Row positions: count from bottom (1) to top (n_total_rows)
# After rev(), top row (1) becomes position n_total_rows
row_positions_all <- n_total_rows:1
sampleRowPositions <- row_positions_all[sampleLocations]
headerRowPositions <- row_positions_all[studyHeadingLocations]

# Use custom forest.col function (only plot data rows, not headers)
forest.col(
  x = data_for_plot$r_display,
  vi = (data_for_plot$ci_upper - data_for_plot$ci_lower)^2 / (2 * 1.96)^2,
  xlim = xlim,
  alim = alim,
  ylim = ylim,
  digits = 2,
  xlab = "Pearson's r",
  slab = data_for_plot$study,
  ilab = ilab,
  ilab.xpos = ilab.xpos,
  ilab.pos = ilab.pos,
  rows = sampleRowPositions,
  efac = 0.4,
  cex = cex,
  refline = 0,
  anno.pos = anno.pos
)

# Add header rows manually (study names for multi-sample studies)
if (length(studyHeadingLocations) > 0) {
  text(
    x = xlim[1],
    y = headerRowPositions,
    labels = studyHeadings,
    cex = cex,
    pos = 4,
    font = 1
  )
  # Don't add N and Weight for header rows - leave them blank
}

# Add column labels
text(
  x = c(ilab.xpos, anno.pos),
  y = n_studies + 2,
  labels = c("N", "Weight", "Pearson's r [95% CI]"),
  cex = cex,
  pos = ilab.pos,
  font = 2
)

# Add pooled effect at bottom
pooled_r <- tanh(rve_model$reg_table$b.r)
pooled_ci_lower <- tanh(rve_model$reg_table$CI.L)
pooled_ci_upper <- tanh(rve_model$reg_table$CI.U)

addpoly.col(
  x = pooled_r,
  ci.lb = pooled_ci_lower,
  ci.ub = pooled_ci_upper,
  rows = -1,
  cex = cex,
  mlab = "Overall RE Model",
  efac = 0.4,
  digits = 2,
  font = 2,
  anno.pos = anno.pos
)

# Add total N and weight
# Sum N only from header rows (multi-sample studies) and data rows for single-sample studies
# Header rows have NA for r_display and contain max N per study
total_n_display <- sum(studyHeadingsN, na.rm = TRUE) + 
                   sum(data_for_plot$n[!data_for_plot$study_group %in% data_forest$study_group[is.na(data_forest$r_display)]], na.rm = TRUE)
text(
  x = c(ilab.xpos[1], ilab.xpos[2]),
  y = -1,
  labels = c(total_n_display, "100.00"),
  cex = cex,
  pos = ilab.pos,
  font = 2
)

dev.off()

cat(sprintf("  - Saved: %s\n", forest_file))

# =============================================================================
# STEP 7: Create funnel plot
# =============================================================================

cat("\n\nSTEP 7: Creating funnel plot...\n")

funnel_file <- file.path(figures_dir, "funnel_2way_sis.tiff")

tiff(funnel_file, width = 2400, height = 2400, res = 300)

# Create funnel plot with publication status
# Separate published and unpublished studies
data_pub <- data %>% filter(published_binary == 1)
data_unpub <- data %>% filter(published_binary == 0)

# Use Pearson's r (not Fisher's z) for x-axis
# SE for r can be approximated as: SE_r ≈ (1-r²)/sqrt(n-3)
se_vals <- sqrt((1 - data$r^2)^2 / (data$n_total - 3))
max_se <- max(se_vals, na.rm = TRUE)
x_range <- range(data$r, na.rm = TRUE)
x_margin <- diff(x_range) * 0.2
plot(
  x = data$r,
  y = se_vals,
  type = "n",
  xlim = c(x_range[1] - x_margin, x_range[2] + x_margin),
  ylim = rev(c(0, max_se * 1.1)),  # Reverse ylim so small SE at top
  xlab = "Pearson's r",
  ylab = "Standard Error",
  main = "Funnel Plot: 2-way S/IS Continuity"
)

# Add funnel (pseudo-confidence region) - convert pooled effect from z to r
pooled_effect <- tanh(rve_model$reg_table$b.r)
se_seq <- seq(0, max_se * 1.1, length.out = 100)
ci_lower <- pooled_effect - 1.96 * se_seq
ci_upper <- pooled_effect + 1.96 * se_seq
polygon(c(ci_lower, rev(ci_upper)), c(se_seq, rev(se_seq)), 
        col = "gray90", border = NA)

# Add reference line at pooled effect
abline(v = pooled_effect, lty = 2, col = "blue", lwd = 2)

# Plot points: Published (filled circles), Unpublished (open circles)
se_pub <- sqrt((1 - data_pub$r^2)^2 / (data_pub$n_total - 3))
se_unpub <- sqrt((1 - data_unpub$r^2)^2 / (data_unpub$n_total - 3))
if (nrow(data_pub) > 0) {
  points(data_pub$r, se_pub, pch = 19, col = "black", cex = 1.2)
}
if (nrow(data_unpub) > 0) {
  points(data_unpub$r, se_unpub, pch = 1, col = "black", cex = 1.2, lwd = 1.5)
}

# Add legend
legend("topright", 
       legend = c(sprintf("Published (n=%d)", nrow(data_pub)), 
                  sprintf("Unpublished (n=%d)", nrow(data_unpub))),
       pch = c(19, 1), 
       col = "black",
       pt.cex = 1.2,
       bty = "n")

dev.off()

cat(sprintf("  - Saved: %s\n", funnel_file))

# =============================================================================
# STEP 8: Export results
# =============================================================================

cat("\n\nSTEP 8: Exporting results...\n")

# Create summary table
summary_table <- data.frame(
  Analysis = "2-way S/IS",
  k_samples = nrow(data),
  k_studies = length(unique(data$studynum)),
  N_total = sum(data$n, na.rm = TRUE),
  Mean_r = mean(data$r, na.rm = TRUE),
  SD_r = sd(data$r, na.rm = TRUE),
  Pooled_z = rve_model$reg_table$b.r,
  Pooled_r = tanh(rve_model$reg_table$b.r),
  SE_z = rve_model$reg_table$SE,
  CI_lower_z = rve_model$reg_table$CI.L,
  CI_upper_z = rve_model$reg_table$CI.U,
  CI_lower_r = tanh(rve_model$reg_table$CI.L),
  CI_upper_r = tanh(rve_model$reg_table$CI.U),
  p_value = rve_model$reg_table$prob,
  tau_sq = rve_model$mod_info$tau.sq,
  I_sq = rve_model$mod_info$I.2,
  Egger_p = ifelse(!is.null(egger_test), egger_test$reg_table$prob[2], NA)
)

# Save summary table
write.csv(summary_table, file.path(output_dir, "summary_statistics.csv"), row.names = FALSE)
cat("  - Saved: summary_statistics.csv\n")

# Note: Moderator results are saved by script 09_ModeratorAnalysis.R
# This ensures consistent methodology across all analyses

# Save full results object
save(rve_model, rho_sensitivity, egger_test, data, 
     file = file.path(output_dir, "analysis_results.RData"))
cat("  - Saved: analysis_results.RData\n")

# Create text summary
sink(file.path(output_dir, "analysis_summary.txt"))

cat("=============================================================================\n")
cat("META-ANALYSIS: 2-WAY SECURE/INSECURE ATTACHMENT CONTINUITY\n")
cat("Infant/Child to Adult Attachment (SSP/AQS → AAI)\n")
cat("=============================================================================\n\n")

cat("SAMPLE CHARACTERISTICS\n")
cat("----------------------\n")
cat(sprintf("Number of samples: %d\n", nrow(data)))
cat(sprintf("Number of independent studies: %d\n", length(unique(data$studynum))))
# Calculate total N using max N from each independent study (to avoid double-counting subsamples)
total_n <- data %>% group_by(studynum) %>% summarise(max_n = max(n_total, na.rm = TRUE)) %>% pull(max_n) %>% sum()
cat(sprintf("Total participants: %d\n", total_n))
cat(sprintf("Mean time lag: %.1f years (SD = %.1f)\n\n", 
            mean(data$time_lag, na.rm = TRUE), 
            sd(data$time_lag, na.rm = TRUE)))

cat("OVERALL EFFECT\n")
cat("--------------\n")
cat(sprintf("Pooled r: %.3f [%.3f, %.3f]\n", 
            tanh(rve_model$reg_table$b.r),
            tanh(rve_model$reg_table$CI.L),
            tanh(rve_model$reg_table$CI.U)))
cat(sprintf("Fisher's z: %.3f (SE = %.3f)\n", 
            rve_model$reg_table$b.r,
            rve_model$reg_table$SE))
cat(sprintf("p-value: %.4f %s\n", 
            rve_model$reg_table$prob,
            ifelse(rve_model$reg_table$prob < 0.001, "***",
                   ifelse(rve_model$reg_table$prob < 0.01, "**",
                          ifelse(rve_model$reg_table$prob < 0.05, "*", "")))))
cat(sprintf("df: %.1f\n\n", rve_model$reg_table$dfs))

cat("HETEROGENEITY\n")
cat("-------------\n")
cat(sprintf("Tau-squared: %.4f\n", rve_model$mod_info$tau.sq))
cat(sprintf("I-squared: %.2f%%\n\n", rve_model$mod_info$I.2))

cat("PUBLICATION BIAS\n")
cat("----------------\n")
if (!is.null(egger_test)) {
  cat(sprintf("Egger's test p-value (published only): %.4f %s\n",
              egger_test$reg_table$prob[2],
              ifelse(egger_test$reg_table$prob[2] < 0.05, "*", "")))
  cat(sprintf("N published studies tested: %d\n\n", nrow(data_published)))
} else {
  cat("Egger's test: Not conducted (insufficient published studies)\n\n")
}

cat("SENSITIVITY ANALYSIS (rho)\n")
cat("--------------------------\n")
print(rho_sensitivity)
cat("\n")

cat("MODERATOR ANALYSES\n")
cat("------------------\n")
cat("Comprehensive moderator analyses (9 moderators) are reported separately in:\n")
cat("  ../Output/moderator_summary.txt\n\n")

cat("=============================================================================\n")

sink()

cat("  - Saved: analysis_summary.txt\n")

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n\n=============================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("=============================================================================\n\n")

cat(sprintf("Pooled effect: r = %.3f [%.3f, %.3f], p = %.4f\n",
            tanh(rve_model$reg_table$b.r),
            tanh(rve_model$reg_table$CI.L),
            tanh(rve_model$reg_table$CI.U),
            rve_model$reg_table$prob))

cat(sprintf("\nInterpretation: %s\n",
            ifelse(rve_model$reg_table$prob < 0.05,
                   "Significant positive continuity",
                   "No significant continuity")))

cat("\nAll outputs saved to:", output_dir, "\n")
cat("=============================================================================\n")
