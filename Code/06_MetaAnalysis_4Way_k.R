# =============================================================================
# 06_MetaAnalysis_4Way_k.R
# =============================================================================
# Purpose: Meta-analysis of 4-way (A/B/C) attachment continuity
#          from infancy to adulthood using Cohen's kappa
# Input:   data_4way_es.csv (with effect sizes)
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
source("Utilities/meta_analysis_helpers.r")

# Define paths
data_dir <- "../Data"
output_dir <- "../Output/4way"
figures_dir <- "../Figures"

# Create output directory if it doesn't exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

cat("=============================================================================\n")
cat("META-ANALYSIS: 4-WAY A/B/C/D (Cohen's kappa)\n")
cat("=============================================================================\n\n")

# =============================================================================
# STEP 1: Load and prepare data
# =============================================================================

cat("STEP 1: Loading data...\n")

# Load data with effect sizes
data <- read.csv(file.path(data_dir, "data_4way_es.csv"), stringsAsFactors = FALSE)

cat(sprintf("  - Loaded %d samples with effect sizes\n", nrow(data)))
cat(sprintf("  - From %d unique studies\n", length(unique(data$study_group))))

# Prepare data for RVE analysis
data <- data %>%
  mutate(
    # Create study number for clustering
    studynum = as.numeric(as.factor(study_group)),
    
    # Effect size and variance (Cohen's kappa)
    es = kappa,
    var_es = var_kappa,
    
    # Sample size
    n_total = n,
    
    # Create moderator variables
    measure_type_binary = ifelse(measure_type == "SSP", 0, 1),
    parent_risk_binary = ifelse(parent_risk == "Yes", 1, 0),
    child_risk_binary = ifelse(child_risk == "Yes", 1, 0),
    any_risk_binary = ifelse(parent_risk == "Yes" | child_risk == "Yes", 1, 0),
    sample_type_binary = ifelse(sample_type == "Risk", 1, 0),
    year_centered = year - mean(year, na.rm = TRUE),
    
    # Time lag between assessments
    time_lag = adult_age_yr - infant_age_mo / 12
  ) %>%
  # Remove rows with missing effect sizes
  filter(!is.na(es) & !is.na(var_es))

cat(sprintf("  - %d samples with complete data for analysis\n", nrow(data)))

# Print descriptive statistics
cat("\nDescriptive Statistics:\n")
cat(sprintf("  Mean kappa: %.3f (SD = %.3f)\n", mean(data$kappa, na.rm = TRUE), sd(data$kappa, na.rm = TRUE)))
cat(sprintf("  Range: [%.3f, %.3f]\n", min(data$kappa, na.rm = TRUE), max(data$kappa, na.rm = TRUE)))
cat(sprintf("  Total participants: %d\n", sum(data$n_total, na.rm = TRUE)))
cat(sprintf("  Mean time lag: %.1f years\n", mean(data$time_lag, na.rm = TRUE)))

# =============================================================================
# STEP 2: Primary RVE meta-analysis
# =============================================================================

cat("\n\nSTEP 2: Running primary RVE meta-analysis...\n")

rve_model <- robu(
  formula = es ~ 1,
  data = data,
  studynum = studynum,
  var.eff.size = var_es,
  rho = 0.8,
  small = TRUE,
  modelweights = "CORR"
)

cat("\n")
print(rve_model)

# =============================================================================
# STEP 3: Sensitivity analysis for rho
# =============================================================================

cat("\n\nSTEP 3: Sensitivity analysis for rho...\n")

rho_sensitivity <- sensitivity(rve_model)
print(rho_sensitivity)

# =============================================================================
# STEP 4: Publication bias test
# =============================================================================

cat("\n\nSTEP 4: Testing for publication bias...\n")

# Egger's regression test - PUBLISHED STUDIES ONLY
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

# Prepare data using helper function
data_forest <- prepare_forest_data(
  data = data,
  effect_col = "kappa",
  var_col = "var_kappa",
  rve_model = rve_model,
  convert_to_r = FALSE  # Don't convert kappa
)

# Create forest plot
forest_file <- file.path(figures_dir, "forest_4way.tiff")

create_forest_plot(
  data_forest = data_forest,
  output_file = forest_file,
  main_title = "4-Way A/B/C/D Attachment Continuity (Infancy to Adulthood)",
  xlab = "Cohen's kappa",
  xlim = c(-5, 2.5),
  alim = c(-0.5, 0.5)
)

# Add pooled effect - need to manually add since we saved the plot
# Re-open the file and add pooled effect
tiff(forest_file, width = 18, height = 3 + 30 * nrow(data_forest) / 75, 
     units = "cm", res = 600, compression = "lzw", pointsize = 6)

# Recreate the base plot first
n_total_rows <- nrow(data_forest)
data_for_plot <- data_forest[!is.na(data_forest$es_display), ]
sampleLocations <- (1:n_total_rows)[!is.na(data_forest$es_display)]
row_positions_all <- n_total_rows:1
sampleRowPositions <- row_positions_all[sampleLocations]

ilab <- cbind(data_for_plot$n, format(round(data_for_plot$weight, 2), trim = FALSE))
cex <- 1.5

forest.col(
  x = data_for_plot$es_display,
  vi = (data_for_plot$ci_upper - data_for_plot$ci_lower)^2 / (2 * 1.96)^2,
  xlim = c(-5, 2.5),
  alim = c(-0.5, 0.5),
  ylim = c(-1, n_total_rows + 3),
  digits = 2,
  xlab = "Cohen's kappa",
  main = "4-Way A/B/C/D Attachment Continuity (Infancy to Adulthood)",
  slab = data_for_plot$study,
  ilab = ilab,
  ilab.xpos = c(-1.8, -1.2),
  ilab.pos = 2,
  rows = sampleRowPositions,
  efac = 0.4,
  cex = cex,
  refline = 0,
  anno.pos = 2.5
)

# Add column headers
text(x = c(-1.9, -1.2, 2.5), y = n_total_rows + 2,
     labels = c("N", "Weight", "Cohen's kappa [95% CI]"),
     cex = cex, pos = 2, font = 2)

# Add headers
studyHeadingLocations <- (1:n_total_rows)[is.na(data_forest$es_display)]
if (length(studyHeadingLocations) > 0) {
  headerRowPositions <- row_positions_all[studyHeadingLocations]
  text(x = -5, y = headerRowPositions, 
       labels = data_forest$study[is.na(data_forest$es_display)],
       cex = cex, pos = 4, font = 1)
}

# Add pooled effect
addpoly.col(
  x = rve_model$reg_table$b.r,
  ci.lb = rve_model$reg_table$CI.L,
  ci.ub = rve_model$reg_table$CI.U,
  rows = -1,
  cex = cex,
  mlab = "Overall RE Model",
  efac = 0.4,
  digits = 2,
  font = 2,
  anno.pos = 2.5
)

# Add total N and weight
# Calculate total N using max N per independent study to avoid double-counting
total_n_overall <- data %>% 
  group_by(study_group) %>% 
  summarise(max_n = max(n, na.rm = TRUE)) %>% 
  pull(max_n) %>% 
  sum()
text(x = c(-1.9, -1.2), y = -1,
     labels = c(total_n_overall, "100.00"),
     cex = cex, pos = 2, font = 2)

dev.off()

cat(sprintf("  - Saved: %s\n", forest_file))

# =============================================================================
# STEP 7: Create funnel plot
# =============================================================================

cat("\n\nSTEP 7: Creating funnel plot...\n")

funnel_file <- file.path(figures_dir, "funnel_4way.tiff")

tiff(funnel_file, width = 2400, height = 2400, res = 300)

# Create funnel plot with publication status
data_pub <- data %>% filter(published_binary == 1)
data_unpub <- data %>% filter(published_binary == 0)

se_vals <- sqrt(data$var_es)
max_se <- max(se_vals)
x_range <- range(data$es)
x_margin <- diff(x_range) * 0.2
plot(
  x = data$es,
  y = se_vals,
  type = "n",
  xlim = c(x_range[1] - x_margin, x_range[2] + x_margin),
  ylim = rev(c(0, max_se * 1.1)),
  xlab = "kappa",
  ylab = "Standard Error",
  main = "Funnel Plot: 4-way A/B/C/D Continuity"
)

pooled_effect <- rve_model$reg_table$b.r
se_seq <- seq(0, max_se * 1.1, length.out = 100)
ci_lower <- pooled_effect - 1.96 * se_seq
ci_upper <- pooled_effect + 1.96 * se_seq
polygon(c(ci_lower, rev(ci_upper)), c(se_seq, rev(se_seq)),
        col = "gray90", border = NA)

abline(v = pooled_effect, lty = 2, col = "blue", lwd = 2)

if (nrow(data_pub) > 0) {
  points(data_pub$es, sqrt(data_pub$var_es), pch = 19, col = "black", cex = 1.2)
}
if (nrow(data_unpub) > 0) {
  points(data_unpub$es, sqrt(data_unpub$var_es), pch = 1, col = "black", cex = 1.2, lwd = 1.5)
}

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

# Summary table
summary_table <- data.frame(
  Analysis = "4-way",
  k_samples = nrow(data),
  k_studies = length(unique(data$studynum)),
  N_total = sum(data$n, na.rm = TRUE),
  Mean_kappa = mean(data$kappa, na.rm = TRUE),
  SD_kappa = sd(data$kappa, na.rm = TRUE),
  Pooled_kappa = rve_model$reg_table$b.r,
  SE_kappa = rve_model$reg_table$SE,
  CI_lower = rve_model$reg_table$CI.L,
  CI_upper = rve_model$reg_table$CI.U,
  p_value = rve_model$reg_table$prob,
  tau_sq = rve_model$mod_info$tau.sq,
  I_sq = rve_model$mod_info$I.2,
  Egger_p = ifelse(!is.null(egger_test), egger_test$reg_table$prob[2], NA)
)

write.csv(summary_table, file.path(output_dir, "summary_statistics.csv"), row.names = FALSE)
cat("  - Saved: summary_statistics.csv\n")

# Note: Moderator results are saved by script 09_ModeratorAnalysis.R
# This ensures consistent methodology across all analyses

# Full results
save(rve_model, rho_sensitivity, egger_test, data, 
     file = file.path(output_dir, "analysis_results.RData"))
cat("  - Saved: analysis_results.RData\n")

# Text summary
sink(file.path(output_dir, "analysis_summary.txt"))
cat("=============================================================================\n")
cat("META-ANALYSIS: 4-WAY A/B/C/D ATTACHMENT CONTINUITY\n")
cat("Infant/Child to Adult Attachment (SSP/AQS → AAI)\n")
cat("=============================================================================\n\n")

cat("SAMPLE CHARACTERISTICS\n")
cat("----------------------\n")
cat(sprintf("Number of samples: %d\n", nrow(data)))
cat(sprintf("Number of independent studies: %d\n", length(unique(data$studynum))))
cat(sprintf("Total participants: %d\n", sum(data$n_total, na.rm = TRUE)))
cat(sprintf("Mean time lag: %.1f years (SD = %.1f)\n\n", 
            mean(data$time_lag, na.rm = TRUE), 
            sd(data$time_lag, na.rm = TRUE)))

cat("OVERALL EFFECT\n")
cat("--------------\n")
cat(sprintf("Pooled kappa: %.3f [%.3f, %.3f]\n", 
            rve_model$reg_table$b.r,
            rve_model$reg_table$CI.L,
            rve_model$reg_table$CI.U))
cat(sprintf("SE: %.3f\n", rve_model$reg_table$SE))
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

cat("\n=============================================================================\n")
sink()
