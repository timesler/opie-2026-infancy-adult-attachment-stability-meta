# =============================================================================
# 11_CreateModeratorSummary.R
# =============================================================================
# Purpose: Create comprehensive moderator summary file with interpretable results
# Input:   Moderator analysis CSV and raw data files
# Output:  moderator_summary.txt with all moderator results and group counts
# =============================================================================

library(dplyr)

# Set working directory
if (interactive()) {
  script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
  setwd(script_dir)
}

data_dir <- "../Data"
output_dir <- "../Output"

cat("=============================================================================\n")
cat("CREATING COMPREHENSIVE MODERATOR SUMMARY\n")
cat("=============================================================================\n\n")

# Load moderator results
mod_summary_file <- file.path(output_dir, "moderator_analysis_summary.csv")
if (!file.exists(mod_summary_file)) {
  cat("Error: Moderator analysis summary file not found!\n")
  cat("Please run 09_ModeratorAnalysis.R first.\n")
  stop()
}

mod_all <- read.csv(mod_summary_file, stringsAsFactors = FALSE)

# Load raw data files to get group counts for binary moderators
load_data_with_counts <- function(analysis_type) {
  data_file <- file.path(data_dir, paste0("data_", analysis_type, "_es.csv"))
  if (!file.exists(data_file)) {
    return(NULL)
  }
  
  data <- read.csv(data_file, stringsAsFactors = FALSE)
  
  # Calculate counts for binary moderators
  counts <- list(
    year = list(n = nrow(data), range = paste0(min(data$year_continuous, na.rm=T), "-", max(data$year_continuous, na.rm=T))),
    country_usa = sum(data$country_binary == 0, na.rm = TRUE),
    country_other = sum(data$country_binary == 1, na.rm = TRUE),
    age_t1_range = paste0(min(data$age_t1, na.rm=T), "-", max(data$age_t1, na.rm=T), " mo"),
    age_t2_range = paste0(min(data$age_t2, na.rm=T), "-", max(data$age_t2, na.rm=T), " yr"),
    risk_community = sum(data$sample_risk_binary == 0, na.rm = TRUE),
    risk_atrisk = sum(data$sample_risk_binary == 1, na.rm = TRUE),
    prior_no = sum(data$prior_meta_binary == 0, na.rm = TRUE),
    prior_yes = sum(data$prior_meta_binary == 1, na.rm = TRUE),
    interval_range = paste0(min(data$assessment_interval, na.rm=T), "-", max(data$assessment_interval, na.rm=T), " mo"),
    published_no = sum(data$published_binary == 0, na.rm = TRUE),
    published_yes = sum(data$published_binary == 1, na.rm = TRUE)
  )
  
  return(counts)
}

# Create summary file
output_file <- file.path(output_dir, "moderator_summary.txt")
sink(output_file)

cat("=============================================================================\n")
cat("COMPREHENSIVE MODERATOR ANALYSIS SUMMARY\n")
cat("Infant/Child to Adult Attachment Meta-Analysis\n")
cat("=============================================================================\n\n")

# Get unique analyses
analyses <- unique(mod_all$Analysis)
analysis_names <- list(
  "2way_sis" = "2-Way Secure/Insecure (Fisher's z, r)",
  "2way_od" = "2-Way Organized/Disorganized (Fisher's z, r)",
  "3way" = "3-Way A/B/C (Cohen's kappa)",
  "4way" = "4-Way A/B/C/D (Cohen's kappa)",
  "ssp_only" = "SSP-Only (Fisher's z, r)"
)

for (analysis_type in analyses) {
  
  # Get data counts for this analysis
  data_counts <- load_data_with_counts(analysis_type)
  
  # Get moderators for this analysis
  mod_data <- mod_all %>% filter(Analysis == analysis_type)
  
  # Get sample count from data
  n_samples <- if (!is.null(data_counts)) data_counts$year$n else nrow(mod_data)
  
  cat(sprintf("\n%s\n", toupper(analysis_names[[analysis_type]])))
  cat(paste(rep("=", 80), collapse = ""), "\n")
  cat(sprintf("Samples: %d | Studies: (see main analysis summary)\n\n", n_samples))
  
  # Print each moderator with interpretation
  moderators_order <- c("year", "country", "age_t1", "age_t2", "t2_measure", 
                        "sample_risk", "prior_meta", "assessment_interval", "published")
  
  for (mod_name in moderators_order) {
    mod_row <- mod_data %>% filter(Moderator == mod_name)
    
    if (nrow(mod_row) == 0) next
    
    b <- mod_row$b_coefficient
    se <- mod_row$SE
    p <- mod_row$p_value
    ci_l <- mod_row$CI_lower
    ci_u <- mod_row$CI_upper
    
    # Significance marker
    sig <- ifelse(p < 0.001, "***",
                  ifelse(p < 0.01, "**",
                         ifelse(p < 0.05, "*",
                                ifelse(p < 0.10, "+", ""))))
    
    # Print moderator name and counts
    cat(sprintf("\n%d. %s%s\n", which(moderators_order == mod_name), 
                mod_name, ifelse(sig != "", paste0(" ", sig), "")))
    
    # Print group counts for binary moderators
    if (mod_name == "country" && !is.null(data_counts)) {
      cat(sprintf("   USA: n=%d, Other: n=%d\n", data_counts$country_usa, data_counts$country_other))
    } else if (mod_name == "sample_risk" && !is.null(data_counts)) {
      cat(sprintf("   Community: n=%d, At-risk: n=%d\n", data_counts$risk_community, data_counts$risk_atrisk))
    } else if (mod_name == "prior_meta" && !is.null(data_counts)) {
      cat(sprintf("   Not included: n=%d, Included: n=%d\n", data_counts$prior_no, data_counts$prior_yes))
    } else if (mod_name == "published" && !is.null(data_counts)) {
      cat(sprintf("   Unpublished/Other: n=%d, Published article: n=%d\n", 
                  data_counts$published_no, data_counts$published_yes))
    } else if (mod_name == "age_t1" && !is.null(data_counts)) {
      cat(sprintf("   Range: %s\n", data_counts$age_t1_range))
    } else if (mod_name == "age_t2" && !is.null(data_counts)) {
      cat(sprintf("   Range: %s\n", data_counts$age_t2_range))
    } else if (mod_name == "assessment_interval" && !is.null(data_counts)) {
      cat(sprintf("   Range: %s\n", data_counts$interval_range))
    }
    
    # Print statistics
    cat(sprintf("   b = %.4f (SE = %.4f)\n", b, se))
    cat(sprintf("   95%% CI: [%.4f, %.4f]\n", ci_l, ci_u))
    cat(sprintf("   p-value: %.4f %s\n", p, sig))
    
    # Direction of effect interpretation
    # Build interpretation string based on coefficient value and moderator type
    direction_text <- ""
    if (mod_name == "year") {
      if (b < 0) {
        direction_text <- sprintf("Recent publications show LOWER continuity (b=%.4f)", b)
      } else {
        direction_text <- sprintf("Recent publications show HIGHER continuity (b=%.4f)", b)
      }
    } else if (mod_name == "country") {
      if (b < 0) {
        direction_text <- sprintf("USA (coded 0) shows HIGHER continuity than Other countries (b=%.4f)", b)
      } else {
        direction_text <- sprintf("Other countries (coded 1) show HIGHER continuity than USA (b=%.4f)", b)
      }
    } else if (mod_name == "age_t1") {
      if (b < 0) {
        direction_text <- sprintf("Younger infants at T1 show HIGHER continuity (b=%.4f)", b)
      } else {
        direction_text <- sprintf("Older infants at T1 show HIGHER continuity (b=%.4f)", b)
      }
    } else if (mod_name == "age_t2") {
      if (b < 0) {
        direction_text <- sprintf("Younger adults at T2 show HIGHER continuity (b=%.4f)", b)
      } else {
        direction_text <- sprintf("Older adults at T2 show HIGHER continuity (b=%.4f)", b)
      }
    } else if (mod_name == "t2_measure") {
      if (b < 0) {
        direction_text <- sprintf("SSP (coded 0) shows HIGHER continuity than AQS (b=%.4f)", b)
      } else {
        direction_text <- sprintf("AQS (coded 1) shows HIGHER continuity than SSP (b=%.4f)", b)
      }
    } else if (mod_name == "sample_risk") {
      if (b < 0) {
        direction_text <- sprintf("Community samples (coded 0) show HIGHER continuity than at-risk (b=%.4f)", b)
      } else {
        direction_text <- sprintf("At-risk samples (coded 1) show HIGHER continuity than community (b=%.4f)", b)
      }
    } else if (mod_name == "prior_meta") {
      if (b < 0) {
        direction_text <- sprintf("Studies NOT in prior meta-analyses show HIGHER continuity (b=%.4f)", b)
      } else {
        direction_text <- sprintf("Studies included in prior meta-analyses show HIGHER continuity (b=%.4f)", b)
      }
    } else if (mod_name == "assessment_interval") {
      if (b < 0) {
        direction_text <- sprintf("Shorter follow-up intervals show HIGHER continuity (b=%.4f)", b)
      } else {
        direction_text <- sprintf("Longer follow-up intervals show HIGHER continuity (b=%.4f)", b)
      }
    } else if (mod_name == "published") {
      if (b < 0) {
        direction_text <- sprintf("Unpublished/other sources (coded 0) show HIGHER continuity (b=%.4f)", b)
      } else {
        direction_text <- sprintf("Published articles (coded 1) show HIGHER continuity (b=%.4f)", b)
      }
    }
    
    # Print direction text
    if (direction_text != "") {
      cat(sprintf("   Effect: %s\n", direction_text))
    }
    
    # Interpretation
    if (p < 0.05) {
      if (b > 0) {
        cat(sprintf("   ✓ SIGNIFICANT positive effect\n"))
      } else {
        cat(sprintf("   ✓ SIGNIFICANT negative effect\n"))
      }
    } else if (p < 0.10) {
      cat(sprintf("   ~ Marginal effect (p < .10)\n"))
    } else {
      cat(sprintf("   Non-significant\n"))
    }
  }
}

# Summary statistics
cat("\n\n")
cat("=============================================================================\n")
cat("SIGNIFICANCE LEGEND\n")
cat("=============================================================================\n")
cat("*** p < .001 (highly significant)\n")
cat("**  p < .01  (very significant)\n")
cat("*   p < .05  (significant)\n")
cat("+   p < .10  (marginally significant)\n\n")

cat("MODERATOR DEFINITIONS\n")
cat("=============================================================================\n")
cat("1. Publication Year: Continuous measure of publication year (centered)\n")
cat("                     Higher values = more recent publications\n\n")
cat("2. Country: Binary (0=USA, 1=Other countries)\n")
cat("            Negative b indicates USA shows higher continuity\n\n")
cat("3. Age T1: Continuous measure of infant/child age at T1 (months)\n")
cat("           Higher values = older infants at first assessment\n\n")
cat("4. Age T2: Continuous measure of adult age at T2 (years)\n")
cat("           Higher values = older adults at follow-up\n\n")
cat("5. T2 Measure: Binary (0=SSP, 1=AQS)\n")
cat("               Note: Limited variation - most use AAI for adult measure\n\n")
cat("6. Sample Risk: Binary (0=Community, 1=At-risk)\n")
cat("                At-risk includes parent or child risk factors\n")
cat("                Positive b indicates at-risk samples show higher continuity\n\n")
cat("7. Prior Meta: Binary (0=Not included, 1=Included in prior meta-analyses)\n")
cat("               Negative b indicates previously-included studies show lower effects\n\n")
cat("8. Assessment Interval: Continuous measure of time between T1 and T2 (months)\n")
cat("                       Higher values = longer follow-up period\n")
cat("                       Positive b indicates longer intervals show higher continuity\n\n")
cat("9. Published: Binary (0=Unpublished/Other, 1=Published article)\n")
cat("              Positive b indicates published studies show higher continuity\n\n")

cat("=============================================================================\n")
cat("End of Report\n")
cat("=============================================================================\n")

sink()

cat("✓ Moderator summary created: ", output_file, "\n\n")
