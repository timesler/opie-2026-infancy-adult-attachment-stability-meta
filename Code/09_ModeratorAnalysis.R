# =============================================================================
# 09_ModeratorAnalysis.R
# =============================================================================
# Purpose: Comprehensive moderator analysis aligned with Opie2019 approach
# Tests 9 moderators across all meta-analyses
# Input:   data_*_es.csv files with effect sizes and moderator variables
# Output:  Moderator results tables (CSV), summary statistics
# =============================================================================

library(robumeta)
library(dplyr)
library(tidyr)
library(rlang)

# Set working directory to script location
if (interactive()) {
  script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
  setwd(script_dir)
}

data_dir <- "../Data"
output_dir <- "../Output"

cat("=============================================================================\n")
cat("COMPREHENSIVE MODERATOR ANALYSES\n")
cat("=============================================================================\n\n")

# Define moderators to test
# 1. Publication year (continuous)
# 2. Country (binary: USA vs. Other)
# 3. Age at T1 (continuous, in months)
# 4. Age at T2 (continuous, in years)
# 5. T2 measure type (binary: SSP vs AQS, if enough variation)
# 6. Sample risk (binary: at-risk vs. community)
# 7. Included in prior meta-analysis (binary: yes vs. no)
# 8. Assessment interval (continuous, in months)
# 9. Publication status (binary: published article vs. other)

moderators_list <- list(
  year = "year_continuous",
  country = "country_binary",
  age_t1 = "age_t1",
  age_t2 = "age_t2",
  t2_measure = "t2_measure",
  sample_risk = "sample_risk_binary",
  prior_meta = "prior_meta_binary",
  assessment_interval = "assessment_interval",
  published = "published_binary"
)

# =============================================================================
# Function to run moderator analyses for a dataset
# =============================================================================

run_moderator_analyses <- function(data, analysis_name, effect_size_var = "es", 
                                     variance_var = "var_es") {
  
  cat(sprintf("\n%s MODERATOR ANALYSES\n", toupper(analysis_name)))
  cat(paste(rep("-", 80), collapse = ""), "\n")
  
  # Ensure studynum is numeric
  if (!is.numeric(data$studynum)) {
    data$studynum <- as.numeric(as.factor(data$study_group))
  }
  
  # Initialize results matrix
  results <- data.frame(
    Moderator = character(),
    Variable = character(),
    N = integer(),
    b_coefficient = numeric(),
    SE = numeric(),
    df = numeric(),
    p_value = numeric(),
    CI_lower = numeric(),
    CI_upper = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Test each moderator
  for (mod_name in names(moderators_list)) {
    mod_var <- moderators_list[[mod_name]]
    
    # Check if moderator variable exists in data
    if (!mod_var %in% names(data)) {
      cat(sprintf("  ⚠ %s: Variable '%s' not found in dataset\n", mod_name, mod_var))
      next
    }
    
    # Check for missing values
    mod_data <- data %>% filter(!is.na(!!sym(mod_var)))
    
    if (nrow(mod_data) == 0) {
      cat(sprintf("  ⚠ %s: No valid data\n", mod_name))
      next
    }
    
    # For binary variables, check if there's variation
    if (mod_var %in% c("country_binary", "t2_measure", "sample_risk_binary", 
                        "prior_meta_binary", "published_binary")) {
      unique_vals <- length(unique(mod_data[[mod_var]]))
      if (unique_vals < 2) {
        cat(sprintf("  ⚠ %s: No variation in data (n_categories = %d)\n", mod_name, unique_vals))
        next
      }
    }
    
    # Handle t2_measure special case - convert to binary if needed
    if (mod_var == "t2_measure") {
      if (is.character(mod_data[[mod_var]])) {
        mod_data <- mod_data %>% 
          mutate(t2_measure = ifelse(t2_measure == "SSP", 0, 1))
      }
    }
    
    # Run RVE model with moderator
    tryCatch({
      # Build formula string
      formula_str <- paste(effect_size_var, "~ 1 +", mod_var)
      
      # Fit model
      mod_model <- robu(
        formula = as.formula(formula_str),
        data = mod_data,
        studynum = mod_data$studynum,
        var.eff.size = mod_data[[variance_var]],
        rho = 0.8,
        small = TRUE,
        modelweights = "CORR"
      )
      
      # Extract results
      mod_coef <- mod_model$reg_table$b.r[2]
      mod_se <- mod_model$reg_table$SE[2]
      mod_df <- mod_model$reg_table$dfs[2]
      mod_p <- mod_model$reg_table$prob[2]
      mod_ci_l <- mod_model$reg_table$CI.L[2]
      mod_ci_u <- mod_model$reg_table$CI.U[2]
      
      # Significance indicators
      sig_marker <- ifelse(mod_p < 0.001, "***",
                          ifelse(mod_p < 0.01, "**",
                                ifelse(mod_p < 0.05, "*", "")))
      
      cat(sprintf("  %s %-25s b=%.4f, SE=%.4f, p=%.4f %s\n", 
                  sig_marker, mod_name, mod_coef, mod_se, mod_p, 
                  ifelse(sig_marker != "", "←", "")))
      
      # Add to results
      results <- rbind(results, data.frame(
        Moderator = mod_name,
        Variable = mod_var,
        N = nrow(mod_data),
        b_coefficient = round(mod_coef, 6),
        SE = round(mod_se, 6),
        df = round(mod_df, 4),
        p_value = round(mod_p, 6),
        CI_lower = round(mod_ci_l, 6),
        CI_upper = round(mod_ci_u, 6),
        stringsAsFactors = FALSE
      ))
      
    }, error = function(e) {
      cat(sprintf("  ✗ %s: Error - %s\n", mod_name, e$message))
    })
  }
  
  return(results)
}

# =============================================================================
# Run moderator analyses for each analysis type
# =============================================================================

all_results <- list()

# 2-way S/IS
if (file.exists(file.path(data_dir, "data_2way_sis_es.csv"))) {
  data_2way_sis <- read.csv(file.path(data_dir, "data_2way_sis_es.csv"))
  data_2way_sis <- data_2way_sis %>%
    mutate(
      studynum = as.numeric(as.factor(study_group)),
      es = fisher_z,
      var_es = var_z
    )
  all_results$"2way_sis" <- run_moderator_analyses(data_2way_sis, "2-way S/IS", 
                                                     effect_size_var = "es",
                                                     variance_var = "var_es")
}

# 2-way O/D
if (file.exists(file.path(data_dir, "data_2way_od_es.csv"))) {
  data_2way_od <- read.csv(file.path(data_dir, "data_2way_od_es.csv"))
  data_2way_od <- data_2way_od %>%
    mutate(
      studynum = as.numeric(as.factor(study_group)),
      es = fisher_z,
      var_es = var_z
    )
  all_results$"2way_od" <- run_moderator_analyses(data_2way_od, "2-way O/D",
                                                    effect_size_var = "es",
                                                    variance_var = "var_es")
}

# 3-way
if (file.exists(file.path(data_dir, "data_3way_es.csv"))) {
  data_3way <- read.csv(file.path(data_dir, "data_3way_es.csv"))
  data_3way <- data_3way %>%
    mutate(
      studynum = as.numeric(as.factor(study_group)),
      es = kappa,
      var_es = var_kappa
    )
  all_results$"3way" <- run_moderator_analyses(data_3way, "3-way",
                                                effect_size_var = "es",
                                                variance_var = "var_es")
}

# 4-way
if (file.exists(file.path(data_dir, "data_4way_es.csv"))) {
  data_4way <- read.csv(file.path(data_dir, "data_4way_es.csv"))
  data_4way <- data_4way %>%
    mutate(
      studynum = as.numeric(as.factor(study_group)),
      es = kappa,
      var_es = var_kappa
    )
  all_results$"4way" <- run_moderator_analyses(data_4way, "4-way",
                                                effect_size_var = "es",
                                                variance_var = "var_es")
}

# SSP-only
if (file.exists(file.path(data_dir, "data_ssp_only_es.csv"))) {
  data_ssp <- read.csv(file.path(data_dir, "data_ssp_only_es.csv"))
  data_ssp <- data_ssp %>%
    mutate(
      studynum = as.numeric(as.factor(study_group)),
      es = fisher_z,
      var_es = var_z
    )
  all_results$"ssp_only" <- run_moderator_analyses(data_ssp, "SSP-only",
                                                     effect_size_var = "es",
                                                     variance_var = "var_es")
}

# AQS-only
if (file.exists(file.path(data_dir, "data_aqs_only_es.csv"))) {
  data_aqs <- read.csv(file.path(data_dir, "data_aqs_only_es.csv"))
  data_aqs <- data_aqs %>%
    mutate(
      studynum = as.numeric(as.factor(study_group)),
      es = fisher_z,
      var_es = var_z
    )
  all_results$"aqs_only" <- run_moderator_analyses(data_aqs, "AQS-only",
                                                     effect_size_var = "es",
                                                     variance_var = "var_es")
}

# =============================================================================
# Export results
# =============================================================================

cat("\n\n=============================================================================\n")
cat("EXPORTING RESULTS\n")
cat("=============================================================================\n\n")

for (analysis_type in names(all_results)) {
  results_df <- all_results[[analysis_type]]
  
  if (nrow(results_df) > 0) {
    output_file <- file.path(output_dir, analysis_type, 
                            "moderator_analysis_9moderators.csv")
    dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
    write.csv(results_df, output_file, row.names = FALSE)
    cat(sprintf("  ✓ %s: %d moderators tested → %s\n", 
                analysis_type, nrow(results_df), basename(output_file)))
  }
}

# Create summary table across all analyses
cat("\n\nCreating cross-analysis summary...\n")

summary_all <- bind_rows(
  lapply(names(all_results), function(analysis_type) {
    if (nrow(all_results[[analysis_type]]) > 0) {
      all_results[[analysis_type]] %>%
        mutate(Analysis = analysis_type) %>%
        select(Analysis, everything())
    }
  })
)

if (nrow(summary_all) > 0) {
  summary_file <- file.path(output_dir, "moderator_analysis_summary.csv")
  write.csv(summary_all, summary_file, row.names = FALSE)
  cat(sprintf("  ✓ Summary: %s\n", summary_file))
}

cat("\n=============================================================================\n")
cat("✓ MODERATOR ANALYSIS COMPLETE\n")
cat("=============================================================================\n\n")
