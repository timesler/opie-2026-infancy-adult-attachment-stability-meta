# Load required packages
library(readxl)
library(dplyr)
library(tidyr)

# Set working directory to script location
if (interactive()) {
  script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
  setwd(script_dir)
}

# Define paths
data_dir <- "../../data"
output_dir <- "../Data"

cat("=============================================================================\n")
cat("DATA PREPARATION FOR INFANT-ADULT ATTACHMENT META-ANALYSIS\n")
cat("=============================================================================\n\n")

# =============================================================================
# Load cleaned data extraction table
# =============================================================================

cat("Loading cleaned data extraction table...\n")
df <- read.csv(file.path(data_dir, "data_extraction_table.csv"), stringsAsFactors = FALSE)

cat(sprintf("  - Loaded %d rows of data\n", nrow(df)))
cat(sprintf("  - %d columns\n", ncol(df)))

# =============================================================================
# Basic data cleaning
# =============================================================================

cat("\nBasic data cleaning:\n")

# Convert numeric columns
df <- df %>%
  mutate(
    year = as.numeric(year),
    n = as.numeric(n),
    infant_age_months = as.numeric(infant_age_months),
    adult_age_years = as.numeric(adult_age_years)
  )

cat("  - Data types converted\n")

# =============================================================================
# Create analysis-specific datasets
# =============================================================================

cat("\nCreating analysis-specific datasets:\n")

# 2-way S/IS analysis
df_2way_sis <- df %>%
  filter(!is.na(sis_s_to_s) | !is.na(sis_s_to_is) | !is.na(sis_is_to_s) | !is.na(sis_is_to_is) | !is.na(sis_pearsons_r_reported)) %>%
  mutate(
    cell_11 = sis_s_to_s,
    cell_12 = sis_s_to_is,
    cell_21 = sis_is_to_s,
    cell_22 = sis_is_to_is,
    cell_sum = cell_11 + cell_12 + cell_21 + cell_22,
    valid = ifelse(cell_sum > 0, TRUE, FALSE),
    r = sis_pearsons_r_reported
  ) %>%
  select(author, subsample, year, country, source, prior_meta_inclusion, n, infant_age_months, adult_age_years, 
         measures_used, parent_risk, child_risk, sample_type, 
         cell_11, cell_12, cell_21, cell_22, cell_sum, valid, r)

cat(sprintf("  - 2-way S/IS: %d rows\n", nrow(df_2way_sis)))

# 2-way O/D analysis
df_2way_od <- df %>%
  filter(!is.na(odo_o_to_o) | !is.na(odo_o_to_d) | !is.na(odo_d_to_o) | !is.na(odo_d_to_d) | !is.na(odo_pearsons_r_reported)) %>%
  mutate(
    cell_11 = odo_o_to_o,
    cell_12 = odo_o_to_d,
    cell_21 = odo_d_to_o,
    cell_22 = odo_d_to_d,
    cell_sum = cell_11 + cell_12 + cell_21 + cell_22,
    valid = ifelse(cell_sum > 0, TRUE, FALSE),
    r = odo_pearsons_r_reported
  ) %>%
  select(author, subsample, year, country, source, prior_meta_inclusion, n, infant_age_months, adult_age_years, 
         measures_used, parent_risk, child_risk, sample_type, 
         cell_11, cell_12, cell_21, cell_22, cell_sum, valid, r)

cat(sprintf("  - 2-way O/D: %d rows\n", nrow(df_2way_od)))

# 3-way analysis - explicit column selection
df_3way <- df %>%
  filter(!is.na(X3way_a_to_a) | !is.na(X3way_a_to_b) | !is.na(X3way_a_to_c) |
         !is.na(X3way_b_to_a) | !is.na(X3way_b_to_b) | !is.na(X3way_b_to_c) |
         !is.na(X3way_c_to_a) | !is.na(X3way_c_to_b) | !is.na(X3way_c_to_c)) %>%
  mutate(
    cell_11 = X3way_a_to_a,
    cell_12 = X3way_a_to_b,
    cell_13 = X3way_a_to_c,
    cell_21 = X3way_b_to_a,
    cell_22 = X3way_b_to_b,
    cell_23 = X3way_b_to_c,
    cell_31 = X3way_c_to_a,
    cell_32 = X3way_c_to_b,
    cell_33 = X3way_c_to_c,
    cell_sum = cell_11 + cell_12 + cell_13 + cell_21 + cell_22 + cell_23 + cell_31 + cell_32 + cell_33,
    valid = ifelse(cell_sum > 0, TRUE, FALSE)
  ) %>%
  select(author, subsample, year, country, source, prior_meta_inclusion, n, infant_age_months, adult_age_years, 
         measures_used, parent_risk, child_risk, sample_type, 
         cell_11, cell_12, cell_13, cell_21, cell_22, cell_23, cell_31, cell_32, cell_33,
         cell_sum, valid)

cat(sprintf("  - 3-way: %d rows\n", nrow(df_3way)))

# 4-way analysis - explicit column selection
df_4way <- df %>%
  filter(!is.na(X4way_a_to_a) | !is.na(X4way_a_to_b) | !is.na(X4way_a_to_c) | !is.na(X4way_a_to_d) |
         !is.na(X4way_b_to_a) | !is.na(X4way_b_to_b) | !is.na(X4way_b_to_c) | !is.na(X4way_b_to_d) |
         !is.na(X4way_c_to_a) | !is.na(X4way_c_to_b) | !is.na(X4way_c_to_c) | !is.na(X4way_c_to_d) |
         !is.na(X4way_d_to_a) | !is.na(X4way_d_to_b) | !is.na(X4way_d_to_c) | !is.na(X4way_d_to_d)) %>%
  mutate(
    cell_11 = X4way_a_to_a, cell_12 = X4way_a_to_b, cell_13 = X4way_a_to_c, cell_14 = X4way_a_to_d,
    cell_21 = X4way_b_to_a, cell_22 = X4way_b_to_b, cell_23 = X4way_b_to_c, cell_24 = X4way_b_to_d,
    cell_31 = X4way_c_to_a, cell_32 = X4way_c_to_b, cell_33 = X4way_c_to_c, cell_34 = X4way_c_to_d,
    cell_41 = X4way_d_to_a, cell_42 = X4way_d_to_b, cell_43 = X4way_d_to_c, cell_44 = X4way_d_to_d,
    cell_sum = cell_11 + cell_12 + cell_13 + cell_14 + cell_21 + cell_22 + cell_23 + cell_24 +
               cell_31 + cell_32 + cell_33 + cell_34 + cell_41 + cell_42 + cell_43 + cell_44,
    valid = ifelse(cell_sum > 0, TRUE, FALSE)
  ) %>%
  select(author, subsample, year, country, source, prior_meta_inclusion, n, infant_age_months, adult_age_years, 
         measures_used, parent_risk, child_risk, sample_type, 
         cell_11, cell_12, cell_13, cell_14, cell_21, cell_22, cell_23, cell_24,
         cell_31, cell_32, cell_33, cell_34, cell_41, cell_42, cell_43, cell_44,
         cell_sum, valid)

cat(sprintf("  - 4-way: %d rows\n", nrow(df_4way)))

# SSP-only analysis (only SSP measures)
df_ssp_only <- df_2way_sis %>%
  filter(grepl("SSP", measures_used, ignore.case = TRUE))

cat(sprintf("  - SSP-only: %d rows\n", nrow(df_ssp_only)))

# AQS-only analysis (only AQS measures)
df_aqs_only <- df_2way_sis %>%
  filter(grepl("AQS", measures_used, ignore.case = TRUE))

cat(sprintf("  - AQS-only: %d rows\n", nrow(df_aqs_only)))

# =============================================================================
# Save analysis datasets
# =============================================================================

cat("\nSaving analysis datasets:\n")

write.csv(df_2way_sis, file.path(output_dir, "data_2way_sis.csv"), row.names = FALSE)
write.csv(df_2way_od, file.path(output_dir, "data_2way_od.csv"), row.names = FALSE)
write.csv(df_3way, file.path(output_dir, "data_3way.csv"), row.names = FALSE)
write.csv(df_4way, file.path(output_dir, "data_4way.csv"), row.names = FALSE)
write.csv(df_ssp_only, file.path(output_dir, "data_ssp_only.csv"), row.names = FALSE)
write.csv(df_aqs_only, file.path(output_dir, "data_aqs_only.csv"), row.names = FALSE)

cat("  - Saved all analysis datasets\n")

# =============================================================================
# Summary statistics
# =============================================================================

cat("\n=============================================================================\n")
cat("DATA PREPARATION SUMMARY\n")
cat("=============================================================================\n\n")

cat("Analysis datasets created:\n")
cat(sprintf("  - 2-way S/IS: %d rows from %d unique studies\n", 
            nrow(df_2way_sis), n_distinct(df_2way_sis$author)))
cat(sprintf("  - 2-way O/D: %d rows from %d unique studies\n", 
            nrow(df_2way_od), n_distinct(df_2way_od$author)))
cat(sprintf("  - 3-way: %d rows from %d unique studies\n", 
            nrow(df_3way), n_distinct(df_3way$author)))
cat(sprintf("  - 4-way: %d rows from %d unique studies\n", 
            nrow(df_4way), n_distinct(df_4way$author)))
cat(sprintf("  - SSP-only: %d rows from %d unique studies\n", 
            nrow(df_ssp_only), n_distinct(df_ssp_only$author)))
cat(sprintf("  - AQS-only: %d rows from %d unique studies\n", 
            nrow(df_aqs_only), n_distinct(df_aqs_only$author)))

cat("\n✓ Data preparation complete!\n\n")
