# =============================================================================
# effect_size_utils.r
# =============================================================================
# Purpose: Calculate effect sizes (Pearson's r and Cohen's kappa) from 
#          contingency tables for meta-analysis
# =============================================================================

# =============================================================================
# FUNCTION: calculate_r_from_2x2
# =============================================================================
# Calculate Pearson's r (phi coefficient) from a 2x2 contingency table
# Also returns Fisher's z transformation and variance
#
# Parameters:
#   a, b, c, d: Cell frequencies in 2x2 table
#               Layout: [a, b]
#                       [c, d]
#
# Returns:
#   List with: r (Pearson's r), z (Fisher's z), var_z (variance of z), 
#              se_z (standard error of z), n (sample size)
# =============================================================================

calculate_r_from_2x2 <- function(a, b, c, d) {
  # Total sample size
  n <- a + b + c + d
  
  # Check for valid input
  if (n == 0) {
    return(list(r = NA, z = NA, var_z = NA, se_z = NA, n = 0))
  }
  
  # Calculate phi coefficient (equivalent to Pearson's r for 2x2 table)
  numerator <- (a * d) - (b * c)
  denominator <- sqrt((a + b) * (c + d) * (a + c) * (b + d))
  
  # Handle zero denominators
  if (denominator == 0) {
    r <- 0
  } else {
    r <- numerator / denominator
  }
  
  # Constrain r to [-1, 1] to avoid Fisher's z issues
  r <- max(-0.9999, min(0.9999, r))
  
  # Fisher's z transformation
  z <- 0.5 * log((1 + r) / (1 - r))
  
  # Variance of Fisher's z
  var_z <- 1 / (n - 3)
  
  # Standard error of Fisher's z
  se_z <- sqrt(var_z)
  
  return(list(
    r = r,
    z = z,
    var_z = var_z,
    se_z = se_z,
    n = n
  ))
}


# =============================================================================
# FUNCTION: calculate_kappa
# =============================================================================
# Calculate Cohen's kappa from a k x k contingency table
#
# Parameters:
#   contingency_matrix: A k x k matrix of observed frequencies
#
# Returns:
#   List with: kappa (Cohen's kappa), var_kappa (variance), 
#              se_kappa (standard error), n (sample size),
#              po (observed agreement), pe (expected agreement)
# =============================================================================

calculate_kappa <- function(contingency_matrix) {
  # Total sample size
  n <- sum(contingency_matrix)
  
  # Check for valid input
  if (n == 0) {
    return(list(kappa = NA, var_kappa = NA, se_kappa = NA, 
                n = 0, po = NA, pe = NA))
  }
  
  # Observed agreement (proportion of diagonal)
  po <- sum(diag(contingency_matrix)) / n
  
  # Expected agreement by chance
  row_sums <- rowSums(contingency_matrix)
  col_sums <- colSums(contingency_matrix)
  pe <- sum((row_sums * col_sums) / n^2)
  
  # Cohen's kappa
  if (pe >= 1) {
    kappa <- 0
  } else {
    kappa <- (po - pe) / (1 - pe)
  }
  
  # Variance of kappa (Fleiss et al., 1969; simplified approximation)
  # More accurate formula from Fleiss (1981)
  if (pe >= 1) {
    var_kappa <- NA
  } else {
    var_kappa <- (po * (1 - po)) / (n * (1 - pe)^2)
  }
  
  # Standard error
  se_kappa <- sqrt(var_kappa)
  
  return(list(
    kappa = kappa,
    var_kappa = var_kappa,
    se_kappa = se_kappa,
    n = n,
    po = po,
    pe = pe
  ))
}


# =============================================================================
# FUNCTION: calculate_effect_sizes_2way
# =============================================================================
# Calculate Pearson's r for all rows in a dataset with 2x2 contingency data
#
# Parameters:
#   data: Data frame with contingency table data
#   cell_cols: Vector of 4 column names for cells [a, b, c, d]
#
# Returns:
#   Data frame with original data plus effect size columns
# =============================================================================

calculate_effect_sizes_2way <- function(data, cell_cols = c("cell_11", "cell_12", "cell_21", "cell_22")) {
  # Initialize result columns
  data$r <- NA
  data$fisher_z <- NA
  data$var_z <- NA
  data$se_z <- NA
  
  # Calculate effect sizes for each row
  for (i in 1:nrow(data)) {
    # Extract cell values
    cells <- as.numeric(data[i, cell_cols])
    
    # Check if all cells are non-missing
    if (all(!is.na(cells)) && sum(cells) > 0) {
      # Calculate effect size
      es <- calculate_r_from_2x2(cells[1], cells[2], cells[3], cells[4])
      
      # Store results
      data$r[i] <- es$r
      data$fisher_z[i] <- es$z
      data$var_z[i] <- es$var_z
      data$se_z[i] <- es$se_z
    }
  }
  
  return(data)
}


# =============================================================================
# FUNCTION: calculate_effect_sizes_3way
# =============================================================================
# Calculate Cohen's kappa for all rows in a dataset with 3x3 contingency data
#
# Parameters:
#   data: Data frame with contingency table data
#   cell_cols: Vector of 9 column names for 3x3 matrix (row-major order)
#              [A-Ds, A-F, A-E, B-Ds, B-F, B-E, C-Ds, C-F, C-E]
#
# Returns:
#   Data frame with original data plus effect size columns
# =============================================================================

calculate_effect_sizes_3way <- function(data, 
                                        cell_cols = c("cell_A_Ds", "cell_A_F", "cell_A_E",
                                                      "cell_B_Ds", "cell_B_F", "cell_B_E",
                                                      "cell_C_Ds", "cell_C_F", "cell_C_E")) {
  # Initialize result columns
  data$kappa <- NA
  data$var_kappa <- NA
  data$se_kappa <- NA
  data$po <- NA
  data$pe <- NA
  
  # Calculate effect sizes for each row
  for (i in 1:nrow(data)) {
    # Extract cell values
    cells <- as.numeric(data[i, cell_cols])
    
    # Check if all cells are non-missing
    if (all(!is.na(cells)) && sum(cells) > 0) {
      # Build 3x3 matrix (row-major order)
      contingency_matrix <- matrix(cells, nrow = 3, ncol = 3, byrow = TRUE)
      
      # Calculate effect size
      es <- calculate_kappa(contingency_matrix)
      
      # Store results
      data$kappa[i] <- es$kappa
      data$var_kappa[i] <- es$var_kappa
      data$se_kappa[i] <- es$se_kappa
      data$po[i] <- es$po
      data$pe[i] <- es$pe
    }
  }
  
  return(data)
}


# =============================================================================
# FUNCTION: calculate_effect_sizes_4way
# =============================================================================
# Calculate Cohen's kappa for all rows in a dataset with 4x4 contingency data
#
# Parameters:
#   data: Data frame with contingency table data
#   cell_cols: Vector of 16 column names for 4x4 matrix (row-major order)
#              [A-Ds, A-F, A-E, A-U, B-Ds, B-F, B-E, B-U, 
#               C-Ds, C-F, C-E, C-U, D-Ds, D-F, D-E, D-U]
#
# Returns:
#   Data frame with original data plus effect size columns
# =============================================================================

calculate_effect_sizes_4way <- function(data,
                                        cell_cols = c("cell_A_Ds", "cell_A_F", "cell_A_E", "cell_A_U",
                                                      "cell_B_Ds", "cell_B_F", "cell_B_E", "cell_B_U",
                                                      "cell_C_Ds", "cell_C_F", "cell_C_E", "cell_C_U",
                                                      "cell_D_Ds", "cell_D_F", "cell_D_E", "cell_D_U")) {
  # Initialize result columns
  data$kappa <- NA
  data$var_kappa <- NA
  data$se_kappa <- NA
  data$po <- NA
  data$pe <- NA
  
  # Calculate effect sizes for each row
  for (i in 1:nrow(data)) {
    # Extract cell values
    cells <- as.numeric(data[i, cell_cols])
    
    # Check if all cells are non-missing
    if (all(!is.na(cells)) && sum(cells) > 0) {
      # Build 4x4 matrix (row-major order)
      contingency_matrix <- matrix(cells, nrow = 4, ncol = 4, byrow = TRUE)
      
      # Calculate effect size
      es <- calculate_kappa(contingency_matrix)
      
      # Store results
      data$kappa[i] <- es$kappa
      data$var_kappa[i] <- es$var_kappa
      data$se_kappa[i] <- es$se_kappa
      data$po[i] <- es$po
      data$pe[i] <- es$pe
    }
  }
  
  return(data)
}


# =============================================================================
# UTILITY: print_effect_size_summary
# =============================================================================
# Print summary statistics for calculated effect sizes
#
# Parameters:
#   data: Data frame with effect size columns
#   es_type: "r" for Pearson's r or "kappa" for Cohen's kappa
# =============================================================================

print_effect_size_summary <- function(data, es_type = "r") {
  if (es_type == "r") {
    cat("Pearson's r Summary:\n")
    cat(sprintf("  Mean r: %.3f\n", mean(data$r, na.rm = TRUE)))
    cat(sprintf("  SD r: %.3f\n", sd(data$r, na.rm = TRUE)))
    cat(sprintf("  Range: [%.3f, %.3f]\n", min(data$r, na.rm = TRUE), max(data$r, na.rm = TRUE)))
    cat(sprintf("  Mean Fisher's z: %.3f\n", mean(data$fisher_z, na.rm = TRUE)))
    cat(sprintf("  Mean SE(z): %.3f\n", mean(data$se_z, na.rm = TRUE)))
  } else if (es_type == "kappa") {
    cat("Cohen's kappa Summary:\n")
    cat(sprintf("  Mean kappa: %.3f\n", mean(data$kappa, na.rm = TRUE)))
    cat(sprintf("  SD kappa: %.3f\n", sd(data$kappa, na.rm = TRUE)))
    cat(sprintf("  Range: [%.3f, %.3f]\n", min(data$kappa, na.rm = TRUE), max(data$kappa, na.rm = TRUE)))
    cat(sprintf("  Mean SE(kappa): %.3f\n", mean(data$se_kappa, na.rm = TRUE)))
    cat(sprintf("  Mean observed agreement: %.3f\n", mean(data$po, na.rm = TRUE)))
  }
  cat(sprintf("  Sample sizes: %d - %d\n", min(data$n, na.rm = TRUE), max(data$n, na.rm = TRUE)))
}
