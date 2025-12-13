# =============================================================================
# 08_MetaAnalysis_AQS_only.R
# =============================================================================
# Purpose: Meta-analysis of 2-way Secure/Insecure attachment continuity
#          from infancy to adulthood using Pearson's r (Fisher's z)
#          AQS measures only
# Input:   data_aqs_only_es.csv (with effect sizes)
# Output:  Note: Analysis not conducted due to insufficient samples
# =============================================================================

cat("=============================================================================\n")
cat("META-ANALYSIS: AQS-ONLY S/IS (Pearson's r)\n")
cat("=============================================================================\n\n")

# Load data
data <- read.csv("../Data/data_aqs_only_es.csv", stringsAsFactors = FALSE)

cat(sprintf("Number of samples: %d\n", nrow(data)))
cat(sprintf("Total participants: %d\n", sum(data$n, na.rm = TRUE)))

cat("\n")
cat("NOTE: Meta-analysis requires at least 3 independent studies.\n")
cat("With only %d samples, formal meta-analysis cannot be conducted.\n", nrow(data))
cat("\n")

if (nrow(data) > 0) {
  cat("Descriptive statistics for available samples:\n")
  cat(sprintf("  Mean r: %.3f\n", mean(data$r, na.rm = TRUE)))
  cat(sprintf("  Range: [%.3f, %.3f]\n", min(data$r, na.rm = TRUE), max(data$r, na.rm = TRUE)))
  cat("\n")
  
  cat("Individual study results:\n")
  for (i in 1:nrow(data)) {
    cat(sprintf("  %d. %s: r = %.3f [%.3f, %.3f], N = %d\n",
                i,
                data$author_year[i],
                data$r[i],
                tanh(data$fisher_z[i] - 1.96 * sqrt(data$var_z[i])),
                tanh(data$fisher_z[i] + 1.96 * sqrt(data$var_z[i])),
                data$n[i]))
  }
}

cat("\n")
cat("RECOMMENDATION: Combine AQS samples with SSP samples for sufficient power\n")
cat("(see Analysis 1: 2-way S/IS with both SSP and AQS measures)\n")
cat("=============================================================================\n")
