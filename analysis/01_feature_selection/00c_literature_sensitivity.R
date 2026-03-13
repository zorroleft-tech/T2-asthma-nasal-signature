# =========================================================================
# Script: 00c_literature_sensitivity.R
# Purpose: Test robustness to different publication count thresholds
# Author: Zuo
# Date: 2026-02-26
#
# Inputs:
#   - data/raw/01_all_genes_pubmed_counts.csv  # All genes with PubMed counts
#
# Outputs:
#   - data/derived/literature_sensitivity_stats.csv  # Sensitivity analysis results
# =========================================================================

library(dplyr)
library(readr)

# Set paths
# Get base directory (project root)
base_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# Verify we're running from project root
if (!all(dir.exists(file.path(base_dir, c("data", "analysis", "data_preparation"))))) {
  stop("Please run this script from the project root directory (the folder containing data/, analysis/, data_preparation/).", call. = FALSE)
}

# Define directories
DATA_RAW_DIR <- file.path(base_dir, "data", "raw")
DATA_DERIVED_DIR <- file.path(base_dir, "data", "derived")

# Create directory if not exist
dir.create(DATA_DERIVED_DIR, recursive = TRUE, showWarnings = FALSE)

# Log base directory and paths
cat(paste0("Base directory: ", base_dir, "\n"))
cat(paste0("Raw directory: ", DATA_RAW_DIR, "\n"))
cat(paste0("Derived directory: ", DATA_DERIVED_DIR, "\n\n"))

# Load all PubMed counts
all_counts <- read_csv(file.path(DATA_RAW_DIR, "01_all_genes_pubmed_counts.csv"))

# Test different thresholds
thresholds <- c(5, 10, 15, 20)

sensitivity_results <- data.frame(
  Threshold = thresholds,
  N_Candidates = sapply(thresholds, function(t) {
    sum(all_counts$PubMed_Count >= t, na.rm = TRUE)
  })
)

cat("\nSensitivity to publication threshold:\n")
print(sensitivity_results)

# Export sensitivity results
write_csv(
  sensitivity_results,
  file.path(DATA_DERIVED_DIR, "literature_sensitivity_stats.csv")
)

cat("\n✓ Sensitivity analysis complete\n")
cat(sprintf("→ Results exported: %s/literature_sensitivity_stats.csv\n\n", DATA_DERIVED_DIR))
