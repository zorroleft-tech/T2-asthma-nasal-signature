#!/usr/bin/env Rscript

# Script 00e: Compare GSM lists from soft file and series_matrix
# Function: Load both pheno files and compare GSM lists

# Load required packages
library(tidyverse)

# Set paths
# Check if we're in the project root or data_preparation directory
current_dir <- getwd()
if (basename(current_dir) == "data_preparation") {
  # Running from data_preparation directory
  base_dir <- dirname(current_dir)
} else {
  # Running from project root directory
  base_dir <- current_dir
}
processed_full_dir <- file.path(base_dir, "data", "processed_full")
logs_dir <- file.path(base_dir, "output", "logs")

# Create output directories if not exist
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

# Log file
log_file <- file.path(logs_dir, "00e_compare_gsm.log")
sink(log_file, append = TRUE)

cat("\n=== Starting GSM comparison ===\n")
cat(paste0("Execution time: ", Sys.time(), "\n"))

# Load pheno files
soft_pheno_file <- file.path(processed_full_dir, "pheno_GSE201955_from_soft__RAW.rds")
series_pheno_file <- file.path(processed_full_dir, "pheno_GSE201955__FULL.rds")

if (!file.exists(soft_pheno_file)) {
  cat(paste0("❌ Soft pheno file not found: ", soft_pheno_file, "\n"))
  sink()
  stop("Soft pheno file not found")
}

if (!file.exists(series_pheno_file)) {
  cat(paste0("❌ Series pheno file not found: ", series_pheno_file, "\n"))
  sink()
  stop("Series pheno file not found")
}

# Read pheno files
cat("Reading soft pheno file...\n")
soft_pheno <- readRDS(soft_pheno_file)

cat("Reading series pheno file...\n")
series_pheno <- readRDS(series_pheno_file)

# Extract GSM lists
soft_gsm <- soft_pheno$GSM
series_gsm <- series_pheno$sample_id  # Assuming sample_id is the column name for GSM

# Clean GSM values (remove leading/trailing whitespace, equals signs, and quotes)
soft_gsm <- gsub("^[[:space:]]*=[[:space:]]*", "", soft_gsm)
soft_gsm <- gsub("[[:space:]]*=[[:space:]]*$", "", soft_gsm)
soft_gsm <- gsub("^\"|\"$", "", soft_gsm)
soft_gsm <- trimws(soft_gsm)

series_gsm <- gsub("^[[:space:]]*=[[:space:]]*", "", series_gsm)
series_gsm <- gsub("[[:space:]]*=[[:space:]]*$", "", series_gsm)
series_gsm <- gsub("^\"|\"$", "", series_gsm)
series_gsm <- trimws(series_gsm)

# Summary
cat("\n=== GSM Comparison Summary ===\n")
cat(paste0("Total GSM from soft file: ", length(soft_gsm), "\n"))
cat(paste0("Total GSM from series_matrix: ", length(series_gsm), "\n"))

# Check if series GSM is subset of soft GSM
is_subset <- all(series_gsm %in% soft_gsm)
cat(paste0("Is series GSM subset of soft GSM: ", ifelse(is_subset, "YES", "NO"), "\n"))

# Find missing GSMs if any
if (!is_subset) {
  missing_gsm <- series_gsm[!series_gsm %in% soft_gsm]
  cat(paste0("Missing GSMs: ", length(missing_gsm), "\n"))
  if (length(missing_gsm) > 0) {
    cat("Missing GSMs list: ", paste(missing_gsm, collapse = ", "), "\n")
  }
}

# Show first few GSMs from both sources
cat("\nFirst 10 GSMs from soft file:\n")
print(head(soft_gsm, 10))

cat("\nFirst 10 GSMs from series_matrix:\n")
print(head(series_gsm, 10))

# Detailed comparison
cat("\n=== Detailed Comparison ===\n")
cat(paste0("Number of GSMs in soft but not in series: ", sum(!soft_gsm %in% series_gsm), "\n"))
cat(paste0("Number of GSMs in series but not in soft: ", sum(!series_gsm %in% soft_gsm), "\n"))

# Save comparison results
comparison_df <- data.frame(
  GSM = unique(c(soft_gsm, series_gsm)),
  in_soft = unique(c(soft_gsm, series_gsm)) %in% soft_gsm,
  in_series = unique(c(soft_gsm, series_gsm)) %in% series_gsm,
  stringsAsFactors = FALSE
)

comparison_output <- file.path(processed_full_dir, "gsm_comparison_GSE201955.csv")
write.csv(comparison_df, comparison_output, row.names = FALSE)
cat(paste0("\nSaved comparison results: ", basename(comparison_output), "\n"))

cat("\n=== GSM comparison completed ===\n")
sink()

# Memory cleanup
rm(list = ls())
gc()
