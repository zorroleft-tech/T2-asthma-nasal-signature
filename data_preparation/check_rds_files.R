#!/usr/bin/env Rscript

# =========================================================================
# Script: check_rds_files.R
# Purpose: Check generated RDS files for quality and integrity (验收门禁)
# Author: Zuo
# Date: 2026-02-27
#
# Inputs:
#   - data/processed/phase3_outputs/*_pheno_raw.rds  # Standardized phenotype tables
#
# Outputs:
#   - Console output: Validation results
# =========================================================================

# Set base directory
base_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# Verify we're running from project root
required_dirs <- c('data', 'analysis', 'data_preparation')
required_paths <- file.path(base_dir, required_dirs)
dirs_exist <- all(dir.exists(required_paths))
if (!dirs_exist) {
  stop("Please run this script from the project root directory (the folder containing data/, analysis/, data_preparation/).", call. = FALSE)
}

# Define directories for our newly generated files
processed_dir <- file.path(base_dir, "data", "processed", "phase3_outputs")

# Get only *_pheno_raw.rds files
pheno_files <- list.files(processed_dir, pattern = "_pheno_raw\\.rds$", full.names = TRUE)

cat(paste0("Found ", length(pheno_files), " pheno_raw files\n"))

# Check each file
all_ok <- TRUE
for (pheno_file in pheno_files) {
  cat(paste0("\nChecking ", basename(pheno_file), "\n"))
  
  tryCatch({
    # Read the file
    data <- readRDS(pheno_file)
    
    # Check if data is a data frame
    if (!is.data.frame(data)) {
      cat("ERROR: File is not a data frame\n")
      all_ok <- FALSE
      next
    }
    
    # Check std_sample_id
    if ("std_sample_id" %in% colnames(data)) {
      if (anyDuplicated(data$std_sample_id) > 0) {
        cat(paste0("ERROR: Duplicate sample IDs in ", pheno_file, "\n"))
        all_ok <- FALSE
        next
      }
      if (sum(is.na(data$std_sample_id)) > 0) {
        cat(paste0("ERROR: NA sample IDs in ", pheno_file, "\n"))
        all_ok <- FALSE
        next
      }
    } else {
      cat("ERROR: std_sample_id column not found\n")
      all_ok <- FALSE
      next
    }
    
    # Check std_group_label
    if ("std_group_label" %in% colnames(data)) {
      if (all(is.na(data$std_group_label))) {
        cat(paste0("ERROR: All group labels are NA in ", pheno_file, "\n"))
        all_ok <- FALSE
        next
      }
      
      # Show summary
      cat(paste0("Samples: ", nrow(data), "\n"))
      cat(paste0("Group distribution: ", paste(names(table(data$std_group_label)), collapse = ", "), "\n"))
      cat("OK\n")
    } else {
      cat("ERROR: std_group_label column not found\n")
      all_ok <- FALSE
      next
    }
  }, error = function(e) {
    cat(paste0("ERROR: ", e$message, "\n"))
    all_ok <- FALSE
  })
}

# Print final result
cat("\n=== Final Result ===\n")
if (all_ok) {
  cat("All pheno files passed QC!\n")
} else {
  cat("Some pheno files failed QC!\n")
  quit(status = 1)
}
