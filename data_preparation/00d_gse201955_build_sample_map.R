#!/usr/bin/env Rscript

# =========================================================================
# Script: 00d_gse201955_build_sample_map.R
# Purpose: Build and validate sample mapping for GSE201955
# Author: Zuo
# Date: 2026-02-26
#
# Inputs:
#   - data/raw/gse201955&GPL20301&GPL6791/GSE201955_RNAseq_118_processeddata.txt.gz  # Processed expression data
#   - data/raw/gse201955&GPL20301&GPL6791/GSE201955_family.soft.gz  # SOFT file with sample metadata
#   - data/raw/gse201955&GPL20301&GPL6791/GSE201955_series_matrix.txt.gz  # Series matrix file
#
# Outputs:
#   - data/processed_full/pheno_GSE201955_from_soft__RAW.rds  # Phenotype data from SOFT file
#   - data/processed_full/pheno_GSE201955_from_soft__RAW.csv  # Phenotype data from SOFT file (CSV)
#   - data/processed_full/expr_colname_to_gsm_mapping_GSE201955.csv  # Expression column to GSM mapping
#   - data/processed_full/gsm_comparison_GSE201955.csv  # GSM comparison between series matrix and SOFT
#   - output/logs/00d_gse201955_build_sample_map.log  # Log file
# =========================================================================

# Load required packages
library(tidyverse)
library(data.table)

# Set paths
# Get base directory (project root)
base_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# Verify we're running from project root
if (!all(dir.exists(file.path(base_dir, c("data", "analysis", "data_preparation"))))) {
  stop("Please run this script from the project root directory (the folder containing data/, analysis/, data_preparation/).", call. = FALSE)
}

# Define directories
raw_dir <- file.path(base_dir, "data", "raw")
processed_full_dir <- file.path(base_dir, "data", "processed_full")
logs_dir <- file.path(base_dir, "output", "logs")

# Create output directories
if (!dir.exists(processed_full_dir)) {
  dir.create(processed_full_dir, recursive = TRUE)
}
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

# Log file
log_file <- file.path(logs_dir, "00d_gse201955_build_sample_map.log")
sink(log_file, append = TRUE)

# Log base directory and paths
cat("\n=== Starting GSE201955 sample mapping build process ===\n")
cat(paste0("Execution time: ", Sys.time(), "\n"))
cat(paste0("Base directory: ", base_dir, "\n"))

# Define files
expr_file <- file.path(raw_dir, "gse201955&GPL20301&GPL6791", "GSE201955_RNAseq_118_processeddata.txt.gz")
soft_file <- file.path(raw_dir, "gse201955&GPL20301&GPL6791", "GSE201955_family.soft.gz")
series_matrix_file <- file.path(raw_dir, "gse201955&GPL20301&GPL6791", "GSE201955_series_matrix.txt.gz")
soft_pheno_output <- file.path(processed_full_dir, "pheno_GSE201955_from_soft__RAW.rds")
soft_pheno_csv_output <- file.path(processed_full_dir, "pheno_GSE201955_from_soft__RAW.csv")
mapping_output <- file.path(processed_full_dir, "expr_colname_to_gsm_mapping_GSE201955.csv")

# Log input and output files
cat(paste0("Input files: \n"))
cat(paste0("  - Expression file: ", expr_file, "\n"))
cat(paste0("  - SOFT file: ", soft_file, "\n"))
cat(paste0("  - Series matrix file: ", series_matrix_file, "\n"))
cat(paste0("Output files: \n"))
cat(paste0("  - Soft pheno RDS: ", soft_pheno_output, "\n"))
cat(paste0("  - Soft pheno CSV: ", soft_pheno_csv_output, "\n"))
cat(paste0("  - Mapping CSV: ", mapping_output, "\n"))

# Debug: Check file paths
cat(paste0("Debug: expr_file path: ", expr_file, "\n"))
cat(paste0("Debug: File exists: ", file.exists(expr_file), "\n"))

# Section 1: Parse SOFT file and extract sample metadata
cat("\n=== Section 1: Parsing SOFT file ===\n")
if (file.exists(soft_file)) {
  cat(paste0("Reading soft file: ", basename(soft_file), "\n"))
  
  # Read soft file line by line
  con <- gzfile(soft_file, "rt")
  lines <- readLines(con)
  close(con)
  
  # Initialize variables
  samples <- list()
  current_sample <- list()
  in_sample <- FALSE
  
  # Count total GSM lines first
  gsm_lines <- grep("^!Sample_geo_accession", lines)
  cat(paste0("Debug: Total GSM lines in SOFT file: ", length(gsm_lines), "\n"))
  
  # Parse soft file
  for (line in lines) {
    # Check for start of sample section
    if (grepl("^!Sample_geo_accession", line)) {
      # If we were in a sample, save it
      if (in_sample && length(current_sample) > 0) {
        samples <- c(samples, list(current_sample))
      }
      # Start new sample
      current_sample <- list()
      in_sample <- TRUE
      # Extract GSM
      gsm <- sub("^!Sample_geo_accession[[:space:]]*=?[[:space:]]*", "", line)
      current_sample$GSM <- gsm
    } else if (in_sample) {
      # Extract title
      if (grepl("^!Sample_title", line)) {
        title <- sub("^!Sample_title[[:space:]]+", "", line)
        current_sample$title <- title
      } 
      # Extract source_name_ch1
      else if (grepl("^!Sample_source_name_ch1", line)) {
        source_name_ch1 <- sub("^!Sample_source_name_ch1[[:space:]]+", "", line)
        current_sample$source_name_ch1 <- source_name_ch1
      }
      # Extract characteristics_ch1
      else if (grepl("^!Sample_characteristics_ch1", line)) {
        characteristic <- sub("^!Sample_characteristics_ch1[[:space:]]+", "", line)
        # If characteristics already exists, append to it
        if ("characteristics_ch1" %in% names(current_sample)) {
          current_sample$characteristics_ch1 <- c(current_sample$characteristics_ch1, characteristic)
        } else {
          current_sample$characteristics_ch1 <- characteristic
        }
      }
      # Check for end of sample section
      else if (grepl("^!", line) && !grepl("^!Sample_", line)) {
        # Save current sample
        if (length(current_sample) > 0) {
          samples <- c(samples, list(current_sample))
        }
        in_sample <- FALSE
      }
    }
  }
  
  # Save the last sample if any
  if (in_sample && length(current_sample) > 0) {
    samples <- c(samples, list(current_sample))
  }
  
  # Debug: Check sample count
  cat(paste0("Debug: Parsed ", length(samples), " samples from SOFT file\n"))
  
  cat(paste0("Found ", length(samples), " samples in soft file\n"))
  
  # Convert to data frame
  # First determine all possible column names
  all_cols <- unique(unlist(lapply(samples, names)))
  
  # Convert each sample to data frame with all columns
  soft_pheno <- do.call(rbind, lapply(samples, function(x) {
    # Handle list columns
    for (col in names(x)) {
      if (is.list(x[[col]])) {
        x[[col]] <- paste(unlist(x[[col]]), collapse = "; ")
      } else if (length(x[[col]]) > 1) {
        x[[col]] <- paste(x[[col]], collapse = "; ")
      }
    }
    # Create a data frame with all possible columns
    df <- as.data.frame(matrix(NA, nrow = 1, ncol = length(all_cols)), stringsAsFactors = FALSE)
    colnames(df) <- all_cols
    # Fill in existing values
    for (col in names(x)) {
      df[1, col] <- x[[col]]
    }
    df
  }))
  
  # Save soft pheno data
  saveRDS(soft_pheno, soft_pheno_output)
  write.csv(soft_pheno, soft_pheno_csv_output, row.names = FALSE)
  cat(paste0("Saved soft pheno data: ", basename(soft_pheno_output), "\n"))
  cat(paste0("Saved soft pheno CSV: ", basename(soft_pheno_csv_output), "\n"))
  
} else {
  cat(paste0("❌ Soft file not found: ", soft_file, "\n"))
  stop("Soft file not found")
}

# Section 2: Compare series_matrix (74) with soft (118) to confirm subset relationship
cat("\n=== Section 2: Comparing series_matrix and soft samples ===\n")
if (file.exists(series_matrix_file) && file.exists(soft_pheno_output)) {
  cat(paste0("Reading series_matrix file: ", basename(series_matrix_file), "\n"))
  
  # Read series_matrix file to extract GSM IDs
  con <- gzfile(series_matrix_file, "rt")
  lines <- readLines(con)
  close(con)
  
  # Extract GSM IDs from series_matrix
  gsm_lines <- grep("^!Sample_geo_accession", lines, value = TRUE)
  series_gsm <- sub("^!Sample_geo_accession[[:space:]]+", "", gsm_lines)
  
  cat(paste0("Found ", length(series_gsm), " samples in series_matrix\n"))
  
  # Read soft pheno data
  soft_pheno <- readRDS(soft_pheno_output)
  soft_gsm <- sub("[[:space:]]*=.*$", "", soft_pheno$GSM)
  soft_gsm <- gsub("[[:space:]]*$", "", soft_gsm)
  
  cat(paste0("Found ", length(soft_gsm), " samples in soft file\n"))
  
  # Check if series GSM is subset of soft GSM
  series_in_soft <- series_gsm %in% soft_gsm
  all_series_in_soft <- all(series_in_soft)
  
  cat(paste0("All series GSM in soft: ", ifelse(all_series_in_soft, "YES", "NO"), "\n"))
  
  if (!all_series_in_soft) {
    missing_series_gsm <- series_gsm[!series_in_soft]
    cat(paste0("Missing series GSM in soft: ", paste(missing_series_gsm, collapse = ", "), "\n"))
  }
  
  # Save comparison results
  comparison_df <- data.frame(
    GSM = series_gsm,
    in_soft = series_in_soft
  )
  comparison_output <- file.path(processed_full_dir, "gsm_comparison_GSE201955.csv")
  write.csv(comparison_df, comparison_output, row.names = FALSE)
  cat(paste0("Saved GSM comparison: ", basename(comparison_output), "\n"))
  
} else {
  if (!file.exists(series_matrix_file)) {
    cat(paste0("❌ Series matrix file not found: ", series_matrix_file, "\n"))
  }
  if (!file.exists(soft_pheno_output)) {
    cat(paste0("❌ Soft pheno file not found: ", soft_pheno_output, "\n"))
  }
}

# Section 3: Check if expression object carries column annotations
cat("\n=== Section 3: Checking column annotations in expression objects ===\n")

# Find RAW expression files
expr_files <- list.files(processed_full_dir, pattern = "expr_.*GSE201955.*__RAW.rds")
if (length(expr_files) > 0) {
  for (expr_file in expr_files) {
    expr_path <- file.path(processed_full_dir, expr_file)
    cat(paste0("Checking: ", basename(expr_path), "\n"))
    
    # Read the object
    expr_obj <- readRDS(expr_path)
    
    # Check for column annotations based on object type
    if (inherits(expr_obj, "ExpressionSet")) {
      # Check pData
      if (ncol(Biobase::pData(expr_obj)) > 0) {
        cat("  ✅ Found pData (column annotations)\n")
        cat(paste0("  Number of annotation columns: ", ncol(Biobase::pData(expr_obj)), "\n"))
        cat(paste0("  Annotation column names: ", paste(colnames(Biobase::pData(expr_obj)), collapse = ", "), "\n"))
      } else {
        cat("  ❌ No pData (column annotations) found\n")
      }
    } else if (inherits(expr_obj, "SummarizedExperiment")) {
      # Check colData
      if (ncol(SummarizedExperiment::colData(expr_obj)) > 0) {
        cat("  ✅ Found colData (column annotations)\n")
        cat(paste0("  Number of annotation columns: ", ncol(SummarizedExperiment::colData(expr_obj)), "\n"))
        cat(paste0("  Annotation column names: ", paste(colnames(SummarizedExperiment::colData(expr_obj)), collapse = ", "), "\n"))
      } else {
        cat("  ❌ No colData (column annotations) found\n")
      }
    } else if (is.matrix(expr_obj) || inherits(expr_obj, "dgCMatrix")) {
      # Check if matrix has column names
      if (!is.null(colnames(expr_obj))) {
        cat("  ✅ Found column names\n")
        cat(paste0("  First 5 column names: ", paste(head(colnames(expr_obj), 5), collapse = ", "), "\n"))
      } else {
        cat("  ❌ No column names found\n")
      }
      cat("  ❌ No additional column annotations found (matrix/dgCMatrix object)\n")
    } else {
      cat(paste0("  ⚠️  Unknown object type: ", class(expr_obj), "\n"))
    }
    
    # Memory cleanup
    rm(expr_obj)
    gc()
  }
} else {
  cat("No RAW expression files found for GSE201955\n")
}

# Section 4: Build expr_colname ↔ GSM mapping with integrity validation
cat("\n=== Section 4: Building expression column to GSM mapping ===\n")
if (file.exists(expr_file) && file.exists(soft_pheno_output)) {
  # Read expression file to get column names
  con <- gzfile(expr_file, "rt")
  # Read first line (header)
  header <- readLines(con, n = 1)
  close(con)
  
  # Split header to get column names
  col_names <- strsplit(header, "\t")[[1]]
  # Remove first column (gene ID)
  expr_col_names <- col_names[-1]
  
  cat(paste0("Found ", length(expr_col_names), " expression columns\n"))
  
  # Read soft pheno data
  soft_pheno <- readRDS(soft_pheno_output)
  
  cat(paste0("Found ", nrow(soft_pheno), " samples in soft pheno\n"))
  
  # Check if counts match
  if (length(expr_col_names) == nrow(soft_pheno)) {
    cat("✓ Number of expression columns matches number of soft samples\n")
    
    # Create mapping table
    mapping_df <- data.frame(
      expr_colname = expr_col_names,
      GSM = soft_pheno$GSM,
      title = soft_pheno$title,
      source_name_ch1 = soft_pheno$source_name_ch1,
      characteristics_ch1 = soft_pheno$characteristics_ch1,
      stringsAsFactors = FALSE
    )
    
    # Save mapping
    write.csv(mapping_df, mapping_output, row.names = FALSE)
    cat(paste0("Saved mapping table: ", basename(mapping_output), "\n"))
    
    # Show first few rows
    cat("\nFirst few rows of mapping table:\n")
    print(head(mapping_df))
    
    # Validate mapping
    cat("\nValidating mapping...\n")
    cat(paste0("Debug: nrow(mapping_df) = ", nrow(mapping_df), "\n"))
    cat(paste0("Debug: length(expr_col_names) = ", length(expr_col_names), "\n"))
    cat(paste0("Debug: nrow(soft_pheno) = ", nrow(soft_pheno), "\n"))
    if (nrow(mapping_df) != 118) {
      stop("Mapping should have exactly 118 rows")
    }
    if (anyDuplicated(mapping_df$expr_colname) != 0) {
      stop("Duplicate expression column names in mapping")
    }
    if (anyDuplicated(mapping_df$GSM) != 0) {
      stop("Duplicate GSM IDs in mapping")
    }
    cat("✓ Mapping validation passed\n")
    
  } else {
    cat("✗ Number of expression columns does not match number of soft samples\n")
    cat(paste0("Expression columns: ", length(expr_col_names), " vs Soft samples: ", nrow(soft_pheno), "\n"))
  }
  
} else {
  if (!file.exists(expr_file)) {
    cat(paste0("❌ Expression file not found: ", expr_file, "\n"))
  }
  if (!file.exists(soft_pheno_output)) {
    cat(paste0("❌ Soft pheno file not found: ", soft_pheno_output, "\n"))
  }
}

# Summary
cat("\n=== GSE201955 sample mapping build process completed ===\n")
# sink()

# Memory cleanup
rm(list = ls())
gc()
