#!/usr/bin/env Rscript

# Script 00g: Build mapping between expression column names and GSM IDs
# Function: Perform three checks to build expr_colname ↔ GSM mapping without re-downloading data

# Load required packages
library(tidyverse)
library(data.table)

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
raw_dir <- file.path(base_dir, "data", "raw")
processed_full_dir <- file.path(base_dir, "data", "processed_full")
logs_dir <- file.path(base_dir, "output", "logs")

# Create output directories if not exist
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

# Log file
log_file <- file.path(logs_dir, "00g_build_mapping.log")
sink(log_file, append = TRUE)

cat("\n=== Starting mapping build process ===\n")
cat(paste0("Execution time: ", Sys.time(), "\n"))

# Define files
expr_file <- file.path(raw_dir, "gse201955&GPL20301&GPL6791", "GSE201955_RNAseq_118_processeddata.txt.gz")
soft_file <- file.path(raw_dir, "gse201955&GPL20301&GPL6791", "GSE201955_family.soft.gz")
soft_pheno_file <- file.path(processed_full_dir, "pheno_GSE201955_from_soft__RAW.rds")

# Check 1: Examine expression file header
cat("\n=== Check 1: Examine expression file header ===\n")
if (file.exists(expr_file)) {
  cat(paste0("Reading expression file header: ", basename(expr_file), "\n"))
  
  # Read first 10 lines
  con <- gzfile(expr_file, "rt")
  header_lines <- readLines(con, n = 10)
  close(con)
  
  cat("First 10 lines of expression file:\n")
  for (i in 1:length(header_lines)) {
    cat(paste0(i, ": ", header_lines[i], "\n"))
  }
  
  # Check for comment lines
  comment_lines <- grep("^#", header_lines)
  if (length(comment_lines) > 0) {
    cat(paste0("Found ", length(comment_lines), " comment lines\n"))
  }
  
  # Check for GSM in header
  gsm_in_header <- any(grepl("GSM", header_lines, ignore.case = TRUE))
  cat(paste0("GSM found in header: ", ifelse(gsm_in_header, "YES", "NO"), "\n"))
  
  # Check for sample titles
  if (length(header_lines) > 1) {
    # Assume second line might contain sample titles
    cat("Checking for sample titles in header...\n")
  }
  
} else {
  cat(paste0("❌ Expression file not found: ", expr_file, "\n"))
}

# Check 2: Extract GSM info from soft file
cat("\n=== Check 2: Extract GSM info from soft file ===\n")
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
    gsm <- sub("^!Sample_geo_accession[[:space:]]+", "", line)
    current_sample$GSM <- gsm
  } else if (in_sample) {
    # Extract title
    if (grepl("^!Sample_title", line)) {
      title <- sub("^!Sample_title[[:space:]]+", "", line)
      current_sample$title <- title
    } 
    # Extract supplementary_file
    else if (grepl("^!Sample_supplementary_file", line)) {
      supplementary_file <- sub("^!Sample_supplementary_file[[:space:]]+", "", line)
      current_sample$supplementary_file <- supplementary_file
    }
    # Extract relation
    else if (grepl("^!Sample_relation", line)) {
      relation <- sub("^!Sample_relation[[:space:]]+", "", line)
      current_sample$relation <- relation
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

cat(paste0("Found ", length(samples), " samples in soft file\n"))

# Check for X46b4_4d67-like tokens in soft file
cat("Checking for X46b4_4d67-like tokens in soft file...\n")
found_tokens <- FALSE
for (sample in samples) {
  # Check all fields for X****-like patterns
  for (field in names(sample)) {
    value <- sample[[field]]
    if (is.character(value)) {
      if (any(grepl("X[0-9a-f]{4}_[0-9a-f]{4}", value, ignore.case = TRUE))) {
        matching_values <- value[grepl("X[0-9a-f]{4}_[0-9a-f]{4}", value, ignore.case = TRUE)]
        for (val in matching_values) {
          cat(paste0("Found token in GSM ", sample$GSM, " (", field, "): ", val, "\n"))
          found_tokens <- TRUE
        }
      }
    } else if (is.list(value)) {
      for (item in value) {
        if (is.character(item) && grepl("X[0-9a-f]{4}_[0-9a-f]{4}", item, ignore.case = TRUE)) {
          cat(paste0("Found token in GSM ", sample$GSM, " (", field, "): ", item, "\n"))
          found_tokens <- TRUE
        }
      }
    }
  }
}

if (!found_tokens) {
  cat("No X46b4_4d67-like tokens found in soft file\n")
}

} else {
  cat(paste0("❌ Soft file not found: ", soft_file, "\n"))
}

# Check 3: Match expression column names with soft samples
cat("\n=== Check 3: Match expression column names with soft samples ===\n")
if (file.exists(expr_file) && file.exists(soft_pheno_file)) {
  cat("Reading expression file column names...\n")
  
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
  cat("First 10 expression column names:\n")
  print(head(expr_col_names, 10))
  
  # Read soft pheno data
  cat("Reading soft pheno data...\n")
  soft_pheno <- readRDS(soft_pheno_file)
  
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
    mapping_output <- file.path(processed_full_dir, "expr_colname_to_gsm_mapping_GSE201955.csv")
    write.csv(mapping_df, mapping_output, row.names = FALSE)
    cat(paste0("Saved mapping table: ", basename(mapping_output), "\n"))
    
    # Show first few rows
    cat("\nFirst few rows of mapping table:\n")
    print(head(mapping_df))
    
  } else {
    cat("✗ Number of expression columns does not match number of soft samples\n")
    cat(paste0("Expression columns: ", length(expr_col_names), " vs Soft samples: ", nrow(soft_pheno), "\n"))
  }
  
} else {
  if (!file.exists(expr_file)) {
    cat(paste0("❌ Expression file not found: ", expr_file, "\n"))
  }
  if (!file.exists(soft_pheno_file)) {
    cat(paste0("❌ Soft pheno file not found: ", soft_pheno_file, "\n"))
  }
}

# Summary
cat("\n=== Mapping build process completed ===\n")
sink()

# Memory cleanup
rm(list = ls())
gc()
