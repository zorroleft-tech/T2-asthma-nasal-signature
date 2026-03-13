#!/usr/bin/env Rscript

# Script 00f: Check for column annotations in RAW files
# Function: Check if RAW.rds files contain column annotations (pData/colData/sample_table/targets)

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
processed_dir <- file.path(base_dir, "data", "processed")
logs_dir <- file.path(base_dir, "output", "logs")

# Create output directories if not exist
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

# Log file
log_file <- file.path(logs_dir, "00f_check_col_annotations.log")
sink(log_file, append = TRUE)

cat("\n=== Starting column annotations check ===\n")
cat(paste0("Execution time: ", Sys.time(), "\n"))

# Define GSE201955 RAW files
raw_files <- c(
  file.path(processed_full_dir, "expr_gse201955&GPL20301&GPL6791__RAW.rds"),
  file.path(processed_dir, "expr_gse201955&GPL20301&GPL6791__RAW.rds")
)

# Check each file
for (file_path in raw_files) {
  if (file.exists(file_path)) {
    cat(paste0("\nChecking file: ", basename(file_path), "\n"))
    
    # Read the file
    obj <- readRDS(file_path)
    
    # Check object type
    cat(paste0("Object class: ", class(obj)[1], "\n"))
    
    # Check for different types of column annotations
    has_col_data <- FALSE
    
    # Check if it's an ExpressionSet
    if (class(obj)[1] == "ExpressionSet") {
      cat("✓ Object is an ExpressionSet\n")
      # Check for pData
      if ("phenoData" %in% slotNames(obj)) {
        pdata <- Biobase::pData(obj)
        cat(paste0("✓ Found pData with ", nrow(pdata), " samples and ", ncol(pdata), " columns\n"))
        cat("  Column names: ", paste(colnames(pdata), collapse = ", "), "\n")
        has_col_data <- TRUE
      }
    }
    
    # Check if it's a SummarizedExperiment
    else if (class(obj)[1] == "SummarizedExperiment") {
      cat("✓ Object is a SummarizedExperiment\n")
      # Check for colData
      if ("colData" %in% slotNames(obj)) {
        coldata <- SummarizedExperiment::colData(obj)
        cat(paste0("✓ Found colData with ", nrow(coldata), " samples and ", ncol(coldata), " columns\n"))
        cat("  Column names: ", paste(colnames(coldata), collapse = ", "), "\n")
        has_col_data <- TRUE
      }
    }
    
    # Check if it's a list that might contain sample information
    else if (is.list(obj)) {
      cat("Object is a list\n")
      # Check for sample_table or targets
      if ("sample_table" %in% names(obj)) {
        sample_table <- obj$sample_table
        cat(paste0("✓ Found sample_table with ", nrow(sample_table), " samples and ", ncol(sample_table), " columns\n"))
        cat("  Column names: ", paste(colnames(sample_table), collapse = ", "), "\n")
        has_col_data <- TRUE
      } else if ("targets" %in% names(obj)) {
        targets <- obj$targets
        cat(paste0("✓ Found targets with ", nrow(targets), " samples and ", ncol(targets), " columns\n"))
        cat("  Column names: ", paste(colnames(targets), collapse = ", "), "\n")
        has_col_data <- TRUE
      } else {
        cat("✗ No sample_table or targets found in list\n")
      }
    }
    
    # Check if it's a matrix
    else if (is.matrix(obj)) {
      cat("Object is a matrix\n")
      cat(paste0("Dimensions: ", nrow(obj), "x", ncol(obj), "\n"))
      if (is.null(rownames(obj))) {
        cat("Row names: None\n")
      } else {
        cat(paste0("Row names: ", paste(c(head(rownames(obj), 3), "..."), collapse = ", "), "\n"))
      }
      if (is.null(colnames(obj))) {
        cat("Column names: None\n")
      } else {
        cat(paste0("Column names: ", paste(c(head(colnames(obj), 3), "..."), collapse = ", "), "\n"))
      }
      cat("✗ Matrix has no column annotations\n")
    }
    
    # Other types
    else {
      cat("✗ Unknown object type\n")
    }
    
    # Summary for this file
    if (has_col_data) {
      cat("✓ Column annotations found\n")
    } else {
      cat("✗ No column annotations found\n")
    }
    
    # Cleanup
    rm(obj)
    gc()
  } else {
    cat(paste0("❌ File not found: ", basename(file_path), "\n"))
  }
}

# Check for any other potential annotation files
cat("\n=== Checking for other annotation files ===\n")

# Look for sample_table or targets files
sample_files <- list.files(c(processed_full_dir, processed_dir), pattern = "sample_table|targets", ignore.case = TRUE, full.names = TRUE)
if (length(sample_files) > 0) {
  cat(paste0("Found ", length(sample_files), " sample annotation files:\n"))
  for (file in sample_files) {
    cat(paste0("  - ", basename(file), "\n"))
  }
} else {
  cat("No sample annotation files found\n")
}

# Summary
cat("\n=== Column annotations check completed ===\n")
sink()

# Memory cleanup
rm(list = ls())
gc()
