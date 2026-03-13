#!/usr/bin/env Rscript

# =========================================================================
# Script: validate_outputs_only.R
# Purpose: Validate output files for quality and integrity
# Author: Zuo
# Date: 2026-02-26
#
# Inputs:
#   - data/processed/*__FULL.rds  # Processed full expression matrices
#
# Outputs:
#   - Console output with validation results
# =========================================================================

# Load required packages
library(tidyverse)
library(data.table)
library(org.Hs.eg.db)

# Set paths
# Get base directory (project root)
base_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# Verify we're running from project root
if (!all(dir.exists(file.path(base_dir, c("data", "analysis", "data_preparation")))) {
  stop("Please run this script from the project root directory (the folder containing data/, analysis/, data_preparation/).", call. = FALSE)
}

# Define directories
processed_dir <- file.path(base_dir, "data", "processed")
logs_dir <- file.path(base_dir, "output", "logs")

# Function to determine ID type
determine_id_type <- function(ids) {
  if (length(ids) == 0) {
    return("UNKNOWN")
  }
  
  # Get a representative sample of IDs
  sample_size <- min(2000, length(ids))
  sample_ids <- sample(ids, sample_size)
  sample_ids <- as.character(sample_ids)
  
  # Count gene-like vs probe-like IDs
  gene_count <- 0
  entrez_count <- 0
  ensg_count <- 0
  probe_count <- 0
  
  is_probe <- function(x) {
    grepl("^ILMN_", x) ||
    grepl("^A_\\d+_", x) ||                 # Agilent
    grepl("_at$|_s_at$|_x_at$|_a_at$", x) || # Affy
    grepl("^\\d{6,}$", x)                   # 纯数字长串（常见探针/探针集）
  }
  
  is_entrez <- function(x) {
    grepl("^[0-9]+", x)  # 纯数字，全匹配
  }
  
  is_ensg <- function(x) {
    grepl("^ENSG", x)  # ENSG IDs
  }
  
  for (id in sample_ids) {
    id <- as.character(id)
    if (is_entrez(id)) {
      entrez_count <- entrez_count + 1
    } else if (is_probe(id)) {
      probe_count <- probe_count + 1
    } else if (is_ensg(id)) {
      ensg_count <- ensg_count + 1
    } else {
      gene_count <- gene_count + 1
    }
  }
  
  # Calculate percentages
  total_count <- length(sample_ids)
  gene_pct <- gene_count / total_count
  entrez_pct <- entrez_count / total_count
  ensg_pct <- ensg_count / total_count
  probe_pct <- probe_count / total_count
  
  # Determine the dominant ID type with threshold
  if (entrez_pct > 0.5) {
    return("ENTREZ")
  } else if (ensg_pct > 0.5) {
    return("ENSG")
  } else if (probe_pct > 0.5) {
    return("PROBE")
  } else if (gene_pct > 0.3) {
    # Use org.Hs.eg.db to validate if these are real HGNC symbols
    suppressPackageStartupMessages(library(AnnotationDbi))
    
    tryCatch({
      # Use select to check if IDs are valid symbols
      result <- select(org.Hs.eg.db, keys = sample_ids, keytype = "SYMBOL", columns = "ENTREZID")
      valid_symbols <- result[!is.na(result$ENTREZID), "SYMBOL"]
      hit_rate <- length(valid_symbols) / length(sample_ids)
      
      if (hit_rate > 0.3) {
        return("SYMBOL")
      } else {
        return("PROBE")
      }
    }, error = function(e) {
      # If AnnotationDbi fails, fall back to probe
      return("PROBE")
    })
  } else {
    # If no type dominates, use the original logic
    max_count <- max(ensg_count, gene_count, probe_count, entrez_count)
    if (entrez_count == max_count && entrez_count > 0) {
      return("ENTREZ")
    } else if (ensg_count == max_count && ensg_count > 0) {
      return("ENSG")
    } else if (gene_count == max_count && gene_count > 0) {
      # Validate gene symbols with org.Hs.eg.db
      suppressPackageStartupMessages(library(AnnotationDbi))
      
      tryCatch({
        result <- select(org.Hs.eg.db, keys = sample_ids, keytype = "SYMBOL", columns = "ENTREZID")
        valid_symbols <- result[!is.na(result$ENTREZID), "SYMBOL"]
        hit_rate <- length(valid_symbols) / length(sample_ids)
        
        if (hit_rate > 0.3) {
          return("SYMBOL")
        } else {
          return("PROBE")
        }
      }, error = function(e) {
        return("PROBE")
      })
    } else if (probe_count > 0) {
      return("PROBE")
    } else {
      return("UNKNOWN")
    }
  }
}

# Function to validate output files (memory-efficient version)
validate_outputs <- function() {
  cat("\n=== Validating output files ===\n")
  
  # Get all FULL files
  full_files <- list.files(processed_dir, pattern = "__FULL\\.rds$", full.names = TRUE)
  
  all_valid <- TRUE
  
  # Limit validation to first 10 files to reduce memory usage
  validation_files <- head(full_files, 10)
  if (length(full_files) > 10) {
    cat(paste0("⚠️  Validating first 10 of ", length(full_files), " files to reduce memory usage\n"))
  }
  
  for (full_file in validation_files) {
    cat(paste0("\nValidating: ", basename(full_file), "\n"))
    
    # Read the file - use a tryCatch to handle errors
    expr <- tryCatch({
      readRDS(full_file)
    }, error = function(e) {
      cat("❌ Failed to read file: ", conditionMessage(e), "\n")
      all_valid <<- FALSE
      return(NULL)
    })
    
    if (is.null(expr)) next
    
    # Check if it's a matrix
    if (!is.matrix(expr)) {
      cat("❌ Not a matrix\n")
      all_valid <- FALSE
    } else {
      # Check dimensions (lightweight check)
      cat(paste0("  Dimensions: ", nrow(expr), "x", ncol(expr), "\n"))
      
      # Check number of rows
      if (nrow(expr) < 2000) {
        cat(paste0("❌ Insufficient rows: ", nrow(expr), " < 2000\n"))
        all_valid <- FALSE
      }
      
      # Check row names (lightweight check - only first 100)
      row_names <- head(rownames(expr), 100)
      final_id_type <- determine_id_type(row_names)
      if (final_id_type != "SYMBOL") {
        cat(paste0("❌ Final ID type is not SYMBOL: ", final_id_type, "\n"))
        all_valid <- FALSE
      }
      
      # Check for duplicate row names in sample
      if (any(duplicated(row_names))) {
        cat("❌ Duplicate row names in sample\n")
        all_valid <- FALSE
      }
      
      # Check for empty or '---' row names in sample
      if (any(row_names == "")) {
        cat("❌ Empty row names in sample\n")
        all_valid <- FALSE
      }
      if (any(row_names == "---")) {
        cat("❌ Row names contain '---' in sample\n")
        all_valid <- FALSE
      }
    }
    
    # Clean up immediately to free memory
    rm(expr)
    gc()
    
    if (all_valid) {
      cat("✓ Valid\n")
    }
  }
  
  if (all_valid) {
    cat("\n=== All validated output files are valid ===\n")
  } else {
    cat("\n=== Some output files are invalid ===\n")
    stop("Validation failed")
  }
}

# Validate outputs
validate_outputs()
