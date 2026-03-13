#!/usr/bin/env Rscript

# =========================================================================
# Script: 02_prep_anchors_65204_201955_45111.R
# Purpose: Extract 11-gene data from GSE65204 (IgE), GSE201955 (FeNO), and GSE45111 (Sputum eosinophils) anchor cohorts
# Author: Zuo
# Date: 2026-02-26
#
# Inputs:
#   - data/processed_full/expr_*__FULL.rds  # Full expression matrices
#   - data/processed_full/pheno_*__FULL.rds  # Phenotype data
#   - data/derived/locked_weights.csv  # Locked weights for 11 genes
#   - data/processed_full/expr_colname_to_gsm_mapping_GSE201955.csv  # GSE201955 sample mapping
#
# Outputs:
#   - data/processed_diet/*_diet_expr.rds  # Diet expression matrices
#   - data/processed_diet/*_pheno.rds  # Diet phenotype data
#   - output/logs/02_prep_anchors.log  # Log file
# =========================================================================

# Load necessary packages
library(tidyverse)

# Set paths
# Get base directory (project root)
base_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# Verify we're running from project root
if (!all(dir.exists(file.path(base_dir, c("data", "analysis", "data_preparation"))))) {
  stop("Please run this script from the project root directory (the folder containing data/, analysis/, data_preparation/).", call. = FALSE)
}

# Define directories
processed_full_dir <- file.path(base_dir, "data", "processed_full")
processed_diet_dir <- file.path(base_dir, "data", "processed_diet")
derived_dir <- file.path(base_dir, "data", "derived")
logs_dir <- file.path(base_dir, "output", "logs")

# Create output directories if they don't exist
dir.create(processed_diet_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

# Log file
log_file <- file.path(logs_dir, "02_prep_anchors.log")
sink(log_file, append = TRUE)

# Log base directory and paths
cat("\n=== Starting extraction of three anchor cohort datasets ===\n")
cat(paste0("Execution time: ", Sys.time(), "\n"))
cat(paste0("Base directory: ", base_dir, "\n"))
cat(paste0("Processed full directory: ", processed_full_dir, "\n"))
cat(paste0("Processed diet directory: ", processed_diet_dir, "\n"))
cat(paste0("Derived directory: ", derived_dir, "\n"))
cat(paste0("Logs directory: ", logs_dir, "\n"))

# Read locked weights to get the 11-gene list
locked_weights_file <- file.path(derived_dir, "locked_weights.csv")
if (!file.exists(locked_weights_file)) {
  stop("locked_weights.csv not found in data/derived/")
}
locked_weights <- read.csv(locked_weights_file)
locked_genes <- as.character(locked_weights$Gene)

# Assert: must have exactly 11 genes
stopifnot(length(locked_genes) == 11)
cat(paste0("Locked genes (11): ", paste(locked_genes, collapse = ", "), "\n"))

# Define cohorts to process
gse_ids <- c("GSE65204", "GSE201955", "GSE45111")

# Clinical phenotype mapping
clinical_mapping <- list(
  "GSE65204" = "lnige",
  "GSE201955" = "feno",
  "GSE45111" = "sputum_eos"
)

for (gse_id in gse_ids) {
  cat(paste0("\nProcessing cohort: ", gse_id, "\n"))
  
  # Read full expression matrix
  # Find the expression file with GPL info (case-insensitive)
  expr_files <- list.files(processed_full_dir, pattern = paste0("expr_.*", tolower(gse_id), ".*__FULL.rds"))
  if (length(expr_files) == 0) {
    # Try with uppercase
    expr_files <- list.files(processed_full_dir, pattern = paste0("expr_.*", gse_id, ".*__FULL.rds"))
    if (length(expr_files) == 0) {
      cat(paste0("⚠️  Missing expression file for ", gse_id, "\n"))
      next
    }
  }
  expr_file <- file.path(processed_full_dir, expr_files[1])
  
  # Read pheno data
  pheno_file <- file.path(processed_full_dir, paste0("pheno_", gse_id, "__FULL.rds"))
  if (!file.exists(pheno_file)) {
    cat(paste0("⚠️  Missing pheno file for ", gse_id, "\n"))
    next
  }
  
  cat(paste0("Using expression file: ", basename(expr_file), "\n"))
  cat(paste0("Using pheno file: ", basename(pheno_file), "\n"))
  
  # Read expression matrix
  expr <- readRDS(expr_file)
  
  # Read pheno data
  pheno <- readRDS(pheno_file)
  
  # Clean sample_id by removing quotes if present
  pheno$sample_id <- gsub('"', '', pheno$sample_id)
  
  # Special handling for GSE201955
  if (gse_id == "GSE201955") {
    # For GSE201955, the expression matrix has different sample IDs than the pheno file
    # Need to use the mapping table to rename sample IDs
    cat('Special handling for GSE201955: using mapping table to rename sample IDs\n')
    
    # Read mapping table
    mapping_file <- file.path(processed_full_dir, "expr_colname_to_gsm_mapping_GSE201955.csv")
    if (!file.exists(mapping_file)) {
      stop("Mapping file not found: expr_colname_to_gsm_mapping_GSE201955.csv")
    }
    mapping <- read.csv(mapping_file)
    
    # Clean up GSM IDs by removing any leading/trailing whitespace and "="
    mapping$GSM <- gsub("^\\s*=\\s*", "", mapping$GSM)
    mapping$GSM <- gsub("\\s*$", "", mapping$GSM)
    mapping$expr_colname <- gsub("\\s*$", "", mapping$expr_colname)
    
    # Verification 1: Check if mapping is complete and one-to-one
    cat('Verification 1: Checking mapping integrity\n')
    cat(paste0('Mapping rows: ', nrow(mapping), '\n'))
    cat(paste0('Expression columns: ', ncol(expr), '\n'))
    
    # Check mapping rows
    if (nrow(mapping) != 118) {
      cat(paste0('Error: Mapping has ', nrow(mapping), ' rows, expected 118\n'))
      stop("Mapping should have exactly 118 rows")
    }
    
    # Check for duplicate expression column names
    if (anyDuplicated(mapping$expr_colname) != 0) {
      cat('Error: Duplicate expression column names found in mapping\n')
      stop("Duplicate expression column names in mapping")
    }
    
    # Check for duplicate GSM IDs
    if (anyDuplicated(mapping$GSM) != 0) {
      cat('Error: Duplicate GSM IDs found in mapping\n')
      stop("Duplicate GSM IDs in mapping")
    }
    
    # Check if all mapping expression column names are in expression matrix (with or without X prefix)
    cat('First 10 expression column names: ', paste(head(colnames(expr), 10), collapse = ', '), '\n')
    cat('First 10 mapping expression column names: ', paste(head(mapping$expr_colname, 10), collapse = ', '), '\n')
    
    # Check for exact matches
    exact_matches <- mapping$expr_colname %in% colnames(expr)
    
    # Check for matches with X prefix
    mapping_colnames_with_x <- paste0('X', mapping$expr_colname)
    x_prefix_matches <- mapping_colnames_with_x %in% colnames(expr)
    
    # Combine matches
    all_matches <- exact_matches | x_prefix_matches
    
    missing_cols <- mapping$expr_colname[!all_matches]
    if (length(missing_cols) > 0) {
      cat(paste0('Error: Missing expression columns in matrix: ', paste(missing_cols, collapse = ', '), '\n'))
      stop("Some mapping expression column names not found in expression matrix")
    }
    
    # Update mapping to use the actual column names from the expression matrix
    mapping$expr_colname <- ifelse(exact_matches, mapping$expr_colname, mapping_colnames_with_x)
    
    cat('✓ Verification 1 passed: Mapping is complete and one-to-one\n')
    
    # Read soft pheno data to get GSM list
    soft_pheno_file <- file.path(processed_full_dir, "pheno_GSE201955_from_soft__RAW.rds")
    if (!file.exists(soft_pheno_file)) {
      stop("Soft pheno file not found: pheno_GSE201955_from_soft__RAW.rds")
    }
    soft_pheno <- readRDS(soft_pheno_file)
    soft_gsm_118 <- gsub("^\\s*=\\s*", "", soft_pheno$GSM)
    soft_gsm_118 <- gsub("\\s*$", "", soft_gsm_118)
    
    # Rename expression matrix column names using mapping
    cat('Renaming expression matrix column names using mapping\n')
    
    # Create a named vector for mapping
    map_vector <- setNames(mapping$GSM, mapping$expr_colname)
    
    # Rename columns
    new_colnames <- map_vector[colnames(expr)]
    cat(paste0('Number of NA values after mapping: ', sum(is.na(new_colnames)), '\n'))
    
    # Check if any NA values
    if (any(is.na(new_colnames))) {
      cat('Error: Some columns could not be mapped\n')
      stop("Failed to map all expression columns to GSM IDs")
    }
    
    colnames(expr) <- new_colnames
    
    # Verification 2: Check if renamed expression matrix matches soft GSM set
    cat('Verification 2: Checking consistency with soft GSM set\n')
    cat(paste0('Number of expression columns after rename: ', ncol(expr), '\n'))
    cat(paste0('Number of soft GSM IDs: ', length(soft_gsm_118), '\n'))
    
    intersection <- intersect(colnames(expr), soft_gsm_118)
    cat(paste0('Number of common GSM IDs: ', length(intersection), '\n'))
    
    # Show first 10 renamed expression columns
    cat('First 10 renamed expression columns: ', paste(head(colnames(expr), 10), collapse = ', '), '\n')
    # Show first 10 soft GSM IDs
    cat('First 10 soft GSM IDs: ', paste(head(soft_gsm_118, 10), collapse = ', '), '\n')
    
    # Find missing GSM IDs
    missing_gsm <- setdiff(colnames(expr), soft_gsm_118)
    if (length(missing_gsm) > 0) {
      cat(paste0('Missing GSM IDs in soft set: ', paste(missing_gsm, collapse = ', '), '\n'))
    }
    
    # Find extra GSM IDs in soft set
    extra_gsm <- setdiff(soft_gsm_118, colnames(expr))
    if (length(extra_gsm) > 0) {
      cat(paste0('Extra GSM IDs in soft set: ', paste(extra_gsm, collapse = ', '), '\n'))
    }
    
    if (length(intersection) != 118) {
      cat('Error: Renamed expression columns do not match soft GSM set\n')
      stop("Renamed expression columns do not match soft GSM set")
    }
    
    cat('✓ Verification 2 passed: Renamed expression columns match soft GSM set\n')
    
    # Check FeNO subset consistency
    cat('Checking FeNO subset consistency\n')
    common74 <- intersect(colnames(expr), pheno$sample_id)
    # Require at least 70 common samples, prefer exactly 74
    if (length(common74) < 70) {
      stop("Too few common samples between expression matrix and pheno file")
    } else if (length(common74) == nrow(pheno)) {
      cat('✓ FeNO subset: Exactly', length(common74), 'common samples found\n')
    } else {
      cat('⚠️  FeNO subset: Found', length(common74), 'common samples (expected', nrow(pheno), ')\n')
    }
  }
  
  # Assert: expression matrix must be matrix or dgCMatrix
  stopifnot(is.matrix(expr) || inherits(expr, "dgCMatrix"))
  
  # Assert: matrix must have row names and column names
  stopifnot(!is.null(rownames(expr)), !is.null(colnames(expr)))
  
  # Check for duplicate gene symbols
  if (any(duplicated(rownames(expr)))) {
    stop("Duplicate gene symbols detected; resolve in 00b mapping")
  }
  
  # Get original sample counts
  original_expr_samples <- ncol(expr)
  original_pheno_samples <- nrow(pheno)
  
  # Sample alignment: find common samples between expr and pheno
  expr_samples <- colnames(expr)
  pheno_samples <- pheno$sample_id
  common_samples <- intersect(expr_samples, pheno_samples)
  
  if (length(common_samples) == 0) {
    cat(paste0("⚠️  No common samples found between expr and pheno for ", gse_id, "\n"))
    next
  }
  
  # Subset expr and pheno to common samples
  expr <- expr[, common_samples]
  pheno <- pheno[pheno$sample_id %in% common_samples, ]
  
  # Reorder pheno to match expr sample order
  pheno <- pheno[match(colnames(expr), pheno$sample_id), ]
  
  # Add clinical anchor column
  clinical_col <- clinical_mapping[[gse_id]]
  if (!is.null(clinical_col) && clinical_col %in% colnames(pheno)) {
    cat(paste0("Using clinical anchor column: ", clinical_col, "\n"))
  } else {
    cat(paste0("⚠️  Clinical anchor column not found for ", gse_id, "\n"))
    clinical_col <- "Clinical_Anchor"
    pheno[[clinical_col]] <- NA
  }
  
  # Extract only locked genes
  present_genes <- intersect(locked_genes, rownames(expr))
  missing_genes <- setdiff(locked_genes, rownames(expr))
  
  # Subset expression matrix to present genes
  expr_diet <- expr[present_genes, ]
  
  # Reorder rows to match locked_genes order (for present genes)
  # Only keep genes that are present
  ordered_genes <- intersect(locked_genes, present_genes)
  expr_diet <- expr_diet[ordered_genes, ]
  
  # Calculate coverage
  coverage <- length(present_genes) / length(locked_genes)
  cat(paste0("Gene coverage: ", length(present_genes), "/", length(locked_genes), " (", round(coverage * 100, 1), "%)\n"))
  
  # Diagnostic logging for missing genes, especially for GSE65204
  if (length(missing_genes) > 0) {
    cat(paste0("Missing genes: ", paste(missing_genes, collapse = ", "), "\n"))
    
    # For GSE65204, perform detailed diagnostic
    if (gse_id == "GSE65204") {
      cat("\nDiagnostic for missing genes in GSE65204:\n")
      for (gene in missing_genes) {
        # Check 1: Exact match
        exact_match <- gene %in% rownames(expr)
        # Check 2: Case-insensitive match
        case_insensitive_match <- tolower(gene) %in% tolower(rownames(expr))
        # Check 3: Find similar names
        similar_names <- rownames(expr)[agrepl(gene, rownames(expr), ignore.case = TRUE)]
        
        cat(paste0("- ", gene, "\n"))
        cat(paste0("  Exact match: ", exact_match, "\n"))
        cat(paste0("  Case-insensitive match: ", case_insensitive_match, "\n"))
        if (length(similar_names) > 0) {
          cat(paste0("  Similar names found: ", paste(similar_names, collapse = ", "), "\n"))
        } else {
          cat("  No similar names found\n")
        }
      }
    }
  }
  
  # Save diet data
  expr_output <- file.path(processed_diet_dir, paste0(gse_id, "_diet_expr.rds"))
  pheno_output <- file.path(processed_diet_dir, paste0(gse_id, "_pheno.rds"))
  
  saveRDS(expr_diet, expr_output)
  saveRDS(pheno, pheno_output)
  
  # Logging
  cat(paste0("Original expression samples: ", original_expr_samples, "\n"))
  cat(paste0("Original pheno samples: ", original_pheno_samples, "\n"))
  cat(paste0("Common samples after alignment: ", length(common_samples), "\n"))
  cat(paste0("Samples dropped: ", original_expr_samples - length(common_samples), "\n"))
  cat(paste0("Diet expression matrix dimensions: ", nrow(expr_diet), " x ", ncol(expr_diet), "\n"))
  cat(paste0("Saved diet expression to: ", basename(expr_output), "\n"))
  cat(paste0("Saved pheno data to: ", basename(pheno_output), "\n"))
  
  # Memory cleanup
  rm(expr, pheno, expr_diet)
  gc()
}

cat("\n=== Extraction of three anchor cohort datasets completed ===\n")
sink()

# Final memory cleanup
rm(list = ls())
gc()