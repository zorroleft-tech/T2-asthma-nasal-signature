#!/usr/bin/env Rscript

# Script 01: Rigorous scoring functions based on locked weights
# Purpose: Calculate gene scores using fixed coefficients from data/derived/locked_weights.csv

# Load necessary packages
library(tidyverse)

# Read locked weights
read_locked_weights <- function(weights_file = NULL) {
  if (is.null(weights_file)) {
    # Default path
    base_dir <- getwd()
    weights_file <- file.path(base_dir, "data", "derived", "locked_weights.csv")
  }
  
  if (!file.exists(weights_file)) {
    stop(paste0("Locked weights file does not exist: ", weights_file))
  }
  
  weights <- read_csv(weights_file)
  cat(paste0("Reading locked weights, number of genes: ", nrow(weights), "\n"))
  
  return(weights)
}

# Calculate gene scores
calculate_score <- function(expr_matrix, weights, gene_col = "Gene_Symbol", weight_col = "weight") {
  # Ensure consistent gene names
  colnames(expr_matrix)[1] <- gene_col
  
  # Keep only genes present in weights file
  expr_matrix_filtered <- expr_matrix %>%
    filter(!!sym(gene_col) %in% weights[[gene_col]])
  
  # Merge with weights
  # Check for duplicate column names before joining
  common_cols <- intersect(colnames(expr_matrix_filtered), colnames(weights))
  common_cols <- setdiff(common_cols, gene_col)  # Exclude the join column
  
  if (length(common_cols) > 0) {
    cat(sprintf("  Renaming duplicate columns: %s\n", paste(common_cols, collapse = ", ")))
    # Rename columns in weights to avoid duplicates
    for (col in common_cols) {
      colnames(weights)[colnames(weights) == col] <- paste0(col, "_weight")
    }
  }
  
  expr_matrix_weighted <- expr_matrix_filtered %>%
    inner_join(weights, by = gene_col)
  
  # Calculate weighted scores
  # Assume expression matrix structure: Gene_Symbol, sample1, sample2, ...
  sample_cols <- setdiff(colnames(expr_matrix_weighted), c(gene_col, weight_col))
  
  scores <- map_dfc(sample_cols, function(sample_col) {
    expr <- expr_matrix_weighted[[sample_col]]
    w <- expr_matrix_weighted[[weight_col]]
    sum(expr * w, na.rm = TRUE)
  })
  
  colnames(scores) <- sample_cols
  
  # Convert to data frame
  scores_df <- data.frame(
    sample_id = sample_cols,
    score = unlist(scores)
  )
  
  return(scores_df)
}

# Normalize scores (optional)
normalize_score <- function(scores_df, method = "z-score") {
  if (method == "z-score") {
    scores_df$score_normalized <- scale(scores_df$score)[, 1]
  } else if (method == "min-max") {
    scores_df$score_normalized <- (scores_df$score - min(scores_df$score, na.rm = TRUE)) / 
      (max(scores_df$score, na.rm = TRUE) - min(scores_df$score, na.rm = TRUE))
  }
  
  return(scores_df)
}

# Batch calculate scores for multiple cohorts
batch_calculate_score <- function(expr_files, weights_file = NULL) {
  # Read weights
  weights <- read_locked_weights(weights_file)
  
  # Store results
  all_scores <- list()
  
  for (expr_file in expr_files) {
    cat(paste0("Processing file: ", basename(expr_file), "\n"))
    
    # Read expression matrix
    expr_matrix <- readRDS(expr_file)
    
    # Calculate scores
    scores <- calculate_score(expr_matrix, weights)
    
    # Add cohort information
    cohort_name <- tools::file_path_sans_ext(basename(expr_file))
    scores$cohort <- cohort_name
    
    all_scores[[cohort_name]] <- scores
  }
  
  # Combine all results
  combined_scores <- bind_rows(all_scores)
  
  return(combined_scores)
}

# Example usage
# source("R/01_scoring_logic.R")
# weights <- read_locked_weights()
# expr_matrix <- readRDS("data/processed/GSE152004_discovery_expr.rds")
# scores <- calculate_score(expr_matrix, weights)
# scores_normalized <- normalize_score(scores)
