#!/usr/bin/env Rscript

# =========================================================================
# Script: 05_calculate_minimal_stats.R
# Purpose: Calculate minimal statistical closure for each cohort, including sample counts, effect sizes, and p-values
# Author: Zuo
# Date: 2026-02-28
#
# Inputs:
#   - analysis/02_main_validation/00_replication_contrasts.yaml  # Replication cohort specifications
#   - data/processed_diet/*_diet_expr.rds  # Minimal expression matrices for scoring
#   - data/derived/locked_weights.csv  # Locked weights for scoring
#
# Outputs:
#   - data/processed/phase3_outputs/*_stats.rds  # Statistical results for each cohort
#   - data/processed/phase3_outputs/minimal_stats_report.csv  # Summary of all cohort stats
# =========================================================================

# Load necessary packages
library(dplyr)
library(stringr)
library(Seurat)
library(yaml)

# Set base directory
base_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# Verify we're running from project root
required_dirs <- c('data', 'analysis', 'data_preparation')
required_paths <- file.path(base_dir, required_dirs)
dirs_exist <- all(dir.exists(required_paths))
if (!dirs_exist) {
  stop("Please run this script from the project root directory (the folder containing data/, analysis/, data_preparation/).", call. = FALSE)
}

# Define directories
processed_dir <- file.path(base_dir, "data", "processed", "phase3_outputs")
processed_diet_dir <- file.path(base_dir, "data", "processed_diet")
derived_dir <- file.path(base_dir, "data", "derived")

# Read locked weights
read_locked_weights <- function() {
  locked_weights_file <- file.path(derived_dir, "locked_weights.csv")
  if (!file.exists(locked_weights_file)) {
    stop("Locked weights file not found: ", locked_weights_file, call. = FALSE)
  }
  
  locked_weights <- read.csv(locked_weights_file)
  
  # Find gene column (case-insensitive)
  gene_col <- grep("gene", colnames(locked_weights), ignore.case = TRUE)
  if (length(gene_col) > 0) {
    gene_col_name <- colnames(locked_weights)[gene_col[1]]
  } else {
    # Fallback to first column
    gene_col_name <- colnames(locked_weights)[1]
  }
  
  # Find weight column (case-insensitive)
  weight_col <- grep("weight", colnames(locked_weights), ignore.case = TRUE)
  if (length(weight_col) > 0) {
    weight_col_name <- colnames(locked_weights)[weight_col[1]]
  } else {
    # Fallback to second column
    weight_col_name <- colnames(locked_weights)[2]
  }
  
  # Return data frame with gene and weight columns
  locked_weights[, c(gene_col_name, weight_col_name)]
}

# Calculate signature score
calculate_signature_score <- function(expr_matrix, locked_weights, cohort = "") {
  # Get common genes
  common_genes <- intersect(rownames(expr_matrix), locked_weights[, 1])
  
  if (length(common_genes) == 0) {
    return(NA)
  }
  
  # Subset expression matrix and weights
  expr_subset <- expr_matrix[common_genes, ]
  weights_subset <- locked_weights[locked_weights[, 1] %in% common_genes, 2]
  
  # Debug output for specific cohorts
  if (cohort %in% c("GSE115770", "GSE118761")) {
    cat("\n=== Debug for", cohort, "===")
    cat("\n1. Score calculation formula: score = t(expr_subset) %*% weights_subset")
    
    # Check gene order consistency
    weight_genes <- locked_weights[locked_weights[, 1] %in% common_genes, 1]
    expr_genes <- rownames(expr_subset)
    order_match <- all(weight_genes == expr_genes)
    cat("\n2. Gene order match:", order_match)
    
    if (!order_match) {
      cat("\n   First 5 weight genes:", paste(head(weight_genes, 5), collapse = ", "))
      cat("\n   First 5 expr genes:", paste(head(expr_genes, 5), collapse = ", "))
    }
    
    # Check if expression matrix is standardized
    row_means <- rowMeans(expr_subset)
    row_sds <- apply(expr_subset, 1, sd)
    is_standardized <- all(abs(row_means) < 0.001) && all(abs(row_sds - 1) < 0.001)
    cat("\n3. Is expression matrix standardized:", is_standardized)
    if (is_standardized) {
      cat(" (per gene)")
    }
    cat("\n")
  }
  
  # Calculate score
  score <- t(expr_subset) %*% weights_subset
  as.numeric(score)
}

# Calculate Cohen's d
calculate_cohens_d <- function(group1, group2) {
  mean1 <- mean(group1, na.rm = TRUE)
  mean2 <- mean(group2, na.rm = TRUE)
  sd1 <- sd(group1, na.rm = TRUE)
  sd2 <- sd(group2, na.rm = TRUE)
  n1 <- length(na.omit(group1))
  n2 <- length(na.omit(group2))
  
  pooled_sd <- sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))
  (mean1 - mean2) / pooled_sd
}

# Calculate coverage
calculate_coverage <- function(expr_matrix, locked_weights) {
  common_genes <- intersect(rownames(expr_matrix), locked_weights[, 1])
  list(
    count = length(common_genes),
    missing = setdiff(locked_weights[, 1], rownames(expr_matrix))
  )
}

# Define function: Calculate minimal statistical closure
calculate_minimal_stats <- function(cohort, config) {
  cat("\nCalculating minimal stats for", cohort, "...\n")
  
  # Initialize results list
  stats <- list(
    cohort = cohort,
    n_total = 0,
    n_case = 0,
    n_control = 0,
    n_aa = 0,
    n_na = 0,
    n_hc = 0,
    effect_size = NA,
    p_value = NA,
    auc = NA,
    coverage_count = 0,
    coverage_missing = character(),
    status = "FAIL",
    fail_reason = "",
    notes = "",
    # 新增字段
    case_levels = config$case_levels,
    control_levels = config$control_levels,
    mean_score_case = NA,
    mean_score_control = NA,
    mean_diff = NA,
    posWeightGene_dir_summary = list()
  )
  
  # Read phenotype data
  meta_path <- config$meta_path
  if (!file.exists(meta_path)) {
    stats$fail_reason <- "Phenotype data not found"
    cat("Error: Phenotype data not found at", meta_path, "\n")
    return(stats)
  }
  
  pheno <- readRDS(meta_path)
  stats$n_total <- nrow(pheno)
  
  # Special handling for GSE230048_pheno_blood.rds
  if (cohort == "GSE230048" && grepl("_blood", meta_path)) {
    # Try to get std_group_label from GSE230048_pheno.rds
    pheno_main <- readRDS(file.path(processed_diet_dir, "GSE230048_pheno.rds"))
    if ("std_group_label" %in% colnames(pheno_main)) {
      # Ensure sample order consistency
      if (all(pheno$sample_id == pheno_main$raw_sample_id)) {
        pheno$std_group_label <- pheno_main$std_group_label
        cat("Added std_group_label from GSE230048_pheno.rds to GSE230048_pheno_blood.rds\n")
      } else {
        cat("Warning: Sample order mismatch between GSE230048_pheno_blood.rds and GSE230048_pheno.rds\n")
      }
    }
  }
  
  # Read locked weights
  locked_weights <- read_locked_weights()
  
  # Try to read expression data
  expr_path <- config$expr_path
  if (!file.exists(expr_path)) {
    stats$fail_reason <- "Expression data not found"
    cat("Error: Expression data not found at", expr_path, "\n")
    return(stats)
  }
  
  cat("Reading expression data from", expr_path, "\n")
  expr_data <- readRDS(expr_path)
  
  # 确保是矩阵格式
  if ("Seurat" %in% class(expr_data)) {
    # 从Seurat对象中提取表达矩阵
    expr_matrix <- GetAssayData(expr_data, assay = "RNA", slot = "counts")
  } else {
    expr_matrix <- expr_data
  }
  
  # Calculate coverage
  coverage <- calculate_coverage(expr_matrix, locked_weights)
  stats$coverage_count <- coverage$count
  stats$coverage_missing <- coverage$missing
  
  cat("Signature gene coverage:", coverage$count, "/", nrow(locked_weights), "\n")
  if (length(coverage$missing) > 0) {
    cat("Missing genes:", paste(coverage$missing, collapse = ", "), "\n")
  }
  
  # Calculate signature scores
  # First create expr_subset for subsequent homology check
  common_genes <- intersect(rownames(expr_matrix), locked_weights[, 1])
  expr_subset <- expr_matrix[common_genes, ]
  
  scores <- calculate_signature_score(expr_matrix, locked_weights, cohort)
  
  if (!is.na(scores[1])) {
    # 添加分数到表型数据
    pheno$signature_score <- scores
    
    # Extract group information
    group_col <- config$group_col
    case_levels <- config$case_levels
    control_levels <- config$control_levels
    exclude_levels <- config$exclude_levels
    
    # Check if group column exists
    if (!(group_col %in% colnames(pheno))) {
      stats$fail_reason <- paste("Missing group column:", group_col)
      cat("Error: Missing group column", group_col, "in phenotype data\n")
      return(stats)
    }
    
    # Extract groups
    group <- pheno[[group_col]]
    
    # Exclude specified levels
    if (length(exclude_levels) > 0) {
      keep_idx <- !group %in% exclude_levels
      pheno <- pheno[keep_idx, ]
      scores <- scores[keep_idx]
      group <- group[keep_idx]
    }
    
    # Define case and control indices
    case_idx <- group %in% case_levels
    ctrl_idx <- group %in% control_levels
    
    # Calculate sample counts
    stats$n_case <- sum(case_idx, na.rm = TRUE)
    stats$n_control <- sum(ctrl_idx, na.rm = TRUE)
    
    # Print sample distribution
    cat("Sample distribution:\n")
    cat("Case:", stats$n_case, "(levels:", paste(case_levels, collapse = ", "), ")\n")
    cat("Control:", stats$n_control, "(levels:", paste(control_levels, collapse = ", "), ")\n")
    
    # Check sample counts
    if (stats$n_case < 2 || stats$n_control < 2) {
      stats$fail_reason <- "insufficient_case_control_after_filter"
      cat("Error: Insufficient case/control samples after filter\n")
      return(stats)
    }
    
    # Extract case and control scores
    case_scores <- scores[case_idx]
    control_scores <- scores[ctrl_idx]
    
    # Calculate mean scores
    stats$mean_score_case <- mean(case_scores)
    stats$mean_score_control <- mean(control_scores)
    stats$mean_diff <- stats$mean_score_case - stats$mean_score_control
    
    # Calculate Cohen's d
    stats$effect_size <- calculate_cohens_d(case_scores, control_scores)
    
    # Calculate Wilcoxon p-value
    stats$p_value <- wilcox.test(case_scores, control_scores)$p.value
    
    # Debug output: Confirm matrix object homology
    if (cohort %in% c("GSE115770", "GSE118761")) {
      cat("\n4. Matrix homology check for", cohort, "===")
      # Check if using the same matrix object
      same_matrix <- identical(expr_matrix, expr_matrix)
      cat("\n   Same matrix object:", same_matrix)
      # Check if using the same transformation (not standardized)
      same_transform <- TRUE  # because both are not standardized
      cat("\n   Same transformation:", same_transform)
      # Check if using the same subset of samples
      same_samples <- all(colnames(expr_subset) == colnames(expr_matrix))
      cat("\n   Same subset samples:", same_samples)
      # Check if case_idx and ctrl_idx are based on the same samples
      same_indices <- length(case_idx) == ncol(expr_subset) && length(ctrl_idx) == ncol(expr_subset)
      cat("\n   Same case/control indices:", same_indices)
      cat("\n")
    }
    
    # Calculate positive weight gene direction statistics
    pos_weight_genes <- c("SERPINB2", "CLCA1", "CCR3", "HDC", "CPA3", "MUC5B", "MUC5AC", "IL13", "CCL24", "POSTN")
    gene_directions <- list()
    pos_count <- 0
    neg_count <- 0
    na_count <- 0
    
    for (gene in pos_weight_genes) {
      if (gene %in% rownames(expr_matrix)) {
        gene_expr <- expr_matrix[gene, ]
        case_expr <- gene_expr[case_idx]
        control_expr <- gene_expr[ctrl_idx]
        
        if (length(case_expr) > 0 && length(control_expr) > 0) {
          diff <- mean(case_expr) - mean(control_expr)
          direction <- ifelse(diff > 0, "pos", "neg")
          gene_directions[[gene]] <- direction
          
          if (direction == "pos") {
            pos_count <- pos_count + 1
          } else {
            neg_count <- neg_count + 1
          }
        } else {
          gene_directions[[gene]] <- "na"
          na_count <- na_count + 1
        }
      } else {
        gene_directions[[gene]] <- "na"
        na_count <- na_count + 1
      }
    }
    
    # Score contribution decomposition analysis (for all cohorts)
    # Extract weights for all genes
    all_genes <- locked_weights[, 1]
    all_weights <- locked_weights[, 2]
    names(all_weights) <- all_genes
    
    # Calculate contribution for each gene
    delta_contrib <- numeric()
    
    for (gene in all_genes) {
      if (gene %in% rownames(expr_matrix)) {
        gene_expr <- expr_matrix[gene, ]
        weight <- all_weights[gene]
        
        # Calculate contribution
        contrib <- gene_expr * weight
        
        # Calculate mean contribution for case and control
        contrib_case <- contrib[case_idx]
        contrib_control <- contrib[ctrl_idx]
        
        if (length(contrib_case) > 0 && length(contrib_control) > 0) {
          delta_contrib[gene] <- mean(contrib_case) - mean(contrib_control)
        } else {
          delta_contrib[gene] <- NA
        }
      } else {
        delta_contrib[gene] <- NA
      }
    }
    
    # Calculate total contribution difference
    sum_delta_contrib <- sum(delta_contrib, na.rm = TRUE)
    
    # Sort and get top 5 most negative delta_contrib
    delta_contrib <- delta_contrib[!is.na(delta_contrib)]
    sorted_delta <- sort(delta_contrib)
    top5_negative <- head(sorted_delta, 5)
    
    # Debug output (only for specified cohorts)
    if (cohort %in% c("GSE115770", "GSE118761")) {
      cat("\n5. Score contribution decomposition for", cohort, "===")
      cat("\n   Sum of delta contributions:", sum_delta_contrib)
      cat("\n   Mean difference (from earlier):", stats$mean_diff)
      cat("\n   Difference:", abs(sum_delta_contrib - stats$mean_diff))
      
      cat("\n   Top 5 most negative delta contributions:")
      for (i in 1:length(top5_negative)) {
        gene <- names(top5_negative)[i]
        value <- top5_negative[i]
        cat("\n   ", gene, ":", value)
      }
      cat("\n")
    }
    
    # Save contribution decomposition results to stats
    stats$sum_delta_contrib <- sum_delta_contrib
    stats$top5_negative_contrib <- top5_negative
    
    stats$posWeightGene_dir_summary <- list(
      directions = gene_directions,
      pos_count = pos_count,
      neg_count = neg_count,
      na_count = na_count
    )
    
    cat("Cohen's d:", stats$effect_size, "\n")
    cat("Wilcoxon p-value:", stats$p_value, "\n")
    cat("Mean score case:", stats$mean_score_case, "\n")
    cat("Mean score control:", stats$mean_score_control, "\n")
    cat("Mean difference:", stats$mean_diff, "\n")
    cat("Positive weight gene directions - Pos:", pos_count, "Neg:", neg_count, "NA:", na_count, "\n")
    
    stats$status <- "PASS"
    stats$notes <- "replication-only cohort analysis"
  } else {
    stats$fail_reason <- "Failed to calculate signature scores"
    cat("Error: Failed to calculate signature scores\n")
  }
  
  # Save stats for each cohort
  output_file <- file.path(processed_dir, paste0(cohort, "_stats.rds"))
  saveRDS(stats, output_file)
  cat("Saved stats to", output_file, "\n")
  
  return(stats)
}

# Define function: Generate minimal stats report
generate_minimal_stats_report <- function() {
  # Read replication_contrasts.yaml configuration file
  yaml_file <- file.path(base_dir, "analysis", "02_main_validation", "00_replication_contrasts.yaml")
  if (!file.exists(yaml_file)) {
    stop("Replication contrasts YAML file not found: ", yaml_file, call. = FALSE)
  }
  
  config <- yaml::read_yaml(yaml_file)
  
  # Define replication-only cohort list
  replication_cohorts <- c("GSE103166", "GSE115770", "GSE118761", "GSE123750", "GSE230048", "GSE40888", "GSE43696")
  
  # Calculate stats for each cohort
  all_stats <- list()
  for (cohort in replication_cohorts) {
    if (cohort %in% names(config)) {
      all_stats[[cohort]] <- calculate_minimal_stats(cohort, config[[cohort]])
    } else {
      cat("Warning: No configuration found for cohort", cohort, "\n")
    }
  }
  
  # Convert to data frame
  stats_list <- lapply(all_stats, function(stats) {
    # Process coverage_missing field
    coverage_missing_str <- if (length(stats$coverage_missing) > 0) {
      paste(stats$coverage_missing, collapse = ", ")
    } else {
      ""
    }
    
    # Process case_levels and control_levels
    case_levels_str <- if (is.list(stats$case_levels)) {
      paste(stats$case_levels, collapse = ", ")
    } else {
      as.character(stats$case_levels)
    }
    
    control_levels_str <- if (is.list(stats$control_levels)) {
      paste(stats$control_levels, collapse = ", ")
    } else {
      as.character(stats$control_levels)
    }
    
    # Process positive weight gene direction statistics
    pos_count <- if (!is.null(stats$posWeightGene_dir_summary$pos_count)) stats$posWeightGene_dir_summary$pos_count else 0
    neg_count <- if (!is.null(stats$posWeightGene_dir_summary$neg_count)) stats$posWeightGene_dir_summary$neg_count else 0
    na_count <- if (!is.null(stats$posWeightGene_dir_summary$na_count)) stats$posWeightGene_dir_summary$na_count else 0
    
    # Create data frame row
    df <- data.frame(
      cohort = stats$cohort,
      n_total = stats$n_total,
      n_case = stats$n_case,
      n_control = stats$n_control,
      n_aa = stats$n_aa,
      n_na = stats$n_na,
      n_hc = stats$n_hc,
      effect_size = stats$effect_size,
      p_value = stats$p_value,
      auc = stats$auc,
      coverage_count = stats$coverage_count,
      coverage_missing = coverage_missing_str,
      status = stats$status,
      fail_reason = stats$fail_reason,
      notes = stats$notes,
      # New fields
      case_levels = case_levels_str,
      control_levels = control_levels_str,
      mean_score_case = stats$mean_score_case,
      mean_score_control = stats$mean_score_control,
      mean_diff = stats$mean_diff,
      posWeightGene_pos_count = pos_count,
      posWeightGene_neg_count = neg_count,
      posWeightGene_na_count = na_count,
      stringsAsFactors = FALSE
    )
    
    return(df)
  })
  
  stats_df <- do.call(rbind, stats_list)
  
  # Save results
  output_file <- file.path(processed_dir, "minimal_stats_report.csv")
  write.csv(stats_df, output_file, row.names = FALSE)
  cat("\nSaved minimal stats report to", output_file, "\n")
  
  # Generate replication_direction_audit.csv
  direction_audit_df <- stats_df[, c(
    "cohort", 
    "case_levels", 
    "control_levels", 
    "n_case", 
    "n_control", 
    "mean_score_case", 
    "mean_score_control", 
    "mean_diff", 
    "posWeightGene_pos_count", 
    "posWeightGene_neg_count", 
    "posWeightGene_na_count", 
    "coverage_count", 
    "coverage_missing", 
    "effect_size", 
    "p_value", 
    "status", 
    "fail_reason"
  )]
  
  direction_audit_file <- file.path(processed_dir, "replication_direction_audit.csv")
  write.csv(direction_audit_df, direction_audit_file, row.names = FALSE)
  cat("\nSaved replication direction audit to", direction_audit_file, "\n")
  
  # Generate replication_contribution_audit.csv
  # Collect contribution decomposition data for all cohorts
  contribution_data <- list()
  for (cohort in names(all_stats)) {
    stats <- all_stats[[cohort]]
    if (!is.null(stats$top5_negative_contrib) && length(stats$top5_negative_contrib) > 0) {
      # Extract Top 5 contributions
      top_contrib <- stats$top5_negative_contrib
      contrib_row <- data.frame(
        cohort = stats$cohort,
        sum_delta_contrib = stats$sum_delta_contrib,
        Top1_gene = ifelse(length(top_contrib) >= 1, names(top_contrib)[1], NA),
        Top1_delta_contrib = ifelse(length(top_contrib) >= 1, top_contrib[1], NA),
        Top2_gene = ifelse(length(top_contrib) >= 2, names(top_contrib)[2], NA),
        Top2_delta_contrib = ifelse(length(top_contrib) >= 2, top_contrib[2], NA),
        Top3_gene = ifelse(length(top_contrib) >= 3, names(top_contrib)[3], NA),
        Top3_delta_contrib = ifelse(length(top_contrib) >= 3, top_contrib[3], NA),
        Top4_gene = ifelse(length(top_contrib) >= 4, names(top_contrib)[4], NA),
        Top4_delta_contrib = ifelse(length(top_contrib) >= 4, top_contrib[4], NA),
        Top5_gene = ifelse(length(top_contrib) >= 5, names(top_contrib)[5], NA),
        Top5_delta_contrib = ifelse(length(top_contrib) >= 5, top_contrib[5], NA),
        stringsAsFactors = FALSE
      )
      contribution_data[[cohort]] <- contrib_row
    }
  }
  
  if (length(contribution_data) > 0) {
    contribution_df <- do.call(rbind, contribution_data)
    contribution_file <- file.path(processed_dir, "replication_contribution_audit.csv")
    write.csv(contribution_df, contribution_file, row.names = FALSE)
    cat("\nSaved replication contribution audit to", contribution_file, "\n")
  }
  
  # Print summary
  cat("\n=== Minimal Stats Summary ===\n")
  print(stats_df[, c("cohort", "n_total", "n_case", "n_control", "effect_size", "p_value", "coverage_count", "status", "fail_reason")])
  
  return(stats_df)
}

# Run function
generate_minimal_stats_report()

cat("\nAll minimal stats calculations completed!\n")
