# =========================================================================
# Script: 04_finalize_model.R
# Purpose: Determine final gene count based on EPV and biological quality
# Author: Zuo
# Date: 2024-01-01
#
# Inputs:
#   - data/derived/stability_selection_results_final.rds  # Stability selection results
#   - data/derived/limma_thresholds_report.csv  # Limma thresholds report
#   - data/processed_diet/GSE152004_diet_expr.rds  # Diet expression data
#   - data/derived/GSE152004_T2_labels.csv  # T2 status labels
#   - data/raw/candidate_genes_annotated_frozen.csv  # Annotated candidate genes
#
# Outputs:
#   - data/derived/final_signature_genes.csv  # Final signature genes
#   - data/derived/final_signature_summary.csv  # Final signature summary
#   - data/derived/epv_stats.csv  # EPV statistics
#   - data/derived/gene_correlation_matrix.rds  # Gene correlation matrix
#   - data/derived/locked_weights.csv  # Locked model weights
#   - data/derived/gene_panel_final.csv  # Final gene panel
#   - output/logs/04_finalize_model.log  # Log file
# =========================================================================

# Set project root using getwd()
base_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# Validate project root by checking for key directories
if (!dir.exists(file.path(base_dir, "data")) || !dir.exists(file.path(base_dir, "analysis")) || !dir.exists(file.path(base_dir, "data_preparation"))) {
  stop("Please run this script from the project root directory (the folder containing data/, analysis/, data_preparation/).", call. = FALSE)
}

# Define paths using file.path()
data_processed_diet_dir <- file.path(base_dir, "data", "processed_diet")
data_derived_dir <- file.path(base_dir, "data", "derived")
data_raw_dir <- file.path(base_dir, "data", "raw")
logs_dir <- file.path(base_dir, "output", "logs")

# Create output directories if not exist
dir.create(data_derived_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

# Set up logging
log_file <- file.path(logs_dir, "04_finalize_model.log")
sink(log_file, append = TRUE)

# Print base directory and paths to log
cat(sprintf("Base directory: %s\n", base_dir))
cat("Input files:\n")
cat(sprintf("  - %s\n", file.path(data_derived_dir, "stability_selection_results_final.rds")))
cat(sprintf("  - %s\n", file.path(data_derived_dir, "limma_thresholds_report.csv")))
cat(sprintf("  - %s\n", file.path(data_processed_diet_dir, "GSE152004_diet_expr.rds")))
cat(sprintf("  - %s\n", file.path(data_derived_dir, "GSE152004_T2_labels.csv")))
cat(sprintf("  - %s\n", file.path(data_raw_dir, "candidate_genes_annotated_frozen.csv")))
cat("Output files:\n")
cat(sprintf("  - %s\n", file.path(data_derived_dir, "final_signature_genes.csv")))
cat(sprintf("  - %s\n", file.path(data_derived_dir, "final_signature_summary.csv")))
cat(sprintf("  - %s\n", file.path(data_derived_dir, "epv_stats.csv")))
cat(sprintf("  - %s\n", file.path(data_derived_dir, "gene_correlation_matrix.rds")))
cat(sprintf("  - %s\n", file.path(data_derived_dir, "locked_weights.csv")))
cat(sprintf("  - %s\n", file.path(data_derived_dir, "gene_panel_final.csv")))
cat(sprintf("  - %s\n", log_file))

library(dplyr)
library(readr)
library(glmnet)
library(org.Hs.eg.db)

# Determine ID type function
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
    grepl("^[0-9]+$", x)  # 纯数字，全匹配
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
      return("PROBE")
    })
  } else {
    return("UNKNOWN")
  }
}

# Define paths using base directory
DATA_PROCESSED_DIET_DIR <- data_processed_diet_dir
DATA_DERIVED_DIR <- data_derived_dir
DATA_RAW_DIR <- data_raw_dir
LOG_DIR <- logs_dir

cat("=================================================================\n")
cat("Phase 6: Final Model Decision\n")
cat("=================================================================\n\n")

# ============================================================================
# Step 1: Load data
# ============================================================================

cat("[Step 1] Loading data...\n")

# Load stability selection results - priority: final > primary/fallback
final_results_file <- file.path(DATA_DERIVED_DIR, "stability_selection_results_final.rds")
thresholds_report <- read_csv(file.path(DATA_DERIVED_DIR, "limma_thresholds_report.csv"))
fallback_triggered <- thresholds_report$fallback_triggered[1]

if (file.exists(final_results_file)) {
  cat("  Using final stability selection results (fixed interface)\n")
  stability_results <- readRDS(final_results_file)
} else if (fallback_triggered) {
  cat("  Using fallback stability selection results\n")
  stability_results <- readRDS(file.path(DATA_DERIVED_DIR, "stability_selection_results_fallback.rds"))
} else {
  cat("  Using primary stability selection results\n")
  stability_results <- readRDS(file.path(DATA_DERIVED_DIR, "stability_selection_results_primary.rds"))
}

final_genes <- stability_results$final_genes
selection_freq <- stability_results$selection_freq

# Load expression data
expr_diet <- readRDS(file.path(DATA_PROCESSED_DIET_DIR, "GSE152004_diet_expr.rds"))

# Ensure gene × sample format
row_id_type <- determine_id_type(rownames(expr_diet))
col_id_type <- determine_id_type(colnames(expr_diet))

if (row_id_type != "SYMBOL" && col_id_type == "SYMBOL") {
  cat("  Transposing matrix to gene × sample format...\n")
  expr_diet <- t(expr_diet)
  
  # Update IDs after transpose
  row_id_type <- determine_id_type(rownames(expr_diet))
  col_id_type <- determine_id_type(colnames(expr_diet))
}

# Assertions
stopifnot(determine_id_type(rownames(expr_diet)) == "SYMBOL")
stopifnot(ncol(expr_diet) == 695)  # Should match full sample count

t2_labels_df <- read_csv(file.path(DATA_DERIVED_DIR, "GSE152004_T2_labels.csv"))

# Extract selected genes expression
expr_selected <- expr_diet[final_genes, ]

# Calculate EPV (Events Per Variable)
total_samples <- nrow(t2_labels_df)
t2_high_samples <- sum(t2_labels_df$T2_Status == "T2-high")
t2_low_samples <- sum(t2_labels_df$T2_Status == "T2-low")
min_events <- min(t2_high_samples, t2_low_samples)

cat(sprintf("\nSample distribution:\n"))
cat(sprintf("  Total samples: %d\n", total_samples))
cat(sprintf("  T2-high: %d\n", t2_high_samples))
cat(sprintf("  T2-low: %d\n", t2_low_samples))
cat(sprintf("  Minimum events: %d\n", min_events))

# ============================================================================
# Step 2: EPV Decision Tree
# ============================================================================

cat("\n[Step 2] EPV-based gene count decision...\n")

# Pre-defined EPV thresholds
EPV_THRESHOLDS <- list(
  minimal = 5,    # Absolute minimum
  optimal = 10,   # Optimal
  ideal = 15      # Ideal
)

# Calculate current EPV
current_epv <- min_events / length(final_genes)

cat(sprintf("\nCurrent status:\n"))
cat(sprintf("  Stability-selected genes: %d\n", length(final_genes)))
cat(sprintf("  Current EPV: %.2f\n", current_epv))

# EPV-based decision
if (current_epv >= EPV_THRESHOLDS$ideal) {
  # Can potentially include more genes
  decision <- "Ideal EPV - can include all stable genes"
  recommended_count <- length(final_genes)
  
} else if (current_epv >= EPV_THRESHOLDS$optimal) {
  # Maintain current gene count
  decision <- "Optimal EPV - maintain current gene count"
  recommended_count <- length(final_genes)
  
} else if (current_epv >= EPV_THRESHOLDS$minimal) {
  # Acceptable, but could be improved
  decision <- "Acceptable EPV - may consider slight reduction"
  recommended_count <- length(final_genes)
  
} else {
  # Need to reduce gene count
  decision <- "Low EPV - need to reduce gene count"
  # Calculate maximum gene count for minimal EPV
  max_genes_min_epv <- floor(min_events / EPV_THRESHOLDS$minimal)
  recommended_count <- max(9, min(max_genes_min_epv, 15))  # 12±3 range
}

cat(sprintf("\nDecision: %s\n", decision))
cat(sprintf("Recommended gene count: %d\n", recommended_count))

# ============================================================================
# Step 3: Biological Quality Checks
# ============================================================================

cat("\n[Step 3] Biological quality checks...\n")

# 1. Pathway coverage analysis
cat("  1. Pathway coverage analysis...\n")

# Load candidate genes with annotations
candidates <- read_csv(file.path(DATA_RAW_DIR, "candidate_genes_annotated_frozen.csv"))
selected_annotations <- candidates %>% 
  filter(Gene %in% final_genes)

# Count pathways
pathway_count <- selected_annotations %>% 
  group_by(Primary_Pathway) %>% 
  summarise(Count = n()) %>% 
  arrange(desc(Count))

cat(sprintf("  Pathway coverage: %d pathways\n", nrow(pathway_count)))

# 2. Drug target check
cat("  2. Drug target check...\n")
drug_targets <- selected_annotations %>% 
  filter(!is.na(Drug_Target) & Drug_Target != "None" & Drug_Target != "No")
cat(sprintf("  Drug targets: %d/%d genes\n", nrow(drug_targets), nrow(selected_annotations)))

# Show detailed drug target information
if (nrow(drug_targets) > 0) {
  cat("  Detailed drug target information:\n")
  drug_info <- drug_targets %>% 
    dplyr::select(Gene, Drug_Target)
  print(drug_info)
}

# 3. Literature support check
cat("  3. Literature support check...\n")
# Use PubMed count as literature support indicator
literature_support <- selected_annotations %>% 
  filter(PubMed_Count > 0)
cat(sprintf("  Literature supported genes: %d/%d genes\n", nrow(literature_support), nrow(selected_annotations)))

# ============================================================================
# Step 4: Gene Correlation Analysis
# ============================================================================

cat("\n[Step 4] Gene correlation analysis...\n")

# Extract selected genes expression
expr_selected_genes <- expr_diet[final_genes, ]

# Calculate gene-gene correlation (transpose to get gene × gene correlation)
cor_matrix <- cor(t(expr_selected_genes), method = "pearson", use = "pairwise.complete.obs")

# High correlation gene pairs
high_cor_threshold <- 0.7
high_cor_pairs <- which(cor_matrix > high_cor_threshold & cor_matrix < 1, arr.ind = TRUE)

# Count high correlation pairs (avoid duplicate counting)
cat(sprintf("  Highly correlated gene pairs (r > %.1f): %d\n", 
            high_cor_threshold, nrow(high_cor_pairs)/2))

# Show detailed high correlation pairs
if (nrow(high_cor_pairs) > 0) {
  cat("  Detailed high correlation pairs:")
  # Remove duplicate pairs
  unique_pairs <- high_cor_pairs[high_cor_pairs[,1] < high_cor_pairs[,2], ]
  for (i in 1:nrow(unique_pairs)) {
    gene1 <- rownames(cor_matrix)[unique_pairs[i,1]]
    gene2 <- rownames(cor_matrix)[unique_pairs[i,2]]
    corr_value <- cor_matrix[gene1, gene2]
    cat(sprintf("\n    %s - %s: %.3f", gene1, gene2, corr_value))
  }
  cat("\n")
}

# ============================================================================
# Step 5: Final Gene Selection
# ============================================================================

cat("\n[Step 5] Final gene selection...\n")

# Start with all stability-selected genes
final_gene_signature <- final_genes

# Apply pre-registered rule for high correlation (|r| > 0.85)
pre_reg_cor_threshold <- 0.85
high_cor_pairs_pre_reg <- which(abs(cor_matrix) > pre_reg_cor_threshold & abs(cor_matrix) < 1, arr.ind = TRUE)

if (nrow(high_cor_pairs_pre_reg) > 0) {
  cat(sprintf("  Applying pre-registered high correlation rule (|r| > %.2f)...\n", pre_reg_cor_threshold))
  
  # Remove duplicate pairs
  unique_pairs_pre_reg <- high_cor_pairs_pre_reg[high_cor_pairs_pre_reg[,1] < high_cor_pairs_pre_reg[,2], ]
  
  # Process each high correlation pair
  genes_to_remove <- c()
  
  for (i in 1:nrow(unique_pairs_pre_reg)) {
    gene1 <- rownames(cor_matrix)[unique_pairs_pre_reg[i,1]]
    gene2 <- rownames(cor_matrix)[unique_pairs_pre_reg[i,2]]
    corr_value <- cor_matrix[gene1, gene2]
    
    # Compare selection frequencies and remove the one with lower frequency
    freq1 <- selection_freq[gene1]
    freq2 <- selection_freq[gene2]
    
    if (freq1 > freq2) {
      genes_to_remove <- c(genes_to_remove, gene2)
      cat(sprintf("    Removing %s (frequency: %.3f) due to high correlation with %s (frequency: %.3f, r=%.3f)\n", 
                  gene2, freq2, gene1, freq1, corr_value))
    } else if (freq2 > freq1) {
      genes_to_remove <- c(genes_to_remove, gene1)
      cat(sprintf("    Removing %s (frequency: %.3f) due to high correlation with %s (frequency: %.3f, r=%.3f)\n", 
                  gene1, freq1, gene2, freq2, corr_value))
    } else {
      # If frequencies are equal, remove the one that appears later in the list
      if (which(final_genes == gene1) > which(final_genes == gene2)) {
        genes_to_remove <- c(genes_to_remove, gene1)
        cat(sprintf("    Removing %s due to high correlation with %s (equal frequency, r=%.3f)\n", 
                    gene1, gene2, corr_value))
      } else {
        genes_to_remove <- c(genes_to_remove, gene2)
        cat(sprintf("    Removing %s due to high correlation with %s (equal frequency, r=%.3f)\n", 
                    gene2, gene1, corr_value))
      }
    }
  }
  
  # Remove duplicate genes to remove
  genes_to_remove <- unique(genes_to_remove)
  
  # Remove the genes from the signature
  if (length(genes_to_remove) > 0) {
    final_gene_signature <- setdiff(final_gene_signature, genes_to_remove)
    cat(sprintf("  Removed %d genes due to high correlation\n", length(genes_to_remove)))
    
    # Export collinearity pruning log
    collinearity_log <- data.frame(
      gene_removed = character(),
      gene_kept = character(),
      corr_value = numeric(),
      freq_removed = numeric(),
      freq_kept = numeric(),
      rule = character(),
      stringsAsFactors = FALSE
    )
    
    # Re-process pairs to build log
    for (i in 1:nrow(unique_pairs_pre_reg)) {
      gene1 <- rownames(cor_matrix)[unique_pairs_pre_reg[i,1]]
      gene2 <- rownames(cor_matrix)[unique_pairs_pre_reg[i,2]]
      corr_value <- cor_matrix[gene1, gene2]
      
      # Compare selection frequencies
      freq1 <- selection_freq[gene1]
      freq2 <- selection_freq[gene2]
      
      if (freq1 > freq2) {
        collinearity_log <- rbind(collinearity_log, data.frame(
          gene_removed = gene2,
          gene_kept = gene1,
          corr_value = corr_value,
          freq_removed = freq2,
          freq_kept = freq1,
          rule = "abs(r)>0.85 keep higher frequency",
          stringsAsFactors = FALSE
        ))
      } else if (freq2 > freq1) {
        collinearity_log <- rbind(collinearity_log, data.frame(
          gene_removed = gene1,
          gene_kept = gene2,
          corr_value = corr_value,
          freq_removed = freq1,
          freq_kept = freq2,
          rule = "abs(r)>0.85 keep higher frequency",
          stringsAsFactors = FALSE
        ))
      } else {
        # If frequencies are equal, remove the one that appears later in the list
        if (which(final_genes == gene1) > which(final_genes == gene2)) {
          collinearity_log <- rbind(collinearity_log, data.frame(
            gene_removed = gene1,
            gene_kept = gene2,
            corr_value = corr_value,
            freq_removed = freq1,
            freq_kept = freq2,
            rule = "abs(r)>0.85 keep higher frequency",
            stringsAsFactors = FALSE
          ))
        } else {
          collinearity_log <- rbind(collinearity_log, data.frame(
            gene_removed = gene2,
            gene_kept = gene1,
            corr_value = corr_value,
            freq_removed = freq2,
            freq_kept = freq1,
            rule = "abs(r)>0.85 keep higher frequency",
            stringsAsFactors = FALSE
          ))
        }
      }
    }
    
    # Save collinearity pruning log
    write_csv(collinearity_log, file.path(DATA_DERIVED_DIR, "collinearity_pruning_log.csv"))
    cat(sprintf("  Collinearity pruning log saved: %s/collinearity_pruning_log.csv\n", DATA_DERIVED_DIR))
  }
}

# If recommended count is less than current, select top N by selection frequency
if (recommended_count < length(final_gene_signature)) {
  # Sort genes by selection frequency
  sorted_genes <- names(sort(selection_freq[final_gene_signature], decreasing = TRUE))
  final_gene_signature <- sorted_genes[1:recommended_count]
  
  cat(sprintf("  Reduced to %d genes based on selection frequency\n", recommended_count))
} else {
  # If more than 15, limit to top 15
  if (length(final_gene_signature) > 15) {
    sorted_genes <- names(sort(selection_freq[final_gene_signature], decreasing = TRUE))
    final_gene_signature <- sorted_genes[1:15]
    cat("  Limited to 15 genes (upper bound of 12±3)\n")
  }
  
  # If less than 9, try to include more from stability selection
  if (length(final_gene_signature) < 9) {
    # Get next most stable genes
    all_genes_sorted <- names(sort(selection_freq, decreasing = TRUE))
    additional_genes <- setdiff(all_genes_sorted, final_gene_signature)[1:(9 - length(final_gene_signature))]
    final_gene_signature <- c(final_gene_signature, additional_genes)
    cat("  Extended to 9 genes (lower bound of 12±3)\n")
  }
}

# ============================================================================
# Step 6: Export final signature
# ============================================================================

cat("\n[Step 6] Exporting final signature...\n")

# Create final signature dataframe
final_signature_df <- data.frame(
  Gene = final_gene_signature,
  Selection_Frequency = selection_freq[final_gene_signature],
  stringsAsFactors = FALSE
) %>% 
  arrange(desc(Selection_Frequency)) %>%
  mutate(Rank = 1:n())

# Export final signature
write_csv(final_signature_df, file.path(DATA_DERIVED_DIR, "final_signature_genes.csv"))

# Export summary report
summary_report <- data.frame(
  Metric = c("Stability-selected genes", "Final signature genes",
             "EPV", "Minimum events",
             "Pathways covered", "Drug targets",
             "Literature supported genes", "Highly correlated pairs"),
  Value = c(length(final_genes), length(final_gene_signature),
            sprintf("%.2f", current_epv), min_events,
            nrow(pathway_count), nrow(drug_targets),
            nrow(literature_support), nrow(high_cor_pairs)/2),
  stringsAsFactors = FALSE
)

write_csv(summary_report, file.path(DATA_DERIVED_DIR, "final_signature_summary.csv"))

# ============================================================================
# Step 7: Export additional analysis data
# ============================================================================

cat("\n[Step 7] Exporting additional analysis data...\n")

# Export EPV stats
epv_df <- data.frame(
  Gene_Count = 5:20,
  EPV = min_events / (5:20)
)

write_csv(epv_df, file.path(DATA_DERIVED_DIR, "epv_stats.csv"))

# Export gene correlation matrix for final signature genes
expr_final_signature <- expr_diet[final_gene_signature, ]
final_cor_matrix <- cor(t(expr_final_signature), method = "pearson", use = "pairwise.complete.obs")

saveRDS(final_cor_matrix, file.path(DATA_DERIVED_DIR, "gene_correlation_matrix.rds"))

cat("✓ Additional analysis data exported\n")
cat(sprintf("  - EPV stats: %s/epv_stats.csv\n", DATA_DERIVED_DIR))
cat(sprintf("  - Gene correlation matrix: %s/gene_correlation_matrix.rds\n", DATA_DERIVED_DIR))

# ============================================================================
# Step 8: Lock model weights
# ============================================================================

cat("\n[Step 8] Locking model weights...\n")
cat("WARNING: This is the ONLY legitimate weight file for downstream validation!\n")

# Extract only the final signature genes
final_expr <- expr_diet[final_gene_signature, ]

# Transpose to sample × gene format for glmnet
X_final <- t(final_expr)

# Apply the same scale normalization as stability selection
cat("  Applying log2(X + 1) transformation...\n")
X_final_log <- log2(X_final + 1)

cat("  Applying column-wise standardization...\n")
X_final_scaled <- scale(X_final_log)

# Prepare response variable
y <- ifelse(t2_labels_df$T2_Status == "T2-high", 1, 0)

# Fit final model with cross-validation to select lambda
ALPHA <- 0.7  # ElasticNet mixing parameter
cat("  Running cross-validation to select lambda...\n")
cv_fit <- cv.glmnet(X_final_scaled, y,
                    family = "binomial",
                    alpha = ALPHA,
                    nfolds = 10,
                    type.measure = "auc",
                    standardize = FALSE)  # Manual standardization already applied

# Use lambda.min for final model
lambda_used <- cv_fit$lambda.min
cat(sprintf("  Selected lambda: %.6f\n", lambda_used))

# Fit final model on full dataset with selected lambda
final_fit <- glmnet(X_final_scaled, y,
                    family = "binomial",
                    alpha = ALPHA,
                    lambda = lambda_used,
                    standardize = FALSE)  # Manual standardization already applied

# Extract coefficients
final_coef <- coef(final_fit)

# Create weights dataframe
locked_weights <- data.frame(
  Gene = rownames(final_coef)[-1],  # Remove intercept
  Weight = as.numeric(final_coef)[-1],
  stringsAsFactors = FALSE
) %>%
  arrange(desc(abs(Weight)))

# Save locked weights
write_csv(locked_weights, file.path(DATA_DERIVED_DIR, "locked_weights.csv"))

# Save final gene panel - both with count and fixed interface name
write_csv(
  data.frame(Gene = final_gene_signature, stringsAsFactors = FALSE),
  file.path(DATA_DERIVED_DIR, sprintf("gene_panel_%d.csv", length(final_gene_signature)))
)

# Save with fixed interface name
write_csv(
  data.frame(Gene = final_gene_signature, stringsAsFactors = FALSE),
  file.path(DATA_DERIVED_DIR, "gene_panel_final.csv")
)

cat(sprintf("\n✓ Weights locked and saved: %s/locked_weights.csv\n", DATA_DERIVED_DIR))
cat(sprintf("✓ Final gene panel saved: %s/gene_panel%d.csv\n", DATA_DERIVED_DIR, length(final_gene_signature)))
cat(sprintf("✓ Final gene panel (fixed interface): %s/gene_panel_final.csv\n", DATA_DERIVED_DIR))
cat("This file is the ONLY legitimate weight file for downstream validation!\n")

# ============================================================================
# Final report
# ============================================================================

cat("\n=================================================================\n")
cat("✓ Phase 6 Complete: Final T2 Gene Signature\n")
cat("=================================================================\n")
cat(sprintf("Final signature: %d genes\n", length(final_gene_signature)))
cat("\nFinal signature genes:\n")
print(final_signature_df)
cat("\n=================================================================\n")
cat("→ Next steps: Validate signature in independent cohorts\n")
cat("=================================================================\n\n")

#严厉警告：特征工程已正式锁定，严禁下游篡改
cat("\n=================================================================\n")
cat("WARNING: FEATURE ENGINEERING HAS BEEN OFFICIALLY LOCKED\n")
cat("=================================================================\n")
cat("The final signature genes and locked weights file are the ONLY\n")
cat("legitimate inputs for downstream validation and application.\n")
cat("Any modification to these files will invalidate the model.\n")
cat("=================================================================\n\n")
