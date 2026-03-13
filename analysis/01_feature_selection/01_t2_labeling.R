# =========================================================================
# Script: 01_t2_labeling.R
# Purpose: Calculate T2 signature and classify samples into T2-high/low
# Author: Zuo
# Date: 2026-02-26
#
# Inputs:
#   - data/processed_diet/GSE152004_diet_expr.rds  # Diet expression matrix
#
# Outputs:
#   - data/derived/GSE152004_T2_labels.csv  # T2 labels and signature scores
# =========================================================================

library(dplyr)
library(readr)
library(org.Hs.eg.db)

# Set paths
# Get base directory (project root)
base_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# Verify we're running from project root
if (!all(dir.exists(file.path(base_dir, c("data", "analysis", "data_preparation"))))) {
  stop("Please run this script from the project root directory (the folder containing data/, analysis/, data_preparation/).", call. = FALSE)
}

# Define directories
DATA_PROCESSED_DIET_DIR <- file.path(base_dir, "data", "processed_diet")
DATA_DERIVED_DIR <- file.path(base_dir, "data", "derived")

# Create directory if not exist
dir.create(DATA_DERIVED_DIR, recursive = TRUE, showWarnings = FALSE)

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

# Log base directory and paths
cat("=================================================================\n")
cat("Phase 3: T2 Labeling for Discovery Cohort\n")
cat("=================================================================\n\n")
cat(paste0("Base directory: ", base_dir, "\n"))
cat(paste0("Processed diet directory: ", DATA_PROCESSED_DIET_DIR, "\n"))
cat(paste0("Derived directory: ", DATA_DERIVED_DIR, "\n\n"))

# ============================================================================
# Step 1: Load diet expression matrix
# ============================================================================

cat("[Step 1] Loading diet expression matrix...\n")

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

cat(sprintf("✓ Loaded diet expression matrix: %d genes × %d samples\n",
            nrow(expr_diet), ncol(expr_diet)))

# ============================================================================
# Step 2: T2 Classification using Woodruff 12-gene signature
# ============================================================================

cat("\n[Step 2] T2-high/T2-low classification...\n")

# Woodruff 2009 12-gene T2 signature
woodruff_genes <- c("CST1", "CST2", "CLCA1", "POSTN", "SERPINB2",
                   "SERPINB10", "CEACAM5", "DNAH5", "ALOX15",
                   "LRRC31", "SLC26A4", "AGER")

# Check gene availability
woodruff_genes_available <- woodruff_genes %in% rownames(expr_diet)
cat(sprintf("  T2 genes available: %d/12\n", sum(woodruff_genes_available)))

if (sum(woodruff_genes_available) < 8) {
  stop("ERROR: <8 T2 signature genes available. Cannot proceed.")
}

# Extract T2 gene expression
available_gene_names <- woodruff_genes[woodruff_genes %in% rownames(expr_diet)]
cat(sprintf("  Available T2 genes: %s\n", paste(available_gene_names, collapse=", ")))

# Extract available T2 gene expression and transpose to sample × gene for calculation
t2_expr <- t(expr_diet[available_gene_names, , drop=FALSE])

# Debug: Check dimensions
cat(sprintf("  Debug: expr_diet dimensions: %d genes × %d samples\n", nrow(expr_diet), ncol(expr_diet)))
cat(sprintf("  Debug: t2_expr dimensions: %d samples × %d genes\n", nrow(t2_expr), ncol(t2_expr)))

# Calculate T2 signature (mean Z-score)
t2_signature <- rowMeans(scale(t2_expr))

# Debug: Check signature length
cat(sprintf("  Debug: t2_signature length: %d\n", length(t2_signature)))

# Median-based classification
t2_labels <- ifelse(t2_signature > median(t2_signature),
                    "T2-high", "T2-low")

# Debug: Check labels length
cat(sprintf("  Debug: t2_labels length: %d\n", length(t2_labels)))

# Summary
t2_summary <- table(t2_labels)
cat(sprintf("\n✓ T2 classification complete:\n"))
cat(sprintf("  - T2-high: %d samples (%.1f%%)\n",
            t2_summary["T2-high"],
            100 * t2_summary["T2-high"] / sum(t2_summary)))
cat(sprintf("  - T2-low:  %d samples (%.1f%%)\n",
            t2_summary["T2-low"],
            100 * t2_summary["T2-low"] / sum(t2_summary)))

# Statistical test for group separation
t_test <- t.test(t2_signature ~ t2_labels)
cat(sprintf("\n  Group separation: t=%.2f, p=%.2e\n",
            t_test$statistic, t_test$p.value))

# ============================================================================
# Step 3: Export T2 labels
# ============================================================================

cat("\n[Step 3] Exporting T2 labels...\n")

# Create labels data frame
t2_labels_df <- data.frame(
  Sample_ID = colnames(expr_diet),
  T2_Status = t2_labels,
  T2_Signature_Score = t2_signature,
  stringsAsFactors = FALSE
)

# Save labels
write_csv(t2_labels_df, file.path(DATA_DERIVED_DIR, "GSE152004_T2_labels.csv"))

cat(sprintf("✓ T2 labels exported: %s/GSE152004_T2_labels.csv\n", DATA_DERIVED_DIR))

cat("\n=================================================================\n")
cat("✓ Phase 3 Complete: T2 labeling finished\n")
cat("→ Next: Differential expression analysis\n")
cat("=================================================================\n\n")
