# =========================================================================
# Script: 02_limma_filtering.R
# Purpose: Identify differentially expressed genes between T2-high and T2-low
# Author: Zuo
# Date: 2026-02-26
#
# Inputs:
#   - data/processed_diet/GSE152004_diet_expr.rds  # Diet expression matrix
#   - data/derived/GSE152004_T2_labels.csv  # T2 labels
#   - data/raw/candidate_genes_annotated_frozen.csv  # Annotated candidate genes
#
# Outputs:
#   - data/derived/limma_full_stats.csv  # Full differential expression stats
#   - data/derived/limma_passed_genes_primary.csv  # Genes passing primary thresholds
#   - data/derived/limma_passed_genes_fallback_fc1p5_fdr0p01.csv  # Genes passing fallback thresholds
#   - data/derived/limma_thresholds_report.csv  # Thresholds report
#   - output/logs/limma_scale_transformation.log  # Scale transformation log
# =========================================================================

library(limma)
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
DATA_RAW_DIR <- file.path(base_dir, "data", "raw")
DATA_DERIVED_DIR <- file.path(base_dir, "data", "derived")
LOG_DIR <- file.path(base_dir, "output", "logs")

# Create directories if not exist
dir.create(DATA_DERIVED_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

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
cat("Phase 4: Differential Expression Analysis\n")
cat("=================================================================\n\n")
cat(paste0("Base directory: ", base_dir, "\n"))
cat(paste0("Processed diet directory: ", DATA_PROCESSED_DIET_DIR, "\n"))
cat(paste0("Raw directory: ", DATA_RAW_DIR, "\n"))
cat(paste0("Derived directory: ", DATA_DERIVED_DIR, "\n"))
cat(paste0("Logs directory: ", LOG_DIR, "\n\n"))

# ============================================================================
# Step 1: Load data
# ============================================================================

cat("[Step 1] Loading data...\n")

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
candidates <- read_csv(file.path(DATA_RAW_DIR, "candidate_genes_annotated_frozen.csv"))

# Extract candidate genes
candidate_genes <- candidates$Gene
candidates_in_data <- candidate_genes[candidate_genes %in% rownames(expr_diet)]

cat(sprintf("  Candidate genes available: %d/%d\n", length(candidates_in_data), length(candidate_genes)))

expr_candidates <- expr_diet[candidates_in_data, ]

# ============================================================================
# Step 1.1: Apply log2 semantic safety lock
# ============================================================================

cat("\n[Step 1.1] Checking expression scale and applying log2 transformation if needed...\n")

# Check current scale
q99_before <- quantile(expr_candidates, 0.99, na.rm=TRUE)
max_before <- max(expr_candidates, na.rm=TRUE)

cat(sprintf("  Before transformation - 99th percentile: %.2f, Max: %.2f\n", q99_before, max_before))

# Apply log2 transformation if needed (if data appears to be counts-scale)
log_transformed <- FALSE
if (q99_before > 100) {
  cat("  Applying log2(x + 1) transformation...\n")
  expr_candidates <- log2(expr_candidates + 1)
  log_transformed <- TRUE
}

# Check transformed scale
q99_after <- quantile(expr_candidates, 0.99, na.rm=TRUE)
max_after <- max(expr_candidates, na.rm=TRUE)

cat(sprintf("  After transformation - 99th percentile: %.2f, Max: %.2f\n", q99_after, max_after))
cat(sprintf("  Log2 transformation applied: %s\n", ifelse(log_transformed, "Yes", "No")))

# Save transformation info to log
log_file <- file.path(LOG_DIR, "limma_scale_transformation.log")
sink(log_file, append=TRUE)
cat("\n", as.character(Sys.time()), "\n")
cat("Limma scale transformation log:\n")
cat(sprintf("  Before transformation - 99th percentile: %.2f, Max: %.2f\n", q99_before, max_before))
cat(sprintf("  After transformation - 99th percentile: %.2f, Max: %.2f\n", q99_after, max_after))
cat(sprintf("  Log2 transformation applied: %s\n", ifelse(log_transformed, "Yes", "No")))
sink()

# ============================================================================
# Step 2: Differential expression with limma
# ============================================================================

cat("\n[Step 2] Running limma differential expression...\n")

# Design matrix
t2_status <- factor(t2_labels_df$T2_Status, levels = c("T2-low", "T2-high"))
design <- model.matrix(~ t2_status)

# Fit linear model (already in gene x sample format for limma)
fit <- lmFit(expr_candidates, design)
fit <- eBayes(fit)

# Extract results
deg_results <- topTable(fit, coef = 2, number = Inf, sort.by = "P")
deg_results$Gene <- rownames(deg_results)

cat(sprintf("✓ Differential expression complete: %d genes tested\n",
            nrow(deg_results)))

# ============================================================================
# Step 3: Apply filtering criteria with primary and fallback thresholds
# ============================================================================

cat("\n[Step 3] Applying filters...\n")

# Primary thresholds (locked)
LOGFC_THRESHOLD_PRIMARY <- 1.0  # FC > 2
FDR_THRESHOLD_PRIMARY <- 0.001

# Fallback thresholds (locked)
LOGFC_THRESHOLD_FALLBACK <- log2(1.5)  # FC > 1.5
FDR_THRESHOLD_FALLBACK <- 0.01

# Apply primary thresholds
deg_pass_primary <- deg_results %>%
  filter(abs(logFC) > LOGFC_THRESHOLD_PRIMARY & adj.P.Val < FDR_THRESHOLD_PRIMARY)

# Apply fallback thresholds
deg_pass_fallback <- deg_results %>%
  filter(abs(logFC) > LOGFC_THRESHOLD_FALLBACK & adj.P.Val < FDR_THRESHOLD_FALLBACK)

# Determine if fallback is triggered
fallback_triggered <- nrow(deg_pass_primary) < 20

# Select which set to use for downstream analysis
if (fallback_triggered) {
  deg_pass <- deg_pass_fallback
  cat("  Fallback triggered: using fallback thresholds due to low primary pass count\n")
} else {
  deg_pass <- deg_pass_primary
  cat("  Using primary thresholds\n")
}

cat(sprintf("\nPrimary thresholds (|logFC|>%.2f, FDR<%.3f):\n",
            LOGFC_THRESHOLD_PRIMARY, FDR_THRESHOLD_PRIMARY))
cat(sprintf("  - Passed: %d genes (%.1f%%)\n",
            nrow(deg_pass_primary), 100 * nrow(deg_pass_primary) / nrow(deg_results)))

cat(sprintf("\nFallback thresholds (|logFC|>%.2f, FDR<%.3f):\n",
            LOGFC_THRESHOLD_FALLBACK, FDR_THRESHOLD_FALLBACK))
cat(sprintf("  - Passed: %d genes (%.1f%%)\n",
            nrow(deg_pass_fallback), 100 * nrow(deg_pass_fallback) / nrow(deg_results)))

cat(sprintf("\nSelected for downstream: %d genes (using %s)\n",
            nrow(deg_pass), ifelse(fallback_triggered, "fallback thresholds", "primary thresholds")))

# ============================================================================
# Step 4: Export results
# ============================================================================

cat("\n[Step 4] Exporting results...\n")

# Export full stats
write_csv(deg_results, file.path(DATA_DERIVED_DIR, "limma_full_stats.csv"))

# Export primary results
write_csv(deg_pass_primary, file.path(DATA_DERIVED_DIR, "limma_passed_genes_primary.csv"))

# Export fallback results
write_csv(deg_pass_fallback, file.path(DATA_DERIVED_DIR, "limma_passed_genes_fallback_fc1p5_fdr0p01.csv"))

# Export thresholds report
thresholds_report <- data.frame(
  fallback_triggered = fallback_triggered,
  n_pass_primary = nrow(deg_pass_primary),
  n_pass_fallback = nrow(deg_pass_fallback),
  logfc_threshold_primary = LOGFC_THRESHOLD_PRIMARY,
  fdr_threshold_primary = FDR_THRESHOLD_PRIMARY,
  logfc_threshold_fallback = LOGFC_THRESHOLD_FALLBACK,
  fdr_threshold_fallback = FDR_THRESHOLD_FALLBACK,
  q99_before_log2 = q99_before,
  max_before_log2 = max_before,
  q99_after_log2 = q99_after,
  max_after_log2 = max_after,
  log2_applied = log_transformed,
  stringsAsFactors = FALSE
)

write_csv(thresholds_report, file.path(DATA_DERIVED_DIR, "limma_thresholds_report.csv"))

cat(sprintf("✓ Results exported:\n"))
cat(sprintf("  - Full stats: %s/limma_full_stats.csv\n", DATA_DERIVED_DIR))
cat(sprintf("  - Primary passed genes: %s/limma_passed_genes_primary.csv\n", DATA_DERIVED_DIR))
cat(sprintf("  - Fallback passed genes: %s/limma_passed_genes_fallback_fc1p5_fdr0p01.csv\n", DATA_DERIVED_DIR))
cat(sprintf("  - Thresholds report: %s/limma_thresholds_report.csv\n", DATA_DERIVED_DIR))

cat("\n=================================================================\n")
cat(sprintf("✓ Phase 4 Complete: %d genes passed filters\n", nrow(deg_pass)))
cat("→ Next: Stability selection\n")
cat("=================================================================\n\n")
