# =========================================================================
# Script: 01_prep_discovery_GSE152004.R
# Purpose: Extract minimal gene set for T2 asthma signature discovery
# Author: Zuo
# Date: 2026-02-26
#
# Inputs:
#   - data/raw/gse152004&GPL11154/GSE152004_695_expr_norm.txt.gz  # Normalized expression data
#   - data/raw/candidate_genes_annotated_frozen.csv  # Candidate genes
#
# Outputs:
#   - data/processed_diet/GSE152004_diet_expr.rds  # Minimal expression matrix
#   - data/processed_diet/GSE152004_pheno.rds  # Sample information
#   - output/logs/GSE152004_expr_norm_fingerprint.tsv  # Expression fingerprint
#   - output/logs/GSE152004_diet_fingerprint.tsv  # Diet fingerprint
# =========================================================================

library(dplyr)
library(readr)
library(data.table)
library(org.Hs.eg.db)

# Set paths
# Get base directory (project root)
base_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# Verify we're running from project root
if (!all(dir.exists(file.path(base_dir, c("data", "analysis", "data_preparation"))))) {
  stop("Please run this script from the project root directory (the folder containing data/, analysis/, data_preparation/).", call. = FALSE)
}

# Define directories
DATA_RAW_DIR <- file.path(base_dir, "data", "raw")
DATA_PROCESSED_FULL_DIR <- file.path(base_dir, "data", "processed_full")
DATA_PROCESSED_DIET_DIR <- file.path(base_dir, "data", "processed_diet")
LOG_DIR <- file.path(base_dir, "output", "logs")

# Create directories if not exist
dir.create(DATA_PROCESSED_FULL_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_PROCESSED_DIET_DIR, recursive = TRUE, showWarnings = FALSE)
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
cat("Phase 2: GSE152004 Data Preparation - The Data Diet\n")
cat("=================================================================\n\n")
cat(paste0("Base directory: ", base_dir, "\n"))
cat(paste0("Raw directory: ", DATA_RAW_DIR, "\n"))
cat(paste0("Processed full directory: ", DATA_PROCESSED_FULL_DIR, "\n"))
cat(paste0("Processed diet directory: ", DATA_PROCESSED_DIET_DIR, "\n"))
cat(paste0("Logs directory: ", LOG_DIR, "\n\n"))

# ============================================================================
# Step 1: Load GSE152004 from local files
# ============================================================================

cat("[Step 1] Loading GSE152004 from local files...\n")

# Load expression matrix from normalized expression file
expr_file <- normalizePath(file.path(DATA_RAW_DIR, "gse152004&GPL11154", "GSE152004_695_expr_norm.txt.gz"), winslash="/", mustWork=TRUE)

# Record input file information
input_file_path <- expr_file
input_file_basename <- basename(input_file_path)
input_file_mtime <- file.info(input_file_path)$mtime
input_file_size <- file.size(input_file_path)

# Print input file information
cat("  Input file information:", "\n")
cat(sprintf("  - File path: %s\n", input_file_path))
cat(sprintf("  - File basename: %s\n", input_file_basename))
cat(sprintf("  - Last modified: %s\n", input_file_mtime))
cat(sprintf("  - File size: %d bytes\n", input_file_size))

# Specialized reading logic for GSE152004_695_expr_norm.txt.gz
cat("  Using specialized reading logic for expr_norm file...\n")

# Step 1: Read first line as samples
con <- file(expr_file, "r")
first_line <- readLines(con, n = 1)
close(con)

# Split first line and check structure
first_line_split <- unlist(strsplit(first_line, "\t"))
cat(sprintf("  - First line elements: %d\n", length(first_line_split)))
cat(sprintf("  - First 5 elements: %s\n", paste(head(first_line_split, 5), collapse=", ")))

# Step 2: Read remaining data with fread
cat("  Reading data with fread...\n")
data <- fread(expr_file, header = FALSE, skip = 1, sep = "\t")

# Check data dimensions
cat(sprintf("  - Data dimensions: %d rows × %d columns\n", nrow(data), ncol(data)))

# Step 3: Manually set column names
# First column is Gene, rest are samples from first line
colnames(data) <- c("Gene", first_line_split)
samples <- first_line_split

cat(sprintf("  - Number of samples: %d\n", length(samples)))

# Step 4: Set Gene column as row names and convert to numeric matrix
cat("  Converting to numeric matrix...\n")
gene_symbols <- data$Gene
data$Gene <- NULL

# Debug: Check data dimensions after removing Gene column
cat(sprintf("  - Data dimensions after removing Gene: %d rows × %d columns\n", nrow(data), ncol(data)))
cat(sprintf("  - Number of samples: %d\n", length(samples)))

# Convert to matrix
mat <- as.matrix(data)
rownames(mat) <- gene_symbols
colnames(mat) <- samples

# Convert to numeric
mat <- apply(mat, 2, as.numeric)
rownames(mat) <- gene_symbols
colnames(mat) <- samples

# Debug: Check matrix dimensions
cat(sprintf("  - Matrix dimensions: %d rows × %d columns\n", nrow(mat), ncol(mat)))

# Step 5: Assertions
cat("  Running assertions...\n")
# Print first few row names for debugging
cat("  First 5 row names: ", paste(head(rownames(mat), 5), collapse=", "), "\n")

# Use stop for assertions to ensure data integrity
if (ncol(mat) != 695) {
  stop("Number of columns must be 695")
}
cat("  ✓ Number of columns is 695\n")

if (!all(colnames(mat) == samples)) {
  stop("Column names must match samples")
}
cat("  ✓ Column names match samples\n")

# Check if row names look like gene symbols (allowing dots and hyphens)
if (!all(grepl("^[A-Z0-9.-]+$", head(rownames(mat), 3)))) {
  stop("First 3 row names must be gene symbols")
}
cat("  ✓ First 3 row names are gene symbols\n")

cat(sprintf("  - Final matrix dimensions: %d genes × %d samples\n", nrow(mat), ncol(mat)))

# Step 6: Print and save scale fingerprint
cat("  Generating scale fingerprint...\n")
mat_values <- as.vector(mat)
mat_quantiles <- quantile(mat_values, c(0.5, 0.9, 0.99, 1))

cat("  Scale fingerprint:", "\n")
cat(sprintf("  - p50: %.2f\n", mat_quantiles[1]))
cat(sprintf("  - p90: %.2f\n", mat_quantiles[2]))
cat(sprintf("  - p99: %.2f\n", mat_quantiles[3]))
cat(sprintf("  - max: %.2f\n", mat_quantiles[4]))

# Save fingerprint
fingerprint_file <- file.path(LOG_DIR, "GSE152004_expr_norm_fingerprint.tsv")
fingerprint_data <- data.frame(
  Metric = c("p50", "p90", "p99", "max"),
  Value = as.numeric(mat_quantiles)
)
write_tsv(fingerprint_data, fingerprint_file)
cat(sprintf("  Fingerprint saved to: %s\n", fingerprint_file))

# Set expr_raw as the transposed matrix (samples x genes format)
expr_raw <- t(mat)

# Debug: Check expr_raw structure
cat(sprintf("  Debug: expr_raw dimensions: %d samples × %d genes\n", nrow(expr_raw), ncol(expr_raw)))
cat(sprintf("  Debug: First few row names: %s\n", paste(head(rownames(expr_raw)), collapse=", ")))
cat(sprintf("  Debug: First few column names: %s\n", paste(head(colnames(expr_raw)), collapse=", ")))

# Create basic phenotype data frame
pheno <- data.frame(
  Sample_ID = rownames(expr_raw),
  stringsAsFactors = FALSE
)
rownames(pheno) <- rownames(expr_raw)

cat(sprintf("✓ Loaded: %d samples × %d genes\n",
            nrow(expr_raw), ncol(expr_raw)))

# ============================================================================
# Step 2: Quality Control
# ============================================================================

cat("\n[Step 2] Quality control...\n")

# QC Criterion 1: Remove samples with >20% missing values
missing_rate <- rowMeans(is.na(expr_raw))
samples_qc1 <- missing_rate < 0.20

cat(sprintf("  QC1: Missing values - %d samples removed\n",
            sum(!samples_qc1)))

# QC Criterion 2: Remove technical replicates (keep first)
# Assuming sample names don't have obvious replicate patterns
# If replicates exist, handle here

# QC Criterion 3: Remove outliers by PCA (>3SD)
expr_for_pca <- expr_raw[samples_qc1, ]
expr_for_pca[is.na(expr_for_pca)] <- rowMeans(expr_for_pca, na.rm = TRUE)

# Remove genes with zero variance before PCA
gene_var <- apply(expr_for_pca, 2, var, na.rm = TRUE)
non_zero_var_genes <- gene_var > 1e-10
expr_for_pca <- expr_for_pca[, non_zero_var_genes]

cat(sprintf("  QC2.1: Removed %d genes with zero variance\n", sum(!non_zero_var_genes)))

# Apply QC filters
samples_pass_qc <- samples_qc1

# Apply filters
expr_qc <- expr_raw[samples_pass_qc, ]
pheno_qc <- pheno[samples_pass_qc, ]

# Debug: Check dimensions after QC
cat(sprintf("  Debug: After QC - expr_qc dimensions: %d samples × %d genes\n", nrow(expr_qc), ncol(expr_qc)))

cat(sprintf("\n✓ After QC: %d samples retained (%.1f%%)\n",
            nrow(expr_qc),
            100 * nrow(expr_qc) / nrow(expr_raw)))

# ============================================================================
# Step 3: Data Diet - Extract only necessary genes
# ============================================================================

cat("\n[Step 3] Applying Data Diet...\n")

# Load candidate genes from frozen annotation
candidates <- read_csv(file.path(DATA_RAW_DIR, "candidate_genes_annotated_frozen.csv"))
candidate_genes <- candidates$Gene

# Woodruff 2009 12-gene T2 signature
woodruff_genes <- c("CST1", "CST2", "CLCA1", "POSTN", "SERPINB2",
                   "SERPINB10", "CEACAM5", "DNAH5", "ALOX15",
                   "LRRC31", "SLC26A4", "AGER")

# Calculate union of both gene sets
target_genes <- unique(c(candidate_genes, woodruff_genes))

# Check which target genes are available in the data
available_genes <- target_genes[target_genes %in% colnames(expr_qc)]

cat(sprintf("  Target genes: %d (candidates: %d + Woodruff: %d)\n",
            length(target_genes), length(candidate_genes), length(woodruff_genes)))
cat(sprintf("  Available in data: %d/ %d\n", length(available_genes), length(target_genes)))

# Extract only the target genes
expr_diet <- expr_qc[, available_genes]

# ============================================================================
# Step 4: Ensure gene × sample format
# ============================================================================

cat("\n[Step 4] Ensuring gene × sample format...\n")

# Check current format
row_id_type <- determine_id_type(rownames(expr_diet))
col_id_type <- determine_id_type(colnames(expr_diet))

cat(sprintf("  Current row IDs: %s\n", row_id_type))
cat(sprintf("  Current column IDs: %s\n", col_id_type))

# Check if transpose is needed
if (row_id_type != "SYMBOL" && col_id_type == "SYMBOL") {
  cat("  Transposing matrix to gene × sample format...\n")
  expr_diet <- t(expr_diet)
  
  # Update IDs after transpose
  row_id_type <- determine_id_type(rownames(expr_diet))
  col_id_type <- determine_id_type(colnames(expr_diet))
  
  cat(sprintf("  After transpose - row IDs: %s\n", row_id_type))
  cat(sprintf("  After transpose - column IDs: %s\n", col_id_type))
}

# Assertions
stopifnot(determine_id_type(rownames(expr_diet)) == "SYMBOL")
stopifnot(ncol(expr_diet) == 695)  # Should match full sample count

# Check for missing genes
missing_genes <- setdiff(target_genes, rownames(expr_diet))
if (length(missing_genes) > 0) {
  cat(sprintf("  Warning: %d target genes not found in data:\n", length(missing_genes)))
  cat(sprintf("  %s\n", paste(missing_genes, collapse=", ")))
}

# Debug: Check dimensions after data diet
cat(sprintf("  Debug: After Data Diet - expr_diet dimensions: %d genes × %d samples\n",
            nrow(expr_diet), ncol(expr_diet)))

# Check data scale
cat("\n[Step 4.5] Checking data scale...\n")
diet_quantiles <- quantile(expr_diet, c(0.5, 0.9, 0.99, 1))
cat("  Data diet quantiles (50%, 90%, 99%, 100%):\n")
cat(sprintf("  50%%: %.2f\n", diet_quantiles[1]))
cat(sprintf("  90%%: %.2f\n", diet_quantiles[2]))
cat(sprintf("  99%%: %.2f\n", diet_quantiles[3]))
cat(sprintf("  100%%: %.2f\n", diet_quantiles[4]))
cat("  Data diet value range: min=", min(expr_diet), "max=", max(expr_diet), "\n")

# Generate scale fingerprint before saving
cat("\n[Step 4.6] Generating scale fingerprint...\n")
diet_dim <- dim(expr_diet)
diet_quantiles <- quantile(as.vector(expr_diet), c(0.5, 0.9, 0.99, 1))
diet_zero_ratio <- mean(expr_diet == 0)
diet_range <- range(expr_diet)

# Print fingerprint
cat("  Scale fingerprint:", "\n")
cat(sprintf("  - Dimensions: %d x %d\n", diet_dim[1], diet_dim[2]))
cat(sprintf("  - 50th quantile: %.2f\n", diet_quantiles[1]))
cat(sprintf("  - 90th quantile: %.2f\n", diet_quantiles[2]))
cat(sprintf("  - 99th quantile: %.2f\n", diet_quantiles[3]))
cat(sprintf("  - 100th quantile: %.2f\n", diet_quantiles[4]))
cat(sprintf("  - Zero ratio: %.4f\n", diet_zero_ratio))
cat(sprintf("  - Range: %.2f to %.2f\n", diet_range[1], diet_range[2]))

# Save fingerprint to file
fingerprint_file <- file.path(LOG_DIR, "GSE152004_diet_fingerprint.tsv")
fingerprint_data <- data.frame(
  Metric = c("Dimensions", "50th_quantile", "90th_quantile", "99th_quantile", "100th_quantile", "Zero_ratio", "Min", "Max"),
  Value = c(paste(diet_dim[1], diet_dim[2], sep="x"), as.numeric(diet_quantiles), diet_zero_ratio, diet_range[1], diet_range[2])
)
write_tsv(fingerprint_data, fingerprint_file)
cat(sprintf("  Fingerprint saved to: %s\n", fingerprint_file))

# Release memory
rm(expr_raw, expr_qc, expr_for_pca)
gc()

cat("\n✓ Memory released: large matrix removed\n")

# ============================================================================
# Step 5: Export results
# ============================================================================

cat("\n[Step 5] Exporting results...\n")

# Save minimal expression matrix
diet_expr_file <- file.path(DATA_PROCESSED_DIET_DIR, "GSE152004_diet_expr.rds")
saveRDS(expr_diet, diet_expr_file)

# Save sample information
diet_pheno_file <- file.path(DATA_PROCESSED_DIET_DIR, "GSE152004_pheno.rds")
saveRDS(pheno_qc, diet_pheno_file)

# Validate saved file by reading it back and comparing fingerprints
cat("\n[Step 5.1] Validating saved file...\n")
expr_diet_readback <- readRDS(diet_expr_file)

# Generate fingerprint for readback data
readback_dim <- dim(expr_diet_readback)
readback_quantiles <- quantile(as.vector(expr_diet_readback), c(0.5, 0.9, 0.99, 1))
readback_zero_ratio <- mean(expr_diet_readback == 0)
readback_range <- range(expr_diet_readback)

# Print readback fingerprint
cat("  Readback scale fingerprint:", "\n")
cat(sprintf("  - Dimensions: %d x %d\n", readback_dim[1], readback_dim[2]))
cat(sprintf("  - 50th quantile: %.2f\n", readback_quantiles[1]))
cat(sprintf("  - 90th quantile: %.2f\n", readback_quantiles[2]))
cat(sprintf("  - 99th quantile: %.2f\n", readback_quantiles[3]))
cat(sprintf("  - 100th quantile: %.2f\n", readback_quantiles[4]))
cat(sprintf("  - Zero ratio: %.4f\n", readback_zero_ratio))
cat(sprintf("  - Range: %.2f to %.2f\n", readback_range[1], readback_range[2]))

# Compare fingerprints
cat("  Comparing fingerprints...\n")
dim_match <- all(diet_dim == readback_dim)
quantiles_match <- all(abs(diet_quantiles - readback_quantiles) < 1e-6)
zero_ratio_match <- abs(diet_zero_ratio - readback_zero_ratio) < 1e-6
range_match <- all(abs(diet_range - readback_range) < 1e-6)

all_match <- dim_match && quantiles_match && zero_ratio_match && range_match

if (all_match) {
  cat("  ✓ All fingerprints match! Saved file is valid.\n")
} else {
  cat("  ✗ Fingerprints do not match! Saved file may be corrupted.\n")
  if (!dim_match) cat("    - Dimensions do not match\n")
  if (!quantiles_match) cat("    - Quantiles do not match\n")
  if (!zero_ratio_match) cat("    - Zero ratio does not match\n")
  if (!range_match) cat("    - Range does not match\n")
}

cat("\n✓ Results exported:\n")
cat(sprintf("  - Minimal expression matrix: %s\n", diet_expr_file))
cat(sprintf("  - Sample information: %s\n", diet_pheno_file))

cat("\n=================================================================\n")
cat("✓ Phase 2 Complete: Discovery data prepared with Data Diet\n")
cat("→ Next: T2 labeling and classification\n")
cat("=================================================================\n\n")
