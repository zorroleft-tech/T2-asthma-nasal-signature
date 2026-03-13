# =========================================================================
# Script: 03_stability_selection.R
# Purpose: Select robust genes and lock final model weights
# Author: Zuo
# Date: 2026-02-26
#
# Inputs:
#   - data/processed_diet/GSE152004_diet_expr.rds  # Diet expression matrix
#   - data/derived/GSE152004_T2_labels.csv  # T2 labels
#   - data/derived/limma_thresholds_report.csv  # Thresholds report
#   - data/derived/limma_passed_genes_primary.csv  # Genes passing primary thresholds
#   - data/derived/limma_passed_genes_fallback_fc1p5_fdr0p01.csv  # Genes passing fallback thresholds
#
# Outputs:
#   - data/derived/stability_selection_results.rds  # Stability selection results
#   - data/derived/stability_selected_genes.csv  # Selected genes
#   - data/derived/bootstrap_frequencies.csv  # Bootstrap frequencies
#   - data/derived/stability_comparison_summary.csv  # Comparison summary
#   - data/derived/stability_selection_results_final.rds  # Final results
#   - output/logs/stability_audit_GSE152004.tsv  # Audit file
# =========================================================================

library(glmnet)
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
cat("Phase 5: Stability Selection\n")
cat("=================================================================\n\n")
cat(paste0("Base directory: ", base_dir, "\n"))
cat(paste0("Processed diet directory: ", DATA_PROCESSED_DIET_DIR, "\n"))
cat(paste0("Derived directory: ", DATA_DERIVED_DIR, "\n"))
cat(paste0("Logs directory: ", LOG_DIR, "\n\n"))

# ============================================================================
# Step 1: Load data
# ============================================================================

cat("[Step 1] Loading data...\n")

expr_diet <- readRDS(file.path(DATA_PROCESSED_DIET_DIR, "GSE152004_diet_expr.rds"))

# Ensure gene × sample format and validate IDs
# Check if rows are gene symbols and columns are sample IDs
row_is_gene <- all(grepl("^[A-Z0-9.-]+", rownames(expr_diet)))
col_is_sample <- all(grepl("^(HR|SJ)", colnames(expr_diet)))

if (!row_is_gene && col_is_sample) {
  cat("  Transposing matrix to gene × sample format...\n")
  expr_diet <- t(expr_diet)
  
  # Recheck after transpose
  row_is_gene <- all(grepl("^[A-Z0-9.-]+", rownames(expr_diet)))
  col_is_sample <- all(grepl("^(HR|SJ)", colnames(expr_diet)))
}

# Assertions
if (!row_is_gene) {
  stop("Row names must be gene symbols")
}
if (!col_is_sample) {
  stop("Column names must be sample IDs")
}
if (ncol(expr_diet) != 695) {
  stop("Number of columns must be 695")
}  # Should match full sample count

t2_labels_df <- read_csv(file.path(DATA_DERIVED_DIR, "GSE152004_T2_labels.csv"))

# Read thresholds report to determine which passed genes to use
thresholds_report <- read_csv(file.path(DATA_DERIVED_DIR, "limma_thresholds_report.csv"))
fallback_triggered <- thresholds_report$fallback_triggered[1]

# Allow override for testing (set to "primary" or "fallback" to force)
override_type <- NULL  # Set to "primary" or "fallback" to force

# Run comparison mode (run both primary and fallback for comparison)
RUN_COMPARISON <- Sys.getenv("RUN_COMPARISON", "FALSE") == "TRUE"

# Read appropriate passed genes file
if (!is.null(override_type)) {
  if (override_type == "primary") {
    cat("  Override: using primary passed genes\n")
    deg_pass <- read_csv(file.path(DATA_DERIVED_DIR, "limma_passed_genes_primary.csv"))
    passed_genes_type <- "primary"
  } else if (override_type == "fallback") {
    cat("  Override: using fallback passed genes\n")
    deg_pass <- read_csv(file.path(DATA_DERIVED_DIR, "limma_passed_genes_fallback_fc1p5_fdr0p01.csv"))
    passed_genes_type <- "fallback"
  } else {
    stop("Invalid override_type. Must be 'primary' or 'fallback'")
  }
} else if (fallback_triggered) {
  cat("  Fallback triggered: using fallback passed genes\n")
  deg_pass <- read_csv(file.path(DATA_DERIVED_DIR, "limma_passed_genes_fallback_fc1p5_fdr0p01.csv"))
  passed_genes_type <- "fallback"
} else {
  cat("  Using primary passed genes\n")
  deg_pass <- read_csv(file.path(DATA_DERIVED_DIR, "limma_passed_genes_primary.csv"))
  passed_genes_type <- "primary"
}

# Extract DEG-passed genes
deg_genes <- deg_pass$Gene

# Ensure all deg_genes are present in expr_diet
deg_genes_present <- intersect(deg_genes, rownames(expr_diet))
if (length(deg_genes_present) != length(deg_genes)) {
  stop("Some DEG genes not found in expression data")
}

expr_deg <- expr_diet[deg_genes_present, , drop=FALSE]

# Store for later comparison
n_features_in_X <- length(deg_genes_present)
cat(sprintf("  Using %s passed genes: %d features\n", passed_genes_type, n_features_in_X))

# Transpose to get samples × genes format for glmnet
X <- t(expr_deg)

# Prepare response variable
y <- ifelse(t2_labels_df$T2_Status == "T2-high", 1, 0)

# ============================================================================
# Step 1.4: Apply universal scale normalization
# ============================================================================

cat("\n[Step 1.4] Applying universal scale normalization...\n")

# Store original X for comparison
X_original <- X

# Step 1: Log transformation (log2(X + 1))
cat("  Applying log2(X + 1) transformation...\n")
X1 <- log2(X + 1)

# Step 2: Scale by column (mean 0, variance 1)
cat("  Applying column-wise standardization...\n")
X2 <- scale(X1)

# Use the normalized matrix for analysis
X <- X2

# Print dimensions for validation
cat(sprintf("  dim(expr_diet) = %d×%d\n", nrow(expr_diet), ncol(expr_diet)))
cat(sprintf("  dim(X) = %d×%d\n", nrow(X), ncol(X)))
cat(sprintf("✓ %d genes ready for stability selection\n", ncol(X)))

# ============================================================================
# Stability Selection Parameters (PRE-DEFINED)
# ============================================================================

BOOTSTRAP_ITER <- 1000
TRAIN_FRAC <- 0.80
ALPHA <- 0.7  # ElasticNet mixing parameter
SELECTION_THRESHOLD <- 0.70  # 70% selection frequency

cat(sprintf("\nStability Selection Parameters:\n"))
cat(sprintf("  - Bootstrap iterations: %d\n", BOOTSTRAP_ITER))
cat(sprintf("  - Training fraction: %.0f%%\n", TRAIN_FRAC * 100))
cat(sprintf("  - ElasticNet alpha: %.1f\n", ALPHA))
cat(sprintf("  - Selection threshold: %.0f%%\n", SELECTION_THRESHOLD * 100))

# ============================================================================
# Input Audit for Reproducibility
# ============================================================================

cat("\n[Input Audit]\n")

# Create logs directory if it doesn't exist
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

# Prepare audit data
audit_data <- list()

# 1) X/y basic information
audit_data$dim_expr_diet <- paste(nrow(expr_diet), ncol(expr_diet), sep="x")
audit_data$dim_X <- paste(nrow(X), ncol(X), sep="x")
audit_data$y_table <- table(y)
audit_data$first_5_genes <- paste(head(colnames(X), 5), collapse=", ")

# 2) Standardization and value ranges - Original X
audit_data$x_original_range <- paste(range(X_original, na.rm=TRUE), collapse=" to ")
audit_data$x_original_mean <- mean(X_original, na.rm=TRUE)
audit_data$x_original_sd <- sd(X_original, na.rm=TRUE)
audit_data$x_original_na_prop <- mean(is.na(X_original))
audit_data$x_original_p50 <- quantile(X_original, 0.5, na.rm=TRUE)
audit_data$x_original_p90 <- quantile(X_original, 0.9, na.rm=TRUE)
audit_data$x_original_p99 <- quantile(X_original, 0.99, na.rm=TRUE)
audit_data$x_original_max <- max(X_original, na.rm=TRUE)

# 3) Standardization and value ranges - After log2 transformation
audit_data$x_log_range <- paste(range(X1, na.rm=TRUE), collapse=" to ")
audit_data$x_log_mean <- mean(X1, na.rm=TRUE)
audit_data$x_log_sd <- sd(X1, na.rm=TRUE)
audit_data$x_log_na_prop <- mean(is.na(X1))
audit_data$x_log_p50 <- quantile(X1, 0.5, na.rm=TRUE)
audit_data$x_log_p90 <- quantile(X1, 0.9, na.rm=TRUE)
audit_data$x_log_p99 <- quantile(X1, 0.99, na.rm=TRUE)
audit_data$x_log_max <- max(X1, na.rm=TRUE)

# 4) Standardization and value ranges - After scaling
audit_data$x_scaled_range <- paste(range(X, na.rm=TRUE), collapse=" to ")
audit_data$x_scaled_mean <- mean(X, na.rm=TRUE)
audit_data$x_scaled_sd <- sd(X, na.rm=TRUE)
audit_data$x_scaled_na_prop <- mean(is.na(X))
audit_data$x_scaled_p50 <- quantile(X, 0.5, na.rm=TRUE)
audit_data$x_scaled_p90 <- quantile(X, 0.9, na.rm=TRUE)
audit_data$x_scaled_p99 <- quantile(X, 0.99, na.rm=TRUE)
audit_data$x_scaled_max <- max(X, na.rm=TRUE)

audit_data$standardization_info <- "Applied log2(X + 1) transformation followed by column-wise standardization (mean=0, sd=1)"

# 5) Glmnet standardization setting
audit_data$glmnet_standardize <- "FALSE (manual standardization applied)"

# 3) Glmnet key parameters
audit_data$alpha <- ALPHA
audit_data$family <- "binomial"
audit_data$standardize <- audit_data$glmnet_standardize
audit_data$lambda_selection <- "lambda.min"
audit_data$use_cv_glmnet <- "TRUE"

# 4) Reproducibility
audit_data$seed_value <- 2026
audit_data$seed_position <- "Before bootstrap loop"

# Write audit to file
audit_file <- file.path(LOG_DIR, "stability_audit_GSE152004.tsv")
audit_df <- data.frame(
  Parameter = names(audit_data),
  Value = sapply(audit_data, function(x) ifelse(is.table(x), paste(names(x), x, collapse=", "), as.character(x)))
)
write.table(audit_df, audit_file, sep="\t", row.names=FALSE, quote=FALSE)

# Print audit summary
cat("  Audit summary:")
cat(sprintf("\n  - Expression matrix dimensions: %s", audit_data$dim_expr_diet))
cat(sprintf("\n  - Design matrix dimensions: %s", audit_data$dim_X))
cat(sprintf("\n  - y distribution: %s", paste(names(audit_data$y_table), audit_data$y_table, collapse=", ")))
cat(sprintf("\n  - First 5 features: %s", audit_data$first_5_genes))

# Original X statistics
cat("\n  Original X:")
cat(sprintf("\n    - Range: %s", audit_data$x_original_range))
cat(sprintf("\n    - Mean: %.4f, SD: %.4f", audit_data$x_original_mean, audit_data$x_original_sd))
cat(sprintf("\n    - p50: %.4f, p90: %.4f, p99: %.4f, Max: %.4f", 
            audit_data$x_original_p50, audit_data$x_original_p90, 
            audit_data$x_original_p99, audit_data$x_original_max))
cat(sprintf("\n    - NA proportion: %.4f", audit_data$x_original_na_prop))

# After log2 transformation
cat("\n  After log2(X + 1):")
cat(sprintf("\n    - Range: %s", audit_data$x_log_range))
cat(sprintf("\n    - Mean: %.4f, SD: %.4f", audit_data$x_log_mean, audit_data$x_log_sd))
cat(sprintf("\n    - p50: %.4f, p90: %.4f, p99: %.4f, Max: %.4f", 
            audit_data$x_log_p50, audit_data$x_log_p90, 
            audit_data$x_log_p99, audit_data$x_log_max))
cat(sprintf("\n    - NA proportion: %.4f", audit_data$x_log_na_prop))

# After scaling
cat("\n  After scaling (mean=0, sd=1):")
cat(sprintf("\n    - Range: %s", audit_data$x_scaled_range))
cat(sprintf("\n    - Mean: %.4f, SD: %.4f", audit_data$x_scaled_mean, audit_data$x_scaled_sd))
cat(sprintf("\n    - p50: %.4f, p90: %.4f, p99: %.4f, Max: %.4f", 
            audit_data$x_scaled_p50, audit_data$x_scaled_p90, 
            audit_data$x_scaled_p99, audit_data$x_scaled_max))
cat(sprintf("\n    - NA proportion: %.4f", audit_data$x_scaled_na_prop))

cat(sprintf("\n  - Glmnet alpha: %.1f", audit_data$alpha))
cat(sprintf("\n  - Glmnet standardize: %s", audit_data$standardize))
cat(sprintf("\n  - Lambda selection: %s", audit_data$lambda_selection))
cat(sprintf("\n  - Seed value: %d", audit_data$seed_value))
cat(sprintf("\n  - Audit file saved to: %s\n", audit_file))

# ============================================================================
# Run Stability Selection
# ============================================================================

# Check data scale before bootstrap
cat("\n[Step 1.5] Checking data scale...\n")

# Check original data scale (before transformation)
xq99_raw <- quantile(X_original, 0.99, na.rm=TRUE)
cat(sprintf("  Original X 99th percentile: %.2f\n", xq99_raw))

# Check if data appears to be counts-scale
if (xq99_raw > 1e5) {
  cat("  Message: Input appears counts-scale before log/scale (expected).\n")
}

# Check scaled data scale
x_quantile_99 <- quantile(X, 0.99, na.rm=TRUE)
cat(sprintf("  Scaled X 99th percentile: %.2f\n", x_quantile_99))

cat(sprintf("\n[Step 2] Running %d bootstrap iterations...\n", BOOTSTRAP_ITER))
cat("This may take 10-15 minutes...\n\n")

set.seed(2026)  # Reproducibility
selection_freq <- setNames(rep(0, ncol(X)), colnames(X))

pb <- txtProgressBar(min = 0, max = BOOTSTRAP_ITER, style = 3)

for (i in 1:BOOTSTRAP_ITER) {

  # Bootstrap sample (stratified)
  t2_high_idx <- which(y == 1)
  t2_low_idx <- which(y == 0)

  train_high <- sample(t2_high_idx, size = floor(TRAIN_FRAC * length(t2_high_idx)))
  train_low <- sample(t2_low_idx, size = floor(TRAIN_FRAC * length(t2_low_idx)))
  train_idx <- c(train_high, train_low)

  # Prepare training data - sample x gene format for glmnet
  X_train <- X[train_idx, ]
  y_train <- y[train_idx]

  # ElasticNet with cross-validation
  cv_fit <- cv.glmnet(X_train, y_train,
                      family = "binomial",
                      alpha = ALPHA,
                      nfolds = 10,
                      type.measure = "auc",
                      standardize = FALSE)  # Manual standardization already applied

  # Extract selected genes (non-zero coefficients at lambda.min)
  coef_mat <- coef(cv_fit, s = "lambda.min")
  selected_genes <- rownames(coef_mat)[which(coef_mat != 0)][-1]  # Remove intercept

  # Update selection frequency
  selection_freq[selected_genes] <- selection_freq[selected_genes] + 1

  setTxtProgressBar(pb, i)
}

close(pb)

# Calculate selection frequency percentage
selection_freq_pct <- selection_freq / BOOTSTRAP_ITER

cat(sprintf("\n\n✓ Stability selection complete!\n"))

# ============================================================================
# Apply selection threshold
# ============================================================================

cat(sprintf("\n[Step 3] Applying selection threshold (%.0f%%)...\n",
            SELECTION_THRESHOLD * 100))

final_genes <- names(selection_freq_pct[selection_freq_pct >= SELECTION_THRESHOLD])

cat(sprintf("\n✓ Final genes selected: %d\n", length(final_genes)))

# Show selection frequencies
final_results <- data.frame(
  Gene = names(selection_freq_pct),
  Selection_Frequency = selection_freq_pct,
  Selected = ifelse(selection_freq_pct >= SELECTION_THRESHOLD, "Yes", "No"),
  stringsAsFactors = FALSE
) %>%
  arrange(desc(Selection_Frequency))

head_results <- head(final_results, 20)
cat("\nTop 20 genes by selection frequency:\n")
print(head_results)

# ============================================================================
# Stability selection complete - weights will be locked in final model script
# ============================================================================

cat("\n[Step 4] Stability selection complete...\n")
cat("Note: Model weights will be locked in the final model script after gene selection is finalized.\n")

# ============================================================================
# Export results
# ============================================================================

cat("\n[Step 5] Exporting results...\n")

# Save full stability selection results
saveRDS(list(
  selection_freq = selection_freq_pct,
  final_genes = final_genes,
  parameters = list(
    bootstrap_iter = BOOTSTRAP_ITER,
    alpha = ALPHA,
    threshold = SELECTION_THRESHOLD
  )
), file.path(DATA_DERIVED_DIR, "stability_selection_results.rds"))

# Save final gene list
write_csv(
  final_results %>% filter(Selected == "Yes"),
  file.path(DATA_DERIVED_DIR, "stability_selected_genes.csv")
)

# Save bootstrap frequencies
write_csv(
  final_results,
  file.path(DATA_DERIVED_DIR, "bootstrap_frequencies.csv")
)

# Save stability selection results with type
results_filename <- ifelse(passed_genes_type == "fallback", 
                          "stability_selection_results_fallback.rds", 
                          "stability_selection_results_primary.rds")

saveRDS(list(
  selection_freq = selection_freq_pct,
  final_genes = final_genes,
  parameters = list(
    bootstrap_iter = BOOTSTRAP_ITER,
    alpha = ALPHA,
    threshold = SELECTION_THRESHOLD
  ),
  passed_genes_type = passed_genes_type,
  n_features_in_X = n_features_in_X,
  n_selected_ge_0p7 = length(final_genes)
), file.path(DATA_DERIVED_DIR, results_filename))

# Create stability comparison summary
stability_comparison <- data.frame(
  passed_genes_type = passed_genes_type,
  n_features_in_X = n_features_in_X,
  n_selected_ge_0p7 = length(final_genes),
  top_genes = paste(head(names(sort(selection_freq_pct, decreasing = TRUE)), 5), collapse=", "),
  stringsAsFactors = FALSE
)

# Try to read the other set if it exists for comparison
if (fallback_triggered) {
  other_results_file <- file.path(DATA_DERIVED_DIR, "stability_selection_results_primary.rds")
  other_type <- "primary"
} else {
  other_results_file <- file.path(DATA_DERIVED_DIR, "stability_selection_results_fallback.rds")
  other_type <- "fallback"
}

if (file.exists(other_results_file)) {
  other_results <- readRDS(other_results_file)
  other_genes <- other_results$final_genes
  intersection_size <- length(intersect(final_genes, other_genes))
  union_size <- length(union(final_genes, other_genes))
  jaccard_index <- intersection_size / union_size
  
  stability_comparison$other_type <- other_type
  stability_comparison$other_n_selected_ge_0p7 <- length(other_genes)
  stability_comparison$intersection_size <- intersection_size
  stability_comparison$jaccard_index <- jaccard_index
}

# Save comparison summary
write_csv(stability_comparison, file.path(DATA_DERIVED_DIR, "stability_comparison_summary.csv"))

cat("✓ Results exported\n")
cat(sprintf("  - Stability results: %s\n", results_filename))
cat(sprintf("  - Comparison summary: %s\n", file.path(DATA_DERIVED_DIR, "stability_comparison_summary.csv")))

# Run comparison mode if enabled
if (RUN_COMPARISON) {
  cat("\n[Comparison Mode] Running alternative threshold strategy...\n")
  
  # Determine which alternative to run
  if (passed_genes_type == "primary") {
    alt_type <- "fallback"
    alt_file <- "limma_passed_genes_fallback_fc1p5_fdr0p01.csv"
  } else {
    alt_type <- "primary"
    alt_file <- "limma_passed_genes_primary.csv"
  }
  
  cat(sprintf("  Running %s threshold strategy for comparison...\n", alt_type))
  
  # Load alternative passed genes
  alt_deg_pass <- read_csv(file.path(DATA_DERIVED_DIR, alt_file))
  alt_deg_genes <- alt_deg_pass$Gene
  
  # Ensure all alt_deg_genes are present in expr_diet
  alt_deg_genes_present <- intersect(alt_deg_genes, rownames(expr_diet))
  if (length(alt_deg_genes_present) != length(alt_deg_genes)) {
    stop("Some DEG genes not found in expression data")
  }
  
  alt_expr_deg <- expr_diet[alt_deg_genes_present, , drop=FALSE]
  alt_n_features_in_X <- length(alt_deg_genes_present)
  
  # Transpose to get samples × genes format for glmnet
  alt_X <- t(alt_expr_deg)
  
  # Apply universal scale normalization
  alt_X_original <- alt_X
  alt_X1 <- log2(alt_X + 1)
  alt_X2 <- scale(alt_X1)
  alt_X <- alt_X2
  
  # Run stability selection for alternative
  cat(sprintf("  Running %d bootstrap iterations for %s...\n", BOOTSTRAP_ITER, alt_type))
  
  set.seed(2026)  # Reproducibility
  alt_selection_freq <- setNames(rep(0, ncol(alt_X)), colnames(alt_X))
  
  for (i in 1:BOOTSTRAP_ITER) {
    # Bootstrap sample (stratified)
    t2_high_idx <- which(y == 1)
    t2_low_idx <- which(y == 0)

    train_high <- sample(t2_high_idx, size = floor(TRAIN_FRAC * length(t2_high_idx)))
    train_low <- sample(t2_low_idx, size = floor(TRAIN_FRAC * length(t2_low_idx)))
    train_idx <- c(train_high, train_low)

    # Prepare training data - sample x gene format for glmnet
    X_train <- alt_X[train_idx, ]
    y_train <- y[train_idx]

    # ElasticNet with cross-validation
    cv_fit <- cv.glmnet(X_train, y_train,
                        family = "binomial",
                        alpha = ALPHA,
                        nfolds = 10,
                        type.measure = "auc",
                        standardize = FALSE)  # Manual standardization already applied

    # Extract selected genes (non-zero coefficients at lambda.min)
    coef_mat <- coef(cv_fit, s = "lambda.min")
    selected_genes <- rownames(coef_mat)[which(coef_mat != 0)][-1]  # Remove intercept

    # Update selection frequency
    alt_selection_freq[selected_genes] <- alt_selection_freq[selected_genes] + 1
  }
  
  # Calculate selection frequency percentage
  alt_selection_freq_pct <- alt_selection_freq / BOOTSTRAP_ITER
  
  # Apply selection threshold
  alt_final_genes <- names(alt_selection_freq_pct[alt_selection_freq_pct >= SELECTION_THRESHOLD])
  
  # Save alternative results
  alt_results_filename <- paste0("stability_selection_results_", alt_type, ".rds")
  saveRDS(list(
    selection_freq = alt_selection_freq_pct,
    final_genes = alt_final_genes,
    parameters = list(
      bootstrap_iter = BOOTSTRAP_ITER,
      alpha = ALPHA,
      threshold = SELECTION_THRESHOLD
    ),
    passed_genes_type = alt_type,
    n_features_in_X = alt_n_features_in_X,
    n_selected_ge_0p7 = length(alt_final_genes)
  ), file.path(DATA_DERIVED_DIR, alt_results_filename))
  
  cat(sprintf("  ✓ Alternative %s results saved: %s\n", alt_type, alt_results_filename))
  
  # Create comprehensive comparison summary
  stability_comparison <- data.frame(
    passed_genes_type = c(passed_genes_type, alt_type),
    n_features_in_X = c(n_features_in_X, alt_n_features_in_X),
    n_selected_ge_0p7 = c(length(final_genes), length(alt_final_genes)),
    top_genes = c(
      paste(head(names(sort(selection_freq_pct, decreasing = TRUE)), 5), collapse=", "),
      paste(head(names(sort(alt_selection_freq_pct, decreasing = TRUE)), 5), collapse=", ")
    ),
    stringsAsFactors = FALSE
  )
  
  # Calculate overlap metrics
  intersection_size <- length(intersect(final_genes, alt_final_genes))
  union_size <- length(union(final_genes, alt_final_genes))
  jaccard_index <- intersection_size / union_size
  
  stability_comparison$intersection_size <- intersection_size
  stability_comparison$jaccard_index <- jaccard_index
  
  # Save comparison summary
  write_csv(stability_comparison, file.path(DATA_DERIVED_DIR, "stability_comparison_summary.csv"))
  
  cat(sprintf("  ✓ Comparison summary updated with overlap metrics\n"))
  cat(sprintf("  - Intersection size: %d\n", intersection_size))
  cat(sprintf("  - Jaccard index: %.4f\n", jaccard_index))
}

# Save final results to fixed interface file
final_results_filename <- "stability_selection_results_final.rds"
saveRDS(list(
  selection_freq = selection_freq_pct,
  final_genes = final_genes,
  parameters = list(
    bootstrap_iter = BOOTSTRAP_ITER,
    alpha = ALPHA,
    threshold = SELECTION_THRESHOLD
  ),
  passed_genes_type = passed_genes_type,
  n_features_in_X = n_features_in_X,
  n_selected_ge_0p7 = length(final_genes)
), file.path(DATA_DERIVED_DIR, final_results_filename))

cat(sprintf("  - Final results (fixed interface): %s\n", final_results_filename))

cat("\n=================================================================\n")
cat(sprintf("✓ Phase 5 Complete: %d genes selected and ready for finalization\n", length(final_genes)))
cat(sprintf("  Using %s passed genes\n", passed_genes_type))
if (RUN_COMPARISON) {
  cat("  Comparison mode: Both primary and fallback strategies evaluated\n")
}
cat("→ Next: Final model decision\n")
cat("=================================================================\n\n")
