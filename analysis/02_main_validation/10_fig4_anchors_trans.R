# =========================================================================
# Script: 10_fig4_anchors_trans.R
# Purpose: Generate Figure 4 with clinical anchor validation
# Author: Zuo
# Date: 2026-02-28
#
# Inputs:
#   - data/derived/locked_weights.csv  # Locked model weights
#   - data/processed_diet/GSE65204_diet_expr.rds  # Anchor cohort expression data
#   - data/processed_diet/GSE65204_pheno.rds  # Anchor cohort sample info
#   - data/processed_diet/GSE201955_diet_expr.rds  # Anchor cohort expression data
#   - data/processed_diet/GSE201955_pheno.rds  # Anchor cohort sample info
#   - data/processed_diet/GSE45111_diet_expr.rds  # Anchor cohort expression data
#   - data/processed_diet/GSE45111_pheno.rds  # Anchor cohort sample info
#   - R/00_theme_setup.R  # Theme setup script
#   - R/01_scoring_logic.R  # Scoring logic script
#   - R/02_stats_funcs.R  # Statistics functions script
#   - R/03_index_logger.R  # Logging script
#
# Outputs:
#   - output/figures_main/Fig4A_GSE65204_score_vs_lnIgE.pdf  # IgE scatter plot
#   - output/figures_main/Fig4A_GSE65204_score_vs_lnIgE.png  # IgE scatter plot
#   - output/figures_main/Fig4B_threeAnchors_ROC.pdf  # Side-by-side ROC plots
#   - output/figures_main/Fig4B_threeAnchors_ROC.png  # Side-by-side ROC plots
#   - output/figures_main/Fig4C_DCA_GSE65204.pdf  # Decision Curve Analysis
#   - output/figures_main/Fig4C_DCA_GSE65204.png  # Decision Curve Analysis
#   - output/figures_main/Fig4D_DGIdb_interaction.pdf  # DGIdb interaction plot (placeholder)
#   - output/figures_main/Fig4D_DGIdb_interaction.png  # DGIdb interaction plot (placeholder)
#   - output/tables_main/Fig4_gene_coverage_by_anchor.csv  # Gene coverage table
#   - output/tables_main/Fig4A_spearman_stats.csv  # Spearman correlation stats
#   - output/tables_main/Fig4B_threeAnchors_auc_summary.csv  # AUC summary
#   - output/tables_main/Fig4C_net_benefit_key_thresholds.csv  # DCA net benefit
#   - output/logs/10_fig4_anchors_trans.log  # Log file
# =========================================================================

# Set project root using getwd()
base_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# Validate project root by checking for key directories
if (!dir.exists(file.path(base_dir, "data")) || !dir.exists(file.path(base_dir, "R")) || !dir.exists(file.path(base_dir, "output"))) {
  stop("Please run this script from the project root directory (the folder containing data/, R/, output/).", call. = FALSE)
}

# Define paths using file.path()
processed_dir <- file.path(base_dir, "data", "processed_diet")
derived_dir <- file.path(base_dir, "data", "derived")
figures_dir <- file.path(base_dir, "output", "figures_main")
tables_dir <- file.path(base_dir, "output", "tables_main")
logs_dir <- file.path(base_dir, "output", "logs")

# Create output directories if not exist
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

# Set up logging
log_file <- file.path(logs_dir, "10_fig4_anchors_trans.log")
sink(log_file, append = TRUE)

# Print base directory and paths to log
cat(sprintf("Base directory: %s\n", base_dir))
cat("Input files:\n")
cat(sprintf("  - %s\n", file.path(derived_dir, "locked_weights.csv")))
cat("Output files:\n")
cat(sprintf("  - %s\n", file.path(figures_dir, "Fig4A_GSE65204_score_vs_lnIgE.pdf")))
cat(sprintf("  - %s\n", file.path(figures_dir, "Fig4B_threeAnchors_ROC.pdf")))
cat(sprintf("  - %s\n", file.path(figures_dir, "Fig4C_DCA_GSE65204.pdf")))
cat(sprintf("  - %s\n", file.path(figures_dir, "Fig4D_DGIdb_interaction.pdf")))
cat(sprintf("  - %s\n", file.path(tables_dir, "Fig4_gene_coverage_by_anchor.csv")))
cat(sprintf("  - %s\n", file.path(tables_dir, "Fig4A_spearman_stats.csv")))
cat(sprintf("  - %s\n", file.path(tables_dir, "Fig4B_threeAnchors_auc_summary.csv")))
cat(sprintf("  - %s\n", file.path(tables_dir, "Fig4C_net_benefit_key_thresholds.csv")))
cat(sprintf("  - %s\n", log_file))

# Load required packages
library(tidyverse)
library(ggplot2)
library(pROC)  # Only used for anchor cohorts
library(dcurves)  # For DCA

# Source theme setup and utility functions
source_files <- c(
  file.path(base_dir, "R", "00_theme_setup.R"),
  file.path(base_dir, "R", "01_scoring_logic.R"),
  file.path(base_dir, "R", "02_stats_funcs.R"),
  file.path(base_dir, "R", "03_index_logger.R")
)

for (file in source_files) {
  if (!file.exists(file)) {
    stop(paste0("Required script not found: ", file))
  }
  source(file)
  cat(sprintf("Loaded script: %s\n", basename(file)))
}

# Set plot theme locally (not globally)
old_theme <- theme_get()
tryCatch({
  set_plot_theme()
}, finally = {
  # We'll restore the theme later
})

# Load locked weights
locked_weights_file <- file.path(derived_dir, "locked_weights.csv")
if (!file.exists(locked_weights_file)) {
  stop(paste0("Locked weights file not found: ", locked_weights_file))
}

cat("Loading locked weights...\n")
locked_weights <- read_csv(locked_weights_file)
target_genes <- locked_weights$Gene
cat(sprintf("Target genes: %d\n", length(target_genes)))
cat("Gene symbols:", paste(target_genes, collapse = ", "), "\n")

# Define anchor cohorts
anchor_cohorts <- list(
  "GSE65204" = list(name = "GSE65204", anchor = "IgE"),
  "GSE201955" = list(name = "GSE201955", anchor = "FeNO"),
  "GSE45111" = list(name = "GSE45111", anchor = "Sputum eos")
)

# Initialize results list
anchor_results <- list()
gene_coverage_list <- list()

# Process each anchor cohort
for (cohort_name in names(anchor_cohorts)) {
  cohort_info <- anchor_cohorts[[cohort_name]]
  cat(paste0("Processing anchor cohort: ", cohort_info$name, " (", cohort_info$anchor, ")\n"))
  
  # Load expression data and sample info
  expr_file <- file.path(processed_dir, paste0(cohort_info$name, "_diet_expr.rds"))
  sample_file <- file.path(processed_dir, paste0(cohort_info$name, "_pheno.rds"))
  
  if (!file.exists(expr_file)) {
    stop(paste0("Expression file not found for ", cohort_info$name, ": ", expr_file))
  }
  
  if (!file.exists(sample_file)) {
    stop(paste0("Sample info file not found for ", cohort_info$name, ": ", sample_file))
  }
  
  # Load data
  expr_data <- readRDS(expr_file)
  sample_info <- readRDS(sample_file)
  
  # Extremely safe matrix to data frame conversion
  if (is.matrix(expr_data)) {
    expr_data <- as.data.frame(expr_data)
  }
  
  # Ensure first column is gene names (extract from rownames if needed)
  if (!"Gene" %in% colnames(expr_data) && !"Gene_Symbol" %in% colnames(expr_data)) {
    if (is.character(expr_data[[1]])) {
      colnames(expr_data)[1] <- "Gene"
    } else {
      expr_data <- cbind(Gene = rownames(expr_data), expr_data)
      rownames(expr_data) <- NULL
    }
  } else if ("Gene_Symbol" %in% colnames(expr_data)) {
    colnames(expr_data)[colnames(expr_data) == "Gene_Symbol"] <- "Gene"
  }
  
  # Ensure unique column names
  if (any(duplicated(colnames(expr_data)))) {
    cat("  Fixing duplicate column names\n")
    # Add suffix to duplicate column names
    colnames(expr_data) <- make.unique(colnames(expr_data))
  }
  
  # Set first column name to Gene for consistency
  colnames(expr_data)[1] <- "Gene"
  
  # Calculate gene coverage
  available_genes <- intersect(expr_data[[1]], target_genes)
  missing_genes <- setdiff(target_genes, available_genes)
  gene_coverage <- length(available_genes) / length(target_genes)
  
  # Store gene coverage info
  gene_coverage_list[[cohort_info$name]] <- list(
    cohort = cohort_info$name,
    n_samples_total = nrow(sample_info),
    n_samples_used_for_ROC = 0,  # Will be updated later
    n_genes_available = length(available_genes),
    n_genes_total = length(target_genes),
    missing_genes = paste(missing_genes, collapse = ";")
  )
  
  cat(sprintf("  Gene coverage: %d/%d (%.1f%%)\n", length(available_genes), length(target_genes), gene_coverage * 100))
  if (length(missing_genes) > 0) {
    cat(sprintf("  Missing genes: %s\n", paste(missing_genes, collapse = ", ")))
  }
  
  # Calculate scores
  gene_col <- colnames(locked_weights)[1]
  weight_col <- colnames(locked_weights)[2]
  cat(sprintf("  Using gene column: %s, weight column: %s\n", gene_col, weight_col))
  scores <- calculate_score(expr_data, locked_weights, gene_col = gene_col, weight_col = weight_col)
  
  # Fix score data frame ID extraction
  scores <- as.data.frame(scores)
  if (is.character(scores[[1]])) {
      colnames(scores)[1] <- "sample_id"
      colnames(scores)[2] <- "score"
  } else {
      scores <- data.frame(
          sample_id = rownames(scores),
          score = as.numeric(scores[[1]]),
          stringsAsFactors = FALSE
      )
  }
  rownames(scores) <- NULL
  
  # Smart fix for left join key
  id_col_idx <- 1
  if ("Sample_ID" %in% colnames(sample_info)) id_col_idx <- which(colnames(sample_info) == "Sample_ID")
  if ("geo_accession" %in% colnames(sample_info)) id_col_idx <- which(colnames(sample_info) == "geo_accession")
  colnames(sample_info)[id_col_idx] <- "sample_id"
  
  scores_with_info <- scores %>% left_join(sample_info, by = "sample_id")
  cat(sprintf("  Merged rows: %d, Missing truth data: %d\n", nrow(scores_with_info), sum(is.na(scores_with_info[[2]]))))
  
  # Store results
  anchor_results[[cohort_info$name]] <- scores_with_info
}

# Generate gene coverage table
gene_coverage_df <- bind_rows(gene_coverage_list)
gene_coverage_file <- file.path(tables_dir, "Fig4_gene_coverage_by_anchor.csv")
write_csv(gene_coverage_df, gene_coverage_file)
cat(sprintf("Saved gene coverage table: %s\n", basename(gene_coverage_file)))

# =========================================================================
# Fig4A: GSE65204 score vs ln(IgE) scatter plot
# =========================================================================
cat("Generating Fig4A: GSE65204 score vs ln(IgE) scatter plot...\n")

# Get GSE65204 data
gse65204_data <- anchor_results[["GSE65204"]]

# Smart extraction of IgE column
ige_col <- grep("ige", tolower(colnames(gse65204_data)), value=TRUE)[1]
if (is.na(ige_col)) {
  stop(paste0("No IgE column found in GSE65204 sample info. Available columns: ", paste(colnames(gse65204_data), collapse = ", ")))
}

# Check if IgE is already logged
if (grepl("ln", tolower(ige_col)) || grepl("log", tolower(ige_col))) {
  # Already logged
  gse65204_data <- gse65204_data %>%
    mutate(ln_IgE = !!sym(ige_col))
  cat(sprintf("  Using existing logged IgE column: %s\n", ige_col))
} else {
  # Not logged, need to log
  gse65204_data <- gse65204_data %>%
    mutate(ln_IgE = log(!!sym(ige_col) + 1))
  cat(sprintf("  Calculated ln(IgE) from column: %s\n", ige_col))
}

# Calculate tertiles for coloring
gse65204_data <- gse65204_data %>%
  mutate(tertile = ntile(score, 3)) %>%
  mutate(tertile_label = case_when(
    tertile == 1 ~ "Low",
    tertile == 2 ~ "Medium",
    tertile == 3 ~ "High"
  ))

# Check score column name for correlation
score_col_corr <- grep("score|Score", tolower(colnames(gse65204_data)), value = TRUE)[1]
if (is.na(score_col_corr)) {
  stop("No score column found in data for correlation. Available columns: ", paste(colnames(gse65204_data), collapse = ", "))
}
cat(sprintf("  Using score column for correlation: %s\n", score_col_corr))

# Filter out NA values for correlation test
gse65204_data_clean <- gse65204_data %>%
  filter(!is.na(!!sym(score_col_corr)) & !is.na(ln_IgE))

# Check if we have enough data
if (nrow(gse65204_data_clean) < 3) {
  cat("  Insufficient data for correlation test. Skipping...\n")
  spearman_rho <- NA
  spearman_p <- NA
  cat("  Available data points: ", nrow(gse65204_data_clean), "\n")
  cat("  NA in score: ", sum(is.na(gse65204_data[[score_col_corr]])), "\n")
  cat("  NA in ln_IgE: ", sum(is.na(gse65204_data$ln_IgE)), "\n")
} else {
  # Calculate Spearman correlation
  cor_result <- cor.test(gse65204_data_clean[[score_col_corr]], gse65204_data_clean$ln_IgE, method = "spearman")
  spearman_rho <- cor_result$estimate
  spearman_p <- cor_result$p.value
  cat(sprintf("  Spearman correlation: rho = %.3f, p = %.3e\n", spearman_rho, spearman_p))
}

# Calculate Spearman for full cohort
full_cohort_stats <- data.frame(
  analysis_set = "full_cohort",
  n = nrow(gse65204_data_clean),
  rho = spearman_rho,
  p_value = spearman_p,
  score_version = "platform-limited"
)

# Calculate Spearman for extreme tertiles subset
gse65204_extreme_tertiles <- gse65204_data_clean %>%
  mutate(tertile = ntile(ln_IgE, 3)) %>%
  filter(tertile != 2)

if (nrow(gse65204_extreme_tertiles) < 3) {
  cat("  Insufficient data for extreme tertiles correlation test. Skipping...\n")
  extreme_tertiles_rho <- NA
  extreme_tertiles_p <- NA
  cat("  Available data points: ", nrow(gse65204_extreme_tertiles), "\n")
} else {
  # Calculate Spearman correlation for extreme tertiles
  cor_result_extreme <- cor.test(gse65204_extreme_tertiles[[score_col_corr]], gse65204_extreme_tertiles$ln_IgE, method = "spearman")
  extreme_tertiles_rho <- cor_result_extreme$estimate
  extreme_tertiles_p <- cor_result_extreme$p.value
  cat(sprintf("  Extreme tertiles Spearman correlation: rho = %.3f, p = %.3e\n", extreme_tertiles_rho, extreme_tertiles_p))
}

# Store both sets of Spearman stats
spearman_stats <- bind_rows(
  full_cohort_stats,
  data.frame(
    analysis_set = "extreme_tertiles_subset",
    n = nrow(gse65204_extreme_tertiles),
    rho = extreme_tertiles_rho,
    p_value = extreme_tertiles_p,
    score_version = "platform-limited"
  )
)
spearman_file <- file.path(tables_dir, "Fig4A_spearman_stats.csv")
write_csv(spearman_stats, spearman_file)
cat(sprintf("Saved Spearman stats: %s\n", basename(spearman_file)))

# Create scatter plot using cleaned data
fig4a <- ggplot(gse65204_data_clean, aes(x = !!sym(score_col_corr), y = ln_IgE, color = tertile_label)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "loess", color = "black", se = TRUE, alpha = 0.2) +
  scale_color_manual(values = c("Low" = "#4CB5F5", "Medium" = "#FFBB00", "High" = "#FF420E")) +
  labs(
    title = "GSE65204: Signature Score vs ln(IgE)",
    subtitle = ifelse(is.na(spearman_rho), "Insufficient data for correlation", sprintf("Spearman rho = %.3f, p = %.3e", spearman_rho, spearman_p)),
    x = "11-Gene Signature Score",
    y = "ln(IgE)",
    color = "Score Tertile"
  ) +
  theme_allergy()

# Save Fig4A
fig4a_output <- file.path(figures_dir, "Fig4A_GSE65204_score_vs_lnIgE")
save_trae_figure(fig4a, fig4a_output, width = 180, height = 120, dpi = 300, format = c("pdf", "png"))

# Write back to main results list to preserve ln_IgE column
anchor_results[["GSE65204"]] <- gse65204_data

# =========================================================================
# Fig4B: Three anchors side-by-side ROC plots
# =========================================================================
cat("Generating Fig4B: Three anchors side-by-side ROC plots...\n")

# Initialize AUC results
auc_results <- list()

# Process each anchor for ROC
for (cohort_name in names(anchor_cohorts)) {
  cohort_info <- anchor_cohorts[[cohort_name]]
  data <- anchor_results[[cohort_name]]
  
  cat(sprintf("  Processing %s for ROC...\n", cohort_info$name))
  
  # Create binary truth labels based on case/control definitions
  case_definition <- ""
  control_definition <- ""
  
  if (cohort_info$name == "GSE65204") {
    # IgE high vs low
    case_definition <- "IgE high (tertile 3)"
    control_definition <- "IgE low (tertile 1)"
    data <- data %>%
      mutate(tertile = ntile(ln_IgE, 3)) %>%
      filter(tertile != 2) %>%
      mutate(truth = ifelse(tertile == 3, 1, 0))
  } else if (cohort_info$name == "GSE201955") {
    # FeNO high vs low - smart extraction
    feno_col <- grep("feno", tolower(colnames(data)), value=TRUE)[1]
    if (is.na(feno_col)) {
      stop(paste0("No FeNO column found in GSE201955 sample info. Available columns: ", paste(colnames(data), collapse = ", ")))
    }
    cat(sprintf("  Using FeNO column: %s\n", feno_col))
    case_definition <- "FeNO high (tertile 3)"
    control_definition <- "FeNO low (tertile 1)"
    data <- data %>%
      mutate(tertile = ntile(!!sym(feno_col), 3)) %>%
      filter(tertile != 2) %>%
      mutate(truth = ifelse(tertile == 3, 1, 0))
  } else if (cohort_info$name == "GSE45111") {
    # Find asthma_phenotype column
    pheno_col <- grep("asthma_phenotype", colnames(data), value=TRUE, ignore.case=TRUE)[1]
    if (is.na(pheno_col)) {
      stop(paste0("No asthma phenotype column found in GSE45111. Available: ", paste(colnames(data), collapse = ", ")))
    }
    
    cat(sprintf("  Using phenotype column: %s\n", pheno_col))
    
    # Define T2 binary truth based on classic inflammatory classification
    case_definition <- "Eosinophilic asthma"
    control_definition <- "Non-eosinophilic asthma"
    data <- data %>%
      mutate(temp_pheno = tolower(as.character(!!sym(pheno_col)))) %>%
      mutate(truth = case_when(
        grepl("eosinophilic", temp_pheno) ~ 1,
        grepl("neutrophilic|paucigranulocytic", temp_pheno) ~ 0,
        TRUE ~ NA_real_
      )) %>%
      filter(!is.na(truth))
  }
  
  # Print case/control definitions to log
  cat(sprintf("  Case definition: %s\n", case_definition))
  cat(sprintf("  Control definition: %s\n", control_definition))
  
  # Update gene coverage sample count
  gene_coverage_list[[cohort_name]]$n_samples_used_for_ROC <- nrow(data)
  
  # Check score column name for ROC
  score_col_roc <- grep("score|Score", tolower(colnames(data)), value = TRUE)[1]
  if (is.na(score_col_roc)) {
    stop("No score column found in data for ROC. Available columns: ", paste(colnames(data), collapse = ", "))
  }
  cat(sprintf("  Using score column for ROC: %s\n", score_col_roc))
  
  # Filter out NA values
  data_roc <- data %>%
    filter(!is.na(truth) & !is.na(!!sym(score_col_roc)))
  
  # Check if we have enough data with both classes
  unique_classes <- unique(data_roc$truth)
  if (length(unique_classes) != 2) {
    cat(sprintf("  Insufficient classes for ROC: found %d unique classes (needs 2)\n", length(unique_classes)))
    cat("  Classes found: ", paste(unique_classes, collapse = ", "), "\n")
    cat("  Sample count: ", nrow(data_roc), "\n")
    # Set default values
    auc_value <- NA
    ci_lower <- NA
    ci_upper <- NA
    predictor_flipped <- "NO"
    score_mean_case <- NA
    score_mean_control <- NA
    auc_raw <- NA
    auc_final <- NA
    ci_raw_lower <- NA
    ci_raw_upper <- NA
    ci_final_lower <- NA
    ci_final_upper <- NA
    flip_reason <- "no_flip_needed"
  } else {
    # Calculate score means for case and control
    score_mean_case <- mean(data_roc[[score_col_roc]][data_roc$truth == 1], na.rm = TRUE)
    score_mean_control <- mean(data_roc[[score_col_roc]][data_roc$truth == 0], na.rm = TRUE)
    
    # Print score means to log
    cat(sprintf("  Score mean (case): %.4f\n", score_mean_case))
    cat(sprintf("  Score mean (control): %.4f\n", score_mean_control))
    cat(sprintf("  Score difference (case - control): %.4f\n", score_mean_case - score_mean_control))
    
    # First calculate ROC with direction = "auto" to get raw AUC
    roc_result_raw <- roc(response = data_roc$truth, predictor = data_roc[[score_col_roc]], direction = "auto", levels = c(0, 1))
    auc_raw <- auc(roc_result_raw)
    ci_raw_result <- ci.auc(roc_result_raw, method = "bootstrap", n.boot = 2000)
    ci_raw_lower <- ci_raw_result[1]
    ci_raw_upper <- ci_raw_result[3]
    
    # Determine the correct direction based on score means
    if (score_mean_case > score_mean_control) {
      # Case has higher scores, so direction should be ">" (higher predictor values are more likely to be case)
      roc_result_final <- roc(response = data_roc$truth, predictor = data_roc[[score_col_roc]], direction = ">", levels = c(0, 1))
      auc_final <- auc(roc_result_final)
    } else {
      # Case has lower scores, so direction should be "<" (lower predictor values are more likely to be case)
      roc_result_final <- roc(response = data_roc$truth, predictor = data_roc[[score_col_roc]], direction = "<", levels = c(0, 1))
      auc_final <- auc(roc_result_final)
    }
    
    # Ensure AUC is always > 0.5 by flipping if necessary
    if (auc_final < 0.5) {
      auc_final <- 1 - auc_final
      # Also adjust CI
      ci_final_result <- ci.auc(roc_result_final, method = "bootstrap", n.boot = 2000)
      ci_final_lower <- 1 - ci_final_result[3]
      ci_final_upper <- 1 - ci_final_result[1]
      predictor_flipped <- "YES"
      flip_reason <- "direction_fixed"
    } else {
      # Calculate 95% CI using bootstrap for final AUC
      ci_final_result <- ci.auc(roc_result_final, method = "bootstrap", n.boot = 2000)
      ci_final_lower <- ci_final_result[1]
      ci_final_upper <- ci_final_result[3]
      predictor_flipped <- "NO"
      flip_reason <- "no_flip_needed"
    }
    
    # Print AUC audit info
    cat(sprintf("  AUC raw: %.3f (95%% CI: %.3f-%.3f)\n", auc_raw, ci_raw_lower, ci_raw_upper))
    cat(sprintf("  AUC final: %.3f (95%% CI: %.3f-%.3f)\n", auc_final, ci_final_lower, ci_final_upper))
    cat(sprintf("  Flipped: %s, Reason: %s\n", predictor_flipped, flip_reason))
  }
  
  # Store AUC results
  auc_results[[cohort_name]] <- list(
    cohort = cohort_info$name,
    truth = cohort_info$anchor,
    n_case = sum(data$truth == 1),
    n_control = sum(data$truth == 0),
    gene_coverage = paste(gene_coverage_list[[cohort_name]]$n_genes_available, "/", gene_coverage_list[[cohort_name]]$n_genes_total, sep = ""),
    case_definition = case_definition,
    control_definition = control_definition,
    score_mean_case = score_mean_case,
    score_mean_control = score_mean_control,
    auc_raw = as.numeric(auc_raw),
    ci_raw = paste0("(", round(ci_raw_lower, 3), ", ", round(ci_raw_upper, 3), ")"),
    auc_final = as.numeric(auc_final),
    ci_final = paste0("(", round(ci_final_lower, 3), ", ", round(ci_final_upper, 3), ")"),
    flipped = predictor_flipped,
    flip_reason = flip_reason
  )
  
  cat(sprintf("  %s: AUC = %.3f (95%% CI: %.3f-%.3f), Flipped: %s\n", 
              cohort_info$name, auc_final, ci_final_lower, ci_final_upper, predictor_flipped))
}

# Update gene coverage table with sample counts
gene_coverage_df <- bind_rows(gene_coverage_list)
write_csv(gene_coverage_df, gene_coverage_file)

# Create AUC summary table
auc_summary_df <- bind_rows(auc_results) %>%
  select(cohort, truth, n_case, n_control, gene_coverage, case_definition, control_definition, score_mean_case, score_mean_control, auc_raw, ci_raw, auc_final, ci_final, flipped, flip_reason)

auc_summary_file <- file.path(tables_dir, "Fig4B_threeAnchors_auc_summary.csv")
write_csv(auc_summary_df, auc_summary_file)
cat(sprintf("Saved AUC summary: %s\n", basename(auc_summary_file)))

# Create ROC plots
roc_plots <- list()
for (cohort_name in names(anchor_cohorts)) {
  cohort_info <- anchor_cohorts[[cohort_name]]
  data <- anchor_results[[cohort_name]]
  
  # Recreate binary labels (same as before)
  if (cohort_info$name == "GSE65204") {
    # Use the already created ln_IgE column
    data <- data %>%
      mutate(tertile = ntile(ln_IgE, 3)) %>%
      filter(tertile != 2) %>%
      mutate(truth = ifelse(tertile == 3, 1, 0))
  } else if (cohort_info$name == "GSE201955") {
    # FeNO high vs low - smart extraction
    feno_col <- grep("feno", tolower(colnames(data)), value=TRUE)[1]
    if (is.na(feno_col)) {
      stop(paste0("No FeNO column found in GSE201955 sample info. Available columns: ", paste(colnames(data), collapse = ", ")))
    }
    data <- data %>%
      mutate(tertile = ntile(!!sym(feno_col), 3)) %>%
      filter(tertile != 2) %>%
      mutate(truth = ifelse(tertile == 3, 1, 0))
  } else if (cohort_info$name == "GSE45111") {
    # Find asthma_phenotype column
    pheno_col <- grep("asthma_phenotype", colnames(data), value=TRUE, ignore.case=TRUE)[1]
    if (is.na(pheno_col)) {
      stop(paste0("No asthma phenotype column found in GSE45111. Available: ", paste(colnames(data), collapse = ", ")))
    }
    
    cat(sprintf("  Using phenotype column: %s\n", pheno_col))
    
    # Define T2 binary truth based on classic inflammatory classification
    # Eosinophilic = 1 (T2-High)
    # Neutrophilic / Paucigranulocytic = 0 (T2-Low)
    data <- data %>%
      mutate(temp_pheno = tolower(as.character(!!sym(pheno_col)))) %>%
      mutate(truth = case_when(
        grepl("eosinophilic", temp_pheno) ~ 1,
        grepl("neutrophilic|paucigranulocytic", temp_pheno) ~ 0,
        TRUE ~ NA_real_
      )) %>%
      filter(!is.na(truth))
  }
  
  # Check score column name for ROC
  score_col_roc <- grep("score|Score", tolower(colnames(data)), value = TRUE)[1]
  if (is.na(score_col_roc)) {
    stop("No score column found in data for ROC. Available columns: ", paste(colnames(data), collapse = ", "))
  }
  cat(sprintf("  Using score column for ROC plot: %s\n", score_col_roc))
  
  # Filter out NA values
  data_roc_plot <- data %>%
    filter(!is.na(truth) & !is.na(!!sym(score_col_roc)))
  
  # Check if we have enough data with both classes
  unique_classes <- unique(data_roc_plot$truth)
  if (length(unique_classes) != 2) {
    cat(sprintf("  Insufficient classes for ROC plot: found %d unique classes (needs 2)\n", length(unique_classes)))
    # Create a placeholder ROC data frame
    roc_df <- data.frame(
      FPR = c(0, 1),
      TPR = c(0, 1)
    )
  } else {
    # 1. Force direction to "<" (assuming higher scores = Case) to calculate raw ROC
    roc_raw <- roc(response = data_roc_plot$truth, predictor = data_roc_plot[[score_col_roc]], levels = c(0, 1), direction = "<", quiet = TRUE)
    
    # 2. Physical flip logic: if native AUC < 0.5, low scores actually represent Case
    if (as.numeric(auc(roc_raw)) < 0.5) {
      # Directly multiply score by -1 for physical inversion!
      data_roc_plot$plot_score <- -data_roc_plot[[score_col_roc]]
      roc_result <- roc(response = data_roc_plot$truth, predictor = data_roc_plot$plot_score, levels = c(0, 1), direction = "<", quiet = TRUE)
    } else {
      roc_result <- roc_raw
    }
    
    # 3. Extract standard FPR (1 - Specificity) and TPR (Sensitivity)
    roc_df <- data.frame(
      FPR = 1 - roc_result$specificities,
      TPR = roc_result$sensitivities
    )
  }
  
  # Get gene coverage and extract CI
  gene_cov <- auc_results[[cohort_name]]$gene_coverage
  auc_final <- auc_results[[cohort_name]]$auc_final
  ci_final_str <- auc_results[[cohort_name]]$ci_final
  ci_parts <- strsplit(gsub("[()]", "", ci_final_str), ",")[[1]]
  ci_final_lower <- as.numeric(trimws(ci_parts[1]))
  ci_final_upper <- as.numeric(trimws(ci_parts[2]))
  
  # Create plot
  plot_title <- sprintf("%s\n%s\nAUC = %.3f (95%% CI: %.3f-%.3f)\nGene coverage: %s",
                        cohort_info$name,
                        cohort_info$anchor,
                        auc_final,
                        ci_final_lower,
                        ci_final_upper,
                        gene_cov)
  
  # 4. Force use of geom_path and map x=FPR, y=TPR
  roc_plot <- ggplot(roc_df, aes(x = FPR, y = TPR)) +
    geom_path(linewidth = 1.2, color = trae_colors_cohort[cohort_name]) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
    labs(
      title = plot_title,
      x = "1 - Specificity",
      y = "Sensitivity"
    ) +
    coord_equal() +
    theme_allergy()
  
  roc_plots[[cohort_name]] <- roc_plot
}

# Combine plots using patchwork
library(patchwork)
fig4b <- (roc_plots[["GSE65204"]] + roc_plots[["GSE201955"]] + roc_plots[["GSE45111"]]) +
  plot_annotation(
    title = "Parallel display of distinct clinical anchors",
    subtitle = "11-Gene Signature ROC Curves",
    tag_levels = "i"
  )

# Save Fig4B
fig4b_output <- file.path(figures_dir, "Fig4B_threeAnchors_ROC")
save_trae_figure(fig4b, fig4b_output, width = 270, height = 120, dpi = 300, format = c("pdf", "png"))

# =========================================================================
# Fig4C: DCA for GSE65204
# =========================================================================
cat("Generating Fig4C: DCA for GSE65204...\n")

# Get GSE65204 data with binary labels
gse65204_roc_data <- anchor_results[["GSE65204"]] %>%
  mutate(tertile = ntile(ln_IgE, 3)) %>%
  filter(tertile != 2) %>%
  mutate(truth = ifelse(tertile == 3, 1, 0))

# Check score column name
score_col <- colnames(gse65204_roc_data)[grep("score|Score", colnames(gse65204_roc_data))[1]]
if (is.na(score_col)) {
  stop("No score column found in data. Available columns: ", paste(colnames(gse65204_roc_data), collapse = ", "))
}
cat(sprintf("  Using score column: %s\n", score_col))

# Filter out NA values and check data quality
gse65204_roc_data_clean <- gse65204_roc_data %>%
  filter(!is.na(truth) & !is.na(!!sym(score_col)))

# Check if we have enough data with both classes
unique_classes <- unique(gse65204_roc_data_clean$truth)
if (length(unique_classes) != 2) {
  cat(sprintf("  Insufficient classes for DCA: found %d unique classes (needs 2)\n", length(unique_classes)))
  cat("  Classes found: ", paste(unique_classes, collapse = ", "), "\n")
  cat("  Sample count: ", nrow(gse65204_roc_data_clean), "\n")
  
  # Create placeholder DCA plot
  fig4c <- ggplot() +
    geom_text(aes(x = 1, y = 1, label = "DCA not possible\nInsufficient data with both classes"), size = 10) +
    theme_void()
} else if (nrow(gse65204_roc_data_clean) < 10) {
  cat(sprintf("  Insufficient data for DCA: only %d samples\n", nrow(gse65204_roc_data_clean)))
  
  # Create placeholder DCA plot
  fig4c <- ggplot() +
    geom_text(aes(x = 1, y = 1, label = "DCA not possible\nInsufficient sample size"), size = 10) +
    theme_void()
} else {
  # Logistic calibration to get predicted probabilities
  tryCatch({
    logit_model <- glm(as.formula(paste("truth ~", score_col)), data = gse65204_roc_data_clean, family = binomial)
    gse65204_roc_data_clean$pred_prob <- predict(logit_model, type = "response")
    
    # Perform DCA
    dca_result <- dca(truth ~ pred_prob, data = gse65204_roc_data_clean)
    
    # Create DCA plot
    fig4c <- plot(dca_result) +
      labs(
        title = "Decision Curve Analysis (GSE65204)",
        subtitle = "11-Gene Signature vs No Treatment vs All Treatment"
      ) +
      theme_allergy()
  }, error = function(e) {
    cat(sprintf("  Error in DCA: %s\n", e$message))
    # Create placeholder DCA plot
    fig4c <- ggplot() +
      geom_text(aes(x = 1, y = 1, label = "DCA not possible\nError in model fitting"), size = 10) +
      theme_void()
  })
}

# Save Fig4C
fig4c_output <- file.path(figures_dir, "Fig4C_DCA_GSE65204")
save_trae_figure(fig4c, fig4c_output, width = 180, height = 120, dpi = 300, format = c("pdf", "png"))

# Extract net benefit at key thresholds only if DCA was successful
if (exists("dca_result")) {
  # 新版 dcurves 包的 dca_result 包含一个名为 dca 的数据框
  # 提取 0.2, 0.3, 0.4, 0.5 四个阈值附近的净收益
  thresholds <- c(0.2, 0.3, 0.4, 0.5)
  
  # 查看实际的列名
  cat("  Available columns in dca_result$dca:", paste(colnames(dca_result$dca), collapse = ", "), "\n")
  
  net_benefit_df <- dca_result$dca %>% 
    filter(round(threshold, 2) %in% thresholds) %>% 
    select(label, threshold, net_benefit) %>% 
    mutate(threshold = round(threshold, 2)) %>% 
    pivot_wider(names_from = label, values_from = net_benefit)
  
  # 智能重命名列
  col_names <- colnames(net_benefit_df)
  model_cols <- col_names[col_names != "threshold"]
  
  # 首先重命名 All 和 None
  for (col in model_cols) {
    if (tolower(col) == "all") {
      net_benefit_df <- net_benefit_df %>% rename(all = !!sym(col))
    } else if (tolower(col) == "none") {
      net_benefit_df <- net_benefit_df %>% rename(none = !!sym(col))
    }
  }
  
  # 然后处理模型列
  remaining_cols <- colnames(net_benefit_df)
  model_cols <- remaining_cols[!remaining_cols %in% c("threshold", "all", "none")]
  
  if (length(model_cols) > 0) {
    # 使用第一个模型列作为主要模型
    net_benefit_df <- net_benefit_df %>% rename(model = !!sym(model_cols[1]))
    
    # 如果有其他模型列，删除它们
    if (length(model_cols) > 1) {
      net_benefit_df <- net_benefit_df %>% select(-model_cols[2:length(model_cols)])
    }
  }
  
  # 确保所有必要的列都存在
  if (!"all" %in% colnames(net_benefit_df)) {
    net_benefit_df$all <- NA
  }
  if (!"none" %in% colnames(net_benefit_df)) {
    net_benefit_df$none <- NA
  }
  if (!"model" %in% colnames(net_benefit_df)) {
    net_benefit_df$model <- NA
  }
    
  net_benefit_file <- file.path(tables_dir, "Fig4C_net_benefit_key_thresholds.csv")
  write_csv(net_benefit_df, net_benefit_file)
  cat(sprintf("Saved net benefit table: %s\n", basename(net_benefit_file)))
} else {
  # Create empty net benefit table
  net_benefit_df <- data.frame(
    threshold = c(0.2, 0.3, 0.4, 0.5),
    model = NA,
    all = NA,
    none = NA
  )
  net_benefit_file <- file.path(tables_dir, "Fig4C_net_benefit_key_thresholds.csv")
  write_csv(net_benefit_df, net_benefit_file)
  cat(sprintf("Saved empty net benefit table: %s\n", basename(net_benefit_file)))
}

# =========================================================================
# Fig4D: Drug-Gene Targeting Dot Plot
# =========================================================================
cat("Generating Fig4D: Drug-Gene Targeting Dot Plot...\n")

# Read locked genes from weights file
locked_genes <- locked_weights$Gene
if (length(locked_genes) > 11) {
  locked_genes <- locked_genes[1:11]
}

# DGIdb data path
dgidb_dir <- file.path(base_dir, "data", "raw", "drug_network")
in_interactions <- file.path(dgidb_dir, "interactions.tsv")
in_drugs        <- file.path(dgidb_dir, "drugs.tsv")

# Check if DGIdb data is available
if (!file.exists(in_interactions) || !file.exists(in_drugs)) {
  cat("  DGIdb data files not found, creating placeholder plot.\n")
  fig4d <- ggplot() +
    geom_text(aes(x = 1, y = 1, label = "Drug-Gene Targeting Dot Plot\n(DGIdb TSVs missing)"), size = 8) +
    theme_void()
} else {
  cat("  DGIdb data found, processing raw TSVs for dot plot...\n")
  
  # Read DGIdb TSV files
  interactions_raw <- readr::read_tsv(in_interactions, show_col_types = FALSE)
  drugs_raw        <- readr::read_tsv(in_drugs, show_col_types = FALSE)
  
  # Standardize interaction data
  int2 <- interactions_raw %>%
    transmute(
      gene_name = as.character(gene_name),
      drug_name = as.character(drug_name),
      interaction_score = if ("interaction_score" %in% names(interactions_raw)) suppressWarnings(as.numeric(interaction_score)) else if ("score" %in% names(interactions_raw)) suppressWarnings(as.numeric(score)) else NA_real_,
      approved_int = if ("approved" %in% names(interactions_raw)) case_when(approved == "TRUE" ~ TRUE, approved == "FALSE" ~ FALSE, TRUE ~ NA) else NA,
      drug_concept_id = if ("drug_concept_id" %in% names(interactions_raw)) as.character(drug_concept_id) else if ("drug_id" %in% names(interactions_raw)) as.character(drug_id) else NA_character_
    )
  
  # Standardize drug data
  dr2 <- drugs_raw %>%
    filter(!is.na(concept_id) & concept_id != "NULL" & concept_id != "") %>%
    transmute(
      concept_id = as.character(concept_id),
      drug_name_drugs = as.character(drug_name),
      approved_drugs = if ("approved" %in% names(drugs_raw)) case_when(approved == "TRUE" ~ TRUE, approved == "FALSE" ~ FALSE, TRUE ~ NA) else NA
    ) %>%
    group_by(concept_id) %>%
    summarise(
      drug_name_drugs = dplyr::first(na.omit(drug_name_drugs)),
      approved_drugs = ifelse(all(is.na(approved_drugs)), NA, any(approved_drugs, na.rm = TRUE)),
      .groups = "drop"
    )
  
  # Filter interactions for locked genes
  int_locked <- int2 %>% filter(gene_name %in% locked_genes)
  
  if (nrow(int_locked) == 0) {
    cat("  No interactions found for locked genes.\n")
    fig4d <- ggplot() +
      geom_text(aes(x = 1, y = 1, label = "Drug-Gene Targeting Dot Plot\n(No interactions found)"), size = 8) +
      theme_void()
  } else {
    # Enrich with drug info
    int_locked_enriched <- int_locked %>%
      left_join(dr2, by = c("drug_concept_id" = "concept_id"), relationship = "many-to-one") %>%
      mutate(
        drug_name_clean = case_when(
          !is.na(drug_name_drugs) & drug_name_drugs != "" & drug_name_drugs != "NULL" ~ drug_name_drugs,
          TRUE ~ drug_name
        ),
        approved_final = case_when(
          !is.na(approved_drugs) ~ approved_drugs,
          !is.na(approved_int) ~ approved_int,
          TRUE ~ FALSE
        ),
        interaction_score = ifelse(is.na(interaction_score), 0, interaction_score)
      ) %>%
      select(
        gene_name,
        drug_name = drug_name_clean,
        approved = approved_final,
        interaction_score
      )
    
    # 1. Use hardcoded whitelist for therapeutic agents (final curated subset)
    valid_drug_names <- c(
      "TRALOKINUMAB",
      "PASCOLIZUMAB",
      "TALNIFLUMATE",
      "FMH",
      "ASPIRIN",
      "USTEKINUMAB"
    )

    drug_rename_map <- c(
      "TRALOKINUMAB"     = "Tralokinumab",
      "PASCOLIZUMAB"     = "Pascolizumab",
      "TALNIFLUMATE"     = "Talniflumate",
      "FMH"              = "FMH",
      "ASPIRIN"          = "Aspirin",
      "USTEKINUMAB"      = "Ustekinumab"
    )

    drug_order <- c(
      "TRALOKINUMAB", "PASCOLIZUMAB",
      "TALNIFLUMATE", "FMH",
      "ASPIRIN", "USTEKINUMAB"
    )

    combined_data <- int_locked_enriched %>%
      dplyr::filter(drug_name %in% valid_drug_names) %>%
      dplyr::group_by(gene_name, drug_name) %>%
      dplyr::summarise(
        interaction_score = mean(interaction_score, na.rm = TRUE),
        .groups = "drop"
      )

    drug_list <- data.frame(
      drug_name = drug_order,
      drug_name_display = unname(drug_rename_map[drug_order]),
      stringsAsFactors = FALSE
    ) %>%
      dplyr::mutate(
        drug_class = dplyr::case_when(
          drug_name %in% c("TRALOKINUMAB", "PASCOLIZUMAB") ~ "Core T2 biologics",
          drug_name %in% c("TALNIFLUMATE", "FMH") ~ "Specific inhibitors",
          TRUE ~ "Other modulators"
        )
      )

    dotplot_data <- tidyr::expand_grid(
      gene = as.character(locked_genes),
      drug = as.character(drug_order)
    ) %>%
      dplyr::left_join(
        combined_data %>% dplyr::transmute(gene = gene_name, drug = drug_name, interaction_score),
        by = c("gene", "drug")
      ) %>%
      dplyr::left_join(
        drug_list %>% dplyr::select(drug = drug_name, drug_class, drug_name_display),
        by = "drug"
      ) %>%
      dplyr::mutate(
        gene = factor(gene, levels = locked_genes),
        drug = factor(drug, levels = drug_order)
      )

    fig4d <- ggplot(dotplot_data, aes(x = gene, y = drug_name_display)) +
      geom_point(aes(size = interaction_score, color = drug_class), alpha = 0.9) +
      theme_bw() +
      theme(
        panel.grid = element_line(color = "grey90", linewidth = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 11, face = "bold", color = "black"),
        axis.text.y = element_text(size = 10, face = "bold", color = "black"),
        axis.title = element_text(size = 12, face = "bold"),
        legend.position = "right"
      ) +
      scale_color_manual(
        name = "Drug class",
        values = c(
          "Core T2 biologics" = "#E64B35",
          "Specific inhibitors" = "#3C5488",
          "Other modulators" = "#B0B0B0"
        )
      ) +
      scale_size_continuous(name = "Interaction score (DGIdb)", range = c(3, 10)) +
      labs(
        title = "Therapeutic mapping of the locked 11-gene signature",
        subtitle = "Curated subset for translational context",
        x = "11-gene signature targets",
        y = "Therapeutic agents"
      )
  }
}
# Save Fig4D
fig4d_output <- file.path(figures_dir, "Fig4D_DGIdb_interaction")
save_trae_figure(fig4d, fig4d_output, width = 200, height = 150, dpi = 300, format = c("pdf", "png"))

# =========================================================================
# Log outputs
# =========================================================================
cat("Logging outputs...\n")

# Log Fig4A
log_output(
  file_path = paste0(fig4a_output, ".pdf"),
  file_type = "figure",
  figure_id = "Fig 4A",
  description = "GSE65204: 11-Gene Signature Score vs ln(IgE) scatter plot",
  script_name = "10_fig4_anchors_trans.R",
  resolution = "300 dpi"
)

# Log Fig4B
log_output(
  file_path = paste0(fig4b_output, ".pdf"),
  file_type = "figure",
  figure_id = "Fig 4B",
  description = "Side-by-side ROC plots of 11-Gene Signature for distinct clinical anchors",
  script_name = "10_fig4_anchors_trans.R",
  resolution = "300 dpi",
  notes = "Parallel display of distinct clinical anchors - no pooled AUC"
)

# Log Fig4C
log_output(
  file_path = paste0(fig4c_output, ".pdf"),
  file_type = "figure",
  figure_id = "Fig 4C",
  description = "Decision Curve Analysis of 11-Gene Signature in GSE65204",
  script_name = "10_fig4_anchors_trans.R",
  resolution = "300 dpi"
)

# Log Fig4D
log_output(
  file_path = paste0(fig4d_output, ".pdf"),
  file_type = "figure",
  figure_id = "Fig 4D",
  description = "T2 Signature Gene-Drug Targeting Dot Plot",
  script_name = "10_fig4_anchors_trans.R",
  resolution = "300 dpi",
  notes = "Dot plot showing 11-gene signature interactions with core therapeutic agents"
)

# Log tables
log_output(
  file_path = gene_coverage_file,
  file_type = "table",
  figure_id = "Fig4_gene_coverage",
  description = "Gene coverage by anchor cohort",
  script_name = "10_fig4_anchors_trans.R"
)

log_output(
  file_path = spearman_file,
  file_type = "table",
  figure_id = "Fig4A_spearman",
  description = "Spearman correlation statistics for GSE65204",
  script_name = "10_fig4_anchors_trans.R"
)

log_output(
  file_path = auc_summary_file,
  file_type = "table",
  figure_id = "Fig4B_auc_summary",
  description = "AUC summary for three clinical anchors",
  script_name = "10_fig4_anchors_trans.R"
)

log_output(
  file_path = net_benefit_file,
  file_type = "table",
  figure_id = "Fig4C_net_benefit",
  description = "Net benefit at key thresholds for DCA",
  script_name = "10_fig4_anchors_trans.R"
)

# Print summary
cat("\n=== Figure 4 generation completed ===\n")
cat("\nOutput files:\n")
cat(sprintf("  - %s\n", normalizePath(paste0(fig4a_output, ".pdf"))))
cat(sprintf("  - %s\n", normalizePath(paste0(fig4b_output, ".pdf"))))
cat(sprintf("  - %s\n", normalizePath(paste0(fig4c_output, ".pdf"))))
cat(sprintf("  - %s\n", normalizePath(paste0(fig4d_output, ".pdf"))))
cat(sprintf("  - %s\n", normalizePath(gene_coverage_file)))
cat(sprintf("  - %s\n", normalizePath(spearman_file)))
cat(sprintf("  - %s\n", normalizePath(auc_summary_file)))
cat(sprintf("  - %s\n", normalizePath(net_benefit_file)))

cat("\nSample counts:\n")
for (cohort_name in names(anchor_cohorts)) {
  cat(sprintf("  %s: Total = %d, Used for ROC = %d\n", 
              cohort_name,
              gene_coverage_list[[cohort_name]]$n_samples_total,
              gene_coverage_list[[cohort_name]]$n_samples_used_for_ROC))
}

cat("\nAUC results:\n")
for (cohort_name in names(anchor_cohorts)) {
  result <- auc_results[[cohort_name]]
  cat(sprintf("  %s: AUC = %.3f (95%% CI: %.3f-%.3f), Flipped: %s\n", 
              cohort_name,
              result$AUC,
              result$CI_lower,
              result$CI_upper,
              result$flipped))
}

# Restore original theme
theme_set(old_theme)

# Memory cleanup
rm(list = ls(pattern = "^(fig|anchor|cohort|expr|score|sample|locked|weight|roc|dca|net|gene_coverage)", all.names = TRUE))
gc()

# Close logging
sink()
