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
target_genes <- locked_weights$Gene_Symbol
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
  scores <- calculate_score(expr_data, locked_weights)
  
  # Merge scores with sample info
  scores_with_info <- scores %>%
    left_join(sample_info, by = "sample_id")
  
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

# Check for ln(IgE) column
if (!"ln_IgE" %in% colnames(gse65204_data)) {
  if ("IgE" %in% colnames(gse65204_data)) {
    # Assume IgE is not logged and log it
    gse65204_data <- gse65204_data %>%
      mutate(ln_IgE = log(IgE + 1))
    cat("  Calculated ln(IgE) from IgE column\n")
  } else {
    stop("Neither 'ln_IgE' nor 'IgE' column found in GSE65204 sample info")
  }
}

# Calculate tertiles for coloring
gse65204_data <- gse65204_data %>%
  mutate(tertile = ntile(score, 3)) %>%
  mutate(tertile_label = case_when(
    tertile == 1 ~ "Low",
    tertile == 2 ~ "Medium",
    tertile == 3 ~ "High"
  ))

# Calculate Spearman correlation
cor_result <- cor.test(gse65204_data$score, gse65204_data$ln_IgE, method = "spearman")
spearman_rho <- cor_result$estimate
spearman_p <- cor_result$p.value

# Store Spearman stats
spearman_stats <- data.frame(
  rho = spearman_rho,
  p_value = spearman_p,
  n = nrow(gse65204_data)
)
spearman_file <- file.path(tables_dir, "Fig4A_spearman_stats.csv")
write_csv(spearman_stats, spearman_file)
cat(sprintf("Saved Spearman stats: %s\n", basename(spearman_file)))

# Create scatter plot
fig4a <- ggplot(gse65204_data, aes(x = score, y = ln_IgE, color = tertile_label)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "loess", color = "black", se = TRUE, alpha = 0.2) +
  scale_color_manual(values = c("Low" = "#4CB5F5", "Medium" = "#FFBB00", "High" = "#FF420E")) +
  labs(
    title = "GSE65204: Signature Score vs ln(IgE)",
    subtitle = sprintf("Spearman rho = %.3f, p = %.3e", spearman_rho, spearman_p),
    x = "11-Gene Signature Score",
    y = "ln(IgE)",
    color = "Score Tertile"
  ) +
  theme_allergy()

# Save Fig4A
fig4a_output <- file.path(figures_dir, "Fig4A_GSE65204_score_vs_lnIgE")
save_trae_figure(fig4a, fig4a_output, width = 180, height = 120, dpi = 300, format = c("pdf", "png"))

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
  
  # Create binary truth labels based on tertiles (drop middle group)
  if (cohort_info$name == "GSE65204") {
    # IgE high vs low
    data <- data %>%
      mutate(tertile = ntile(ln_IgE, 3)) %>%
      filter(tertile != 2) %>%
      mutate(truth = ifelse(tertile == 3, 1, 0))
  } else if (cohort_info$name == "GSE201955") {
    # FeNO high vs low
    if (!"FeNO" %in% colnames(data)) {
      stop("'FeNO' column not found in GSE201955 sample info")
    }
    data <- data %>%
      mutate(tertile = ntile(FeNO, 3)) %>%
      filter(tertile != 2) %>%
      mutate(truth = ifelse(tertile == 3, 1, 0))
  } else if (cohort_info$name == "GSE45111") {
    # Sputum eos high vs low
    # Check for existing phenotype column
    if ("sputum_eos_phenotype" %in% colnames(data)) {
      # Use existing phenotype
      data <- data %>%
        filter(sputum_eos_phenotype %in% c("High", "Low")) %>%
        mutate(truth = ifelse(sputum_eos_phenotype == "High", 1, 0))
    } else if ("sputum_eos_percent" %in% colnames(data)) {
      # Use continuous eos% and create tertiles
      data <- data %>%
        mutate(tertile = ntile(sputum_eos_percent, 3)) %>%
        filter(tertile != 2) %>%
        mutate(truth = ifelse(tertile == 3, 1, 0))
    } else {
      stop("Neither 'sputum_eos_phenotype' nor 'sputum_eos_percent' column found in GSE45111 sample info")
    }
  }
  
  # Update gene coverage sample count
  gene_coverage_list[[cohort_name]]$n_samples_used_for_ROC <- nrow(data)
  
  # Calculate ROC
  roc_result <- roc(response = data$truth, predictor = data$score, direction = "auto")
  auc_value <- auc(roc_result)
  
  # Calculate 95% CI using bootstrap
  ci_result <- ci.auc(roc_result, method = "bootstrap", n.boot = 2000)
  ci_lower <- ci_result[1]
  ci_upper <- ci_result[3]
  
  # Check if predictor was flipped
  predictor_flipped <- ifelse(roc_result$direction == "<", "YES", "NO")
  
  # Store AUC results
  auc_results[[cohort_name]] <- list(
    cohort = cohort_info$name,
    truth = cohort_info$anchor,
    n_high = sum(data$truth == 1),
    n_low = sum(data$truth == 0),
    AUC = auc_value,
    CI_lower = ci_lower,
    CI_upper = ci_upper,
    flipped = predictor_flipped,
    gene_coverage = paste(gene_coverage_list[[cohort_name]]$n_genes_available, "/", gene_coverage_list[[cohort_name]]$n_genes_total, sep = "")
  )
  
  cat(sprintf("  %s: AUC = %.3f (95%% CI: %.3f-%.3f), Flipped: %s\n", 
              cohort_info$name, auc_value, ci_lower, ci_upper, predictor_flipped))
}

# Update gene coverage table with sample counts
gene_coverage_df <- bind_rows(gene_coverage_list)
write_csv(gene_coverage_df, gene_coverage_file)

# Create AUC summary table
auc_summary_df <- bind_rows(auc_results) %>%
  mutate(CI = paste0("(", round(CI_lower, 3), ", ", round(CI_upper, 3), ")")) %>%
  select(cohort, truth, n_high, n_low, AUC, CI, flipped, gene_coverage)

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
    data <- data %>%
      mutate(tertile = ntile(ln_IgE, 3)) %>%
      filter(tertile != 2) %>%
      mutate(truth = ifelse(tertile == 3, 1, 0))
  } else if (cohort_info$name == "GSE201955") {
    data <- data %>%
      mutate(tertile = ntile(FeNO, 3)) %>%
      filter(tertile != 2) %>%
      mutate(truth = ifelse(tertile == 3, 1, 0))
  } else if (cohort_info$name == "GSE45111") {
    if ("sputum_eos_phenotype" %in% colnames(data)) {
      data <- data %>%
        filter(sputum_eos_phenotype %in% c("High", "Low")) %>%
        mutate(truth = ifelse(sputum_eos_phenotype == "High", 1, 0))
    } else if ("sputum_eos_percent" %in% colnames(data)) {
      data <- data %>%
        mutate(tertile = ntile(sputum_eos_percent, 3)) %>%
        filter(tertile != 2) %>%
        mutate(truth = ifelse(tertile == 3, 1, 0))
    }
  }
  
  # Calculate ROC
  roc_result <- roc(response = data$truth, predictor = data$score, direction = "auto")
  
  # Get AUC and CI
  auc_value <- auc(roc_result)
  ci_result <- ci.auc(roc_result, method = "bootstrap", n.boot = 2000)
  
  # Create ROC data frame for plotting
  roc_df <- data.frame(
    specificity = 1 - roc_result$specificities,
    sensitivity = roc_result$sensitivities
  )
  
  # Get gene coverage
  gene_cov <- auc_results[[cohort_name]]$gene_coverage
  
  # Create plot
  plot_title <- sprintf("%s\n%s\nAUC = %.3f (95%% CI: %.3f-%.3f)\nGene coverage: %s",
                        cohort_info$name,
                        cohort_info$anchor,
                        auc_value,
                        ci_result[1],
                        ci_result[3],
                        gene_cov)
  
  roc_plot <- ggplot(roc_df, aes(x = specificity, y = sensitivity)) +
    geom_line(linewidth = 1.2, color = trae_colors_cohort[cohort_name]) +
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
    tag_levels = "A"
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

# Logistic calibration to get predicted probabilities
logit_model <- glm(truth ~ score, data = gse65204_roc_data, family = binomial)
gse65204_roc_data$pred_prob <- predict(logit_model, type = "response")

# Perform DCA
dca_result <- dca(truth ~ pred_prob, data = gse65204_roc_data)

# Create DCA plot
fig4c <- plot(dca_result) +
  labs(
    title = "Decision Curve Analysis (GSE65204)",
    subtitle = "11-Gene Signature vs No Treatment vs All Treatment"
  ) +
  theme_allergy()

# Save Fig4C
fig4c_output <- file.path(figures_dir, "Fig4C_DCA_GSE65204")
save_trae_figure(fig4c, fig4c_output, width = 180, height = 120, dpi = 300, format = c("pdf", "png"))

# Extract net benefit at key thresholds
thresholds <- c(0.2, 0.3, 0.4, 0.5)
net_benefit_df <- data.frame(threshold = thresholds)

for (threshold in thresholds) {
  nb_result <- net.benefit(dca_result, thresholds = threshold)
  net_benefit_df[net_benefit_df$threshold == threshold, "model"] <- nb_result$net.benefit[nb_result$label == "pred_prob"]
  net_benefit_df[net_benefit_df$threshold == threshold, "all"] <- nb_result$net.benefit[nb_result$label == "All"]
  net_benefit_df[net_benefit_df$threshold == threshold, "none"] <- nb_result$net.benefit[nb_result$label == "None"]
}

net_benefit_file <- file.path(tables_dir, "Fig4C_net_benefit_key_thresholds.csv")
write_csv(net_benefit_df, net_benefit_file)
cat(sprintf("Saved net benefit table: %s\n", basename(net_benefit_file)))

# =========================================================================
# Fig4D: DGIdb Drug-Gene Network
# =========================================================================
cat("Generating Fig4D: DGIdb Drug-Gene Network...\n")

# Load required packages for network analysis
suppressPackageStartupMessages({
  library(igraph)
  library(ggraph)
  library(ggrepel)
})

# Read locked genes from weights file (lines 2-12, 11 genes total)
locked_genes <- locked_weights$Gene
# Ensure we only use 11 genes
if (length(locked_genes) > 11) {
  locked_genes <- locked_genes[1:11]
  cat("Note: Using only first 11 genes from locked_weights.csv\n")
}
cat(sprintf("Using %d locked genes for DGIdb analysis\n", length(locked_genes)))
cat("Gene symbols:", paste(locked_genes, collapse = ", "), "\n")

# DGIdb data path
dgidb_dir <- file.path(base_dir, "data", "raw", "drug_network")
in_interactions <- file.path(dgidb_dir, "interactions.tsv")
in_drugs        <- file.path(dgidb_dir, "drugs.tsv")
in_genes        <- file.path(dgidb_dir, "genes.tsv")

# Check if DGIdb data is available
if (!file.exists(in_interactions) || !file.exists(in_drugs) || !file.exists(in_genes)) {
  cat("  DGIdb data files not found, creating placeholder plot.\n")
  cat("  Expected files in:", dgidb_dir, "\n")
  cat("  Required files: interactions.tsv, drugs.tsv, genes.tsv\n")
  
  fig4d <- ggplot() +
    geom_text(aes(x = 1, y = 1, label = "DGIdb Drug-Gene Interaction Plot\n(placeholder)"), size = 10) +
    theme_void()
} else {
  cat("  DGIdb data found, building network...\n")
  
  # Read DGIdb TSV files
  interactions_raw <- readr::read_tsv(in_interactions, show_col_types = FALSE)
  drugs_raw        <- readr::read_tsv(in_drugs, show_col_types = FALSE)
  genes_raw        <- readr::read_tsv(in_genes, show_col_types = FALSE)
  
  cat(sprintf("  - interactions: %d rows\n", nrow(interactions_raw)))
  cat(sprintf("  - drugs:        %d rows\n", nrow(drugs_raw)))
  cat(sprintf("  - genes:        %d rows\n", nrow(genes_raw)))
  
  # Standardize interaction data
  int2 <- interactions_raw %>%
    transmute(
      gene_name = as.character(gene_name),
      drug_name = as.character(drug_name),
      interaction_type = if ("interaction_type" %in% names(interactions_raw)) {
        as.character(interaction_type)
      } else if ("interaction_types" %in% names(interactions_raw)) {
        as.character(interaction_types)
      } else {
        NA_character_
      },
      interaction_score = if ("interaction_score" %in% names(interactions_raw)) {
        suppressWarnings(as.numeric(interaction_score))
      } else if ("score" %in% names(interactions_raw)) {
        suppressWarnings(as.numeric(score))
      } else {
        NA_real_
      },
      source_db = if ("interaction_source_db_name" %in% names(interactions_raw)) {
        as.character(interaction_source_db_name)
      } else if ("source_db_name" %in% names(interactions_raw)) {
        as.character(source_db_name)
      } else {
        NA_character_
      },
      approved_int = if ("approved" %in% names(interactions_raw)) {
        case_when(
          approved == "TRUE" ~ TRUE,
          approved == "FALSE" ~ FALSE,
          TRUE ~ NA
        )
      } else {
        NA
      },
      drug_concept_id = if ("drug_concept_id" %in% names(interactions_raw)) {
        as.character(drug_concept_id)
      } else if ("drug_id" %in% names(interactions_raw)) {
        as.character(drug_id)
      } else {
        NA_character_
      }
    )
  
  # Standardize drug data
  dr2 <- drugs_raw %>%
    filter(!is.na(concept_id) & concept_id != "NULL" & concept_id != "") %>%
    transmute(
      concept_id = as.character(concept_id),
      drug_name_drugs = as.character(drug_name),
      approved_drugs = if ("approved" %in% names(drugs_raw)) {
        case_when(
          approved == "TRUE" ~ TRUE,
          approved == "FALSE" ~ FALSE,
          TRUE ~ NA
        )
      } else NA
    ) %>%
    group_by(concept_id) %>%
    summarise(
      drug_name_drugs = dplyr::first(na.omit(drug_name_drugs)),
      approved_drugs = ifelse(all(is.na(approved_drugs)), NA, any(approved_drugs, na.rm = TRUE)),
      .groups = "drop"
    )
  
  # Filter interactions for locked genes
  int_locked <- int2 %>%
    filter(gene_name %in% locked_genes)
  
  cat(sprintf("  - interactions kept: %d rows\n", nrow(int_locked)))
  
  if (nrow(int_locked) == 0) {
    cat("  No interactions found for locked genes, creating placeholder plot.\n")
    fig4d <- ggplot() +
      geom_text(aes(x = 1, y = 1, label = "DGIdb Drug-Gene Interaction Plot\n(No interactions found)"), size = 10) +
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
        interaction_score = ifelse(is.na(interaction_score), 0, interaction_score),
        source_db = ifelse(is.na(source_db) | source_db == "", "UnknownSource", source_db),
        interaction_type = ifelse(is.na(interaction_type) | interaction_type == "", "UNKNOWN", interaction_type)
      ) %>%
      select(
        gene_name,
        drug_name = drug_name_clean,
        approved = approved_final,
        interaction_type,
        interaction_score,
        source_db,
        drug_concept_id
      )
    
    # Plot filtering parameters
    topN_per_gene <- 3
    whitelist_drugs <- c(
      "Dupilumab", "Mepolizumab", "Benralizumab", "Omalizumab",
      "Dexamethasone", "Prednisone", "Hydrocortisone", "Budesonide", "Fluticasone",
      "Montelukast", "Zileuton",
      "Loratadine", "Cetirizine"
    )
    blacklist_drugs <- c(
      "HEXACHLOROPHENE", "SILVER SULFADIAZINE", "PYRITHIONE ZINC", "SELENIUM SULFIDE",
      "CHOLECALCIFEROL",
      "LINOLEIC ACID", "LINOLENIC ACID", "OLEIC ACID",
      "NULL", "NA"
    )
    
    # Create filtered edge list for plotting
    plot_df <- int_locked_enriched %>%
      group_by(gene_name) %>%
      arrange(desc(approved), desc(interaction_score), drug_name) %>%
      mutate(rank_within_gene = row_number()) %>%
      ungroup()
    
    # Filter logic
    plot_df_filtered <- plot_df %>%
      filter(
        approved | 
        rank_within_gene <= topN_per_gene |
        tolower(drug_name) %in% tolower(whitelist_drugs)
      ) %>%
      filter(
        !tolower(drug_name) %in% tolower(blacklist_drugs),
        !str_detect(tolower(drug_name), "zinc|silver"),
        !str_detect(tolower(drug_name), "acid") | str_detect(tolower(drug_name), "valproic acid")
      ) %>%
      distinct(gene_name, drug_name, .keep_all = TRUE)
    
    cat(sprintf("  - Filtered to %d unique gene-drug pairs\n", nrow(plot_df_filtered)))
    
    if (nrow(plot_df_filtered) == 0) {
      cat("  No interactions after filtering, creating placeholder plot.\n")
      fig4d <- ggplot() +
        geom_text(aes(x = 1, y = 1, label = "DGIdb Drug-Gene Interaction Plot\n(No interactions after filtering)"), size = 10) +
        theme_void()
    } else {
      # Build network
      edges <- plot_df_filtered %>%
        transmute(
          from = gene_name,
          to = drug_name,
          approved = approved,
          interaction_type = interaction_type,
          score = interaction_score
        )
      
      g <- igraph::graph_from_data_frame(edges, directed = FALSE)
      
      V(g)$node_type <- ifelse(V(g)$name %in% locked_genes, "gene", "drug")
      V(g)$degree <- igraph::degree(g)
      
      # Drug approved mapping
      drug_approved_map <- edges %>%
        group_by(to) %>%
        summarise(any_approved = any(approved), .groups = "drop") %>%
        deframe()
      
      # Node styling
      for (v in V(g)$name) {
        if (V(g)[v]$node_type == "gene") {
          V(g)[v]$color <- "#E64B35"  # Gene: dark red
          V(g)[v]$size <- 4.0         # Gene node size
        } else {
          is_whitelist <- tolower(v) %in% tolower(whitelist_drugs) & !tolower(v) %in% c("loratadine")
          is_appr <- FALSE
          if (v %in% names(drug_approved_map)) is_appr <- drug_approved_map[[v]]
          
          if (is_whitelist) {
            V(g)[v]$color <- "#3C5488"  # Key drug: dark blue
            V(g)[v]$size <- 3.0         # Key drug size
          } else {
            V(g)[v]$color <- "#91D1C2"   # Regular drug: light blue
            V(g)[v]$size <- 2.5         # Regular drug size
          }
        }
      }
      
      # Layout
      set.seed(42)
      lay <- igraph::layout_with_fr(g, niter = 1500, start.temp = 0.1)
      
      # Node data for plotting
      node_data <- tibble(
        name = V(g)$name,
        node_type = V(g)$node_type,
        degree = V(g)$degree,
        color = V(g)$color,
        size = V(g)$size,
        x = lay[, 1],
        y = lay[, 2],
        is_whitelist = tolower(name) %in% tolower(whitelist_drugs)
      ) %>%
        left_join(
          tibble(drug_name = names(drug_approved_map), any_approved = as.logical(drug_approved_map)),
          by = c("name" = "drug_name")
        ) %>%
        mutate(any_approved = ifelse(is.na(any_approved), FALSE, any_approved)) %>%
        mutate(
          label_size = case_when(
            node_type == "gene" ~ 3.5,          # Gene: 3.5pt
            is_whitelist ~ 3.0,                  # Key drug: 3.0pt
            TRUE ~ 2.0                           # Regular drug: 2.0pt
          ),
          label_color = case_when(
            node_type == "gene" ~ "black",    # Gene: black
            is_whitelist & !tolower(name) %in% c("loratadine") ~ "#3C5488",  # Key drug: dark blue
            tolower(name) == "loratadine" ~ "grey50",  # Loratadine: grey
            TRUE ~ "grey50"                   # Regular drug: grey
          ),
          label_fontface = case_when(
            node_type == "gene" ~ "bold",     # Gene: bold
            is_whitelist & !tolower(name) %in% c("loratadine") ~ "bold",  # Key drug: bold
            tolower(name) == "loratadine" ~ "plain",  # Loratadine: plain
            TRUE ~ "plain"                    # Regular drug: plain
          ),
          label = case_when(
            node_type == "gene" ~ name,                  # All genes
            is_whitelist ~ name,                          # Whitelist drugs
            TRUE ~ ""
          )
        )
      
      # Edge styling
      edge_data <- edges %>%
        mutate(
          is_highlight = tolower(to) %in% tolower(c("Dupilumab", "Dexamethasone")),
          edge_color = ifelse(is_highlight, "grey40", "grey80"),
          edge_width = ifelse(is_highlight, 0.6, 0.2)
        )
      
      # Update graph edge attributes
      g <- g %>%
        set_edge_attr("color", value = edge_data$edge_color) %>%
        set_edge_attr("width", value = edge_data$edge_width)
      
      # Create plot
      fig4d <- ggraph(g, layout = lay) +
        geom_edge_link(aes(color = color, linewidth = width)) +
        geom_node_point(aes(color = color, size = size), alpha = 0.95) +
        geom_node_text(
          data = node_data %>% filter(label != ""),
          aes(
            x = x, 
            y = y, 
            label = label,
            size = label_size,
            color = label_color,
            fontface = label_fontface
          ),
          repel = TRUE,
          max.overlaps = Inf,
          box.padding = unit(0.6, "lines"),
          point.padding = unit(0.5, "lines"),
          segment.color = "grey50",
          segment.size = 0.5,
          force = 0.5,
          bg.color = "white",
          bg.r = 0.1
        ) +
        scale_size_identity() +
        scale_color_identity() +
        guides(
          color = guide_legend(
            title = "Node Type",
            override.aes = list(
              size = c(4.0, 3.0, 2.5),
              color = c("#E64B35", "#3C5488", "#91D1C2")
            )
          ),
          size = "none"
        ) +
        theme_void() +
        labs(
          title = "DGIdb Drug–Gene Interaction Network (11-Gene Signature)",
          subtitle = paste0("Edges filtered: DGIdb approved flag + top ", topN_per_gene, " drugs per gene"),
          caption = "Data source: DGIdb local TSV (DGIdb approved flag does not necessarily mean FDA approval)"
        ) +
        theme(
          plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
          plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey30"),
          plot.caption = element_text(size = 8, color = "grey50", hjust = 1),
          plot.margin = margin(10, 10, 10, 10)
        )
    }
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
  description = "DGIdb Drug-Gene Interaction Plot (placeholder)",
  script_name = "10_fig4_anchors_trans.R",
  resolution = "300 dpi",
  notes = "Placeholder - data not yet implemented"
)

# Log tables
log_output(
  file_path = gene_coverage_file,
  file_type = "table",
  table_id = "Fig4_gene_coverage",
  description = "Gene coverage by anchor cohort",
  script_name = "10_fig4_anchors_trans.R"
)

log_output(
  file_path = spearman_file,
  file_type = "table",
  table_id = "Fig4A_spearman",
  description = "Spearman correlation statistics for GSE65204",
  script_name = "10_fig4_anchors_trans.R"
)

log_output(
  file_path = auc_summary_file,
  file_type = "table",
  table_id = "Fig4B_auc_summary",
  description = "AUC summary for three clinical anchors",
  script_name = "10_fig4_anchors_trans.R"
)

log_output(
  file_path = net_benefit_file,
  file_type = "table",
  table_id = "Fig4C_net_benefit",
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
