# =========================================================================
# Script: 18_supp_figS4_architecture.R
# Purpose: Generate Supplementary Figure S4 (Signature architecture & gene drivers)
# Author: TRAE Team + Notion QA
# Date: 2026-03-02
#
# Inputs:
#   - data/derived/locked_weights.csv  # Locked gene weights
#   - data/processed/expr_gse65204&GPL14550__FULL.rds  # GSE65204 expression data
#   - data/processed/phase3_outputs/GSE65204_pheno_raw.rds  # GSE65204 phenotype data
#   - data/processed/expr_*&*__FULL.rds  # Replication cohorts expression data
#   - R/03_index_logger.R  # Logging script
#
# Outputs:
#   - output/supplement/FigS4A_GSE65204_genewise_IgE_forest.pdf  # Gene-wise IgE correlation forest plot
#   - output/supplement/FigS4B_cross_cohort_modular_consistency.pdf  # Cross-cohort modular consistency heatmap
#   - output/supplement/FigS4A_Source_Data.csv  # Source data for Fig S4A
#   - output/supplement/FigS4B_Source_Data.csv  # Source data for Fig S4B
#   - output/logs/18_supp_figS4_architecture.log  # Log file
# =========================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
  library(ggplot2)
})

base_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
processed_dir <- file.path(base_dir, "data", "processed")
derived_dir <- file.path(base_dir, "data", "derived")
fig_dir <- file.path(base_dir, "output", "supplement")
tbl_dir <- file.path(base_dir, "output", "supplement")
logs_dir <- file.path(base_dir, "output", "logs")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tbl_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

# 安全日志记录机制 (Notion 建议)
log_file <- file.path(logs_dir, "18_supp_figS4_architecture.log")
sink(log_file, append = FALSE, split = TRUE)
on.exit({ if(sink.number() > 0) sink() }, add = TRUE)

cat("========================================================\n")
cat("=== Generating Figure S4 (Architecture & Drivers) ===\n")
cat("========================================================\n\n")

# Load theme setup
theme_setup_file <- file.path(base_dir, "R", "00_theme_setup.R")
if (file.exists(theme_setup_file)) {
  source(theme_setup_file)
  if (exists("set_plot_theme")) set_plot_theme()
}

source(file.path(base_dir, "R", "03_index_logger.R"))

# 1. Load Locked Genes & Define Modules
weights_file <- file.path(derived_dir, "locked_weights.csv")
if (!file.exists(weights_file)) stop("locked_weights.csv not found!")
sig_genes <- read_csv(weights_file, show_col_types = FALSE)$Gene

# Unified 3-tier Module definition (Notion 建议)
mod_def <- tibble(
  Gene = sig_genes,
  Module = case_when(
    Gene %in% c("MUC5AC", "MUC5B", "CLCA1", "SERPINB2", "POSTN") ~ "Epithelial",
    Gene %in% c("IL4", "IL13", "CCR3") ~ "Immune",
    Gene %in% c("CPA3", "HDC", "CCL24") ~ "Effector",
    TRUE ~ "Other"
  )
) %>% mutate(Module = factor(Module, levels = c("Epithelial", "Immune", "Effector")))

# Structural sub-modules definition (11-gene aligned)
gene_sets <- list(
  full_11 = sig_genes,
  overlap_4 = c("SERPINB2", "MUC5AC", "CPA3", "HDC"),
  unique_7 = setdiff(sig_genes, c("SERPINB2", "MUC5AC", "CPA3", "HDC")),
  epithelial = c("MUC5AC", "MUC5B", "CLCA1", "SERPINB2", "POSTN"),
  immune = c("IL4", "IL13", "CCR3"),
  effector = c("CPA3", "HDC", "CCL24")
)

# ----------------------------------------------------------------------------
# PART A: Gene-wise correlation with serum IgE in GSE65204 (Forest Plot)
# ----------------------------------------------------------------------------
cat("[Part A] Computing Genewise IgE Correlations (B=2000 Bootstrap)...\n")

expr_65204_path <- file.path(processed_dir, "expr_gse65204&GPL14550__FULL.rds")
pheno_65204_path <- file.path(processed_dir, "phase3_outputs", "GSE65204_pheno_raw.rds")

if (file.exists(expr_65204_path) && file.exists(pheno_65204_path)) {
  expr_mat <- readRDS(expr_65204_path)
  pheno_df <- readRDS(pheno_65204_path)
  
  # 从raw_characteristics_all列中解析IgE值
  cat("  Extracting IgE values from raw_characteristics_all column...\n")
  
  # 解析raw_characteristics_all列，提取lnige值
  pheno_df <- pheno_df %>%
    mutate(lnige = as.numeric(str_extract(raw_characteristics_all, "lnige: ([0-9.]+)") %>% str_replace("lnige: ", "")))
  
  ige_col <- "lnige"
  sample_id_col <- "std_sample_id"
  
  if (!is.na(ige_col) && ige_col %in% names(pheno_df)) {
    valid_pheno <- pheno_df %>% filter(!is.na(.data[[ige_col]]))
    common <- intersect(colnames(expr_mat), valid_pheno[[sample_id_col]])
    
    cat(sprintf("  Found %d common samples\n", length(common)))
    
    if (length(common) > 0) {
      expr_sub <- expr_mat[, common, drop=FALSE]
      valid_pheno <- valid_pheno[match(common, valid_pheno[[sample_id_col]]), ]
      
      # IgE值已经是对数形式
      ige_vals <- as.numeric(valid_pheno[[ige_col]])
      cat(sprintf("  IgE values range: %.2f - %.2f\n", min(ige_vals, na.rm = TRUE), max(ige_vals, na.rm = TRUE)))
      
      avail_genes <- intersect(sig_genes, rownames(expr_sub))
      cat(sprintf("  GSE65204 available genes: %d/11\n", length(avail_genes)))
      
      # Bootstrap CI Function
      boot_spearman_ci <- function(x, y, B = 2000) {
        n <- length(x)
        rho_boot <- replicate(B, {
          idx <- sample.int(n, replace = TRUE)
          cor(x[idx], y[idx], method = "spearman")
        })
        quantile(rho_boot, probs = c(0.025, 0.975), na.rm = TRUE)
      }
      
      res_A <- list()
      set.seed(2026)
      
      for (g in avail_genes) {
        g_expr <- as.numeric(expr_sub[g, ])
        ct <- suppressWarnings(cor.test(g_expr, ige_vals, method = "spearman"))
        ci <- boot_spearman_ci(g_expr, ige_vals)
        res_A[[g]] <- data.frame(
          Gene = g, rho = ct$estimate, p = ct$p.value,
          ci_low = ci[1], ci_high = ci[2], stringsAsFactors = FALSE
        )
      }
      
      df_A <- bind_rows(res_A) %>%
        mutate(q = p.adjust(p, method = "BH")) %>%
        left_join(mod_def, by = "Gene") %>%
        arrange(desc(rho)) %>%
        mutate(Gene = factor(Gene, levels = rev(Gene)))
      
      # Plot S4A
      p_s4a <- ggplot(df_A, aes(y = Gene, x = rho, color = Module)) +
        geom_vline(xintercept = 0, color = "gray30", linewidth = 0.5) +
        geom_vline(xintercept = 0.3, linetype = "dashed", color = "gray60", linewidth = 0.4) +
        geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.2, linewidth = 0.8) +
        geom_point(size = 3) +
        geom_text(aes(x = -0.55, label = sprintf("q = %.3g", q)), color = "black", hjust = 0, size = 3.5) +
        scale_color_manual(values = c("Epithelial" = "#E64B35", "Immune" = "#4DBBD5", "Effector" = "#00A087")) +
        scale_x_continuous(limits = c(-0.6, 0.8), breaks = seq(-0.6, 0.8, by = 0.2), labels = c("-0.6", "-0.4", "-0.2", "0", "0.2", "0.4", "0.6", "0.8")) +
        facet_grid(Module ~ ., scales = "free_y", space = "free_y") +
        labs(x = "Spearman Correlation (ρ) with ln(IgE)", y = "", title = "A. Gene-wise correlation with serum IgE (GSE65204)", subtitle = "Error bars: 95% CI (2,000 bootstraps)") +
        theme_minimal(base_size = 12) +
        theme(
          legend.position = "none", 
          strip.background = element_rect(fill = "gray90"), 
          strip.text = element_text(face = "bold"),
          axis.text.x = element_text(face="bold", color="black"), 
          axis.text.y = element_text(face="bold", color="black"),
          panel.grid = element_blank(),
          axis.line = element_line(color = "black"),
          axis.ticks = element_line(color = "black")
        )
      
      ggsave(file.path(fig_dir, "FigS4A_GSE65204_genewise_IgE_forest.pdf"), p_s4a, width = 6.5, height = 5)
      cat("  ✓ Saved FigS4A_GSE65204_genewise_IgE_forest.pdf\n")
    } else {
      cat("  ⚠ No common samples found. Skipping Part A.\n")
    }
  } else {
    cat("  ⚠ IgE column not found in GSE65204. Skipping Part A.\n")
  }
}

# ----------------------------------------------------------------------------
# PART B: Cross-cohort Modular Consistency (Replication-only Proxy)
# ----------------------------------------------------------------------------
cat("\n[Part B] Computing Cross-cohort Modular Consistency...\n")

rep_cohorts <- c("GSE103166", "GSE118761", "GSE43696", "GSE40888", "GSE123750", "GSE230048", "GSE115770")
res_B <- list()

for (cohort in rep_cohorts) {
  # 构建匹配实际文件格式的路径
  expr_files <- list.files(processed_dir, pattern = paste0("expr_", tolower(cohort), "&.*__FULL.rds"), ignore.case = TRUE)
  cat(sprintf("  Checking cohort %s...\n", cohort))
  if (length(expr_files) > 0) {
    expr_path <- file.path(processed_dir, expr_files[1])
    cat(sprintf("  Found file: %s\n", expr_path))
    if (file.exists(expr_path)) {
      expr_mat <- readRDS(expr_path)
      cat(sprintf("  Loaded expression matrix: %d genes × %d samples\n", nrow(expr_mat), ncol(expr_mat)))
      for (mod_name in names(gene_sets)) {
        genes <- gene_sets[[mod_name]]
        avail <- intersect(genes, rownames(expr_mat))
        cat(sprintf("  Module %s: %d/%d genes available\n", mod_name, length(avail), length(genes)))
        if (length(avail) >= 2) {
          cormat <- suppressWarnings(cor(t(expr_mat[avail, ]), method = "spearman"))
          mean_abs_cor <- mean(abs(cormat[upper.tri(cormat)]), na.rm = TRUE)
          cat(sprintf("  Module %s: MeanAbsCor = %.3f\n", mod_name, mean_abs_cor))
        } else {
          mean_abs_cor <- NA
          cat(sprintf("  Module %s: Insufficient genes (skipping)\n", mod_name))
        }
        res_B[[paste(cohort, mod_name)]] <- data.frame(
          Cohort = cohort, Module = mod_name, MeanAbsCor = mean_abs_cor, stringsAsFactors = FALSE
        )
      }
    } else {
      cat(sprintf("  File not found: %s\n", expr_path))
    }
  } else {
    cat(sprintf("  No expression file found for %s\n", cohort))
  }
}

cat(sprintf("\nTotal results collected: %d\n", length(res_B)))

if (length(res_B) > 0) {
  df_B <- bind_rows(res_B) %>% filter(!is.na(MeanAbsCor))
  
  # Order definition (Notion: 11-gene aligned)
  mod_order <- c("full_11", "overlap_4", "unique_7", "epithelial", "immune", "effector")
  coh_order <- c("GSE115770", "GSE43696", "GSE118761", "GSE103166", "GSE40888", "GSE123750", "GSE230048")
  
  df_B <- df_B %>%
    mutate(Module = factor(Module, levels = mod_order),
           Cohort = factor(Cohort, levels = rev(coh_order)))
  
  p_s4b <- ggplot(df_B, aes(x = Module, y = Cohort, fill = MeanAbsCor)) +
    coord_fixed(ratio = 1) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", MeanAbsCor), color = ifelse(MeanAbsCor > 0.45, "white", "black")), size = 3.5) +
    scale_color_identity() +
    scale_fill_viridis_c(option = "mako", direction = -1, limits = c(0, 0.7), oob = scales::squish, name = "Mean\nAbs Cor") +
    labs(x = "Structural Sub-modules", y = "Replication Cohorts", title = "B. Cross-cohort modular consistency") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face="bold", color="black"), 
          axis.text.y = element_text(face="bold", color="black"), panel.grid = element_blank())
  
  ggsave(file.path(fig_dir, "FigS4B_cross_cohort_modular_consistency.pdf"), p_s4b, width = 7, height = 4.5)
  cat("  ✓ Saved FigS4B_cross_cohort_modular_consistency.pdf\n")
  
  # Separate Source Data output (Notion 建议)
  if(exists("df_A")) {
    write_csv(df_A, file.path(tbl_dir, "FigS4A_Source_Data.csv"))
  }
  write_csv(df_B, file.path(tbl_dir, "FigS4B_Source_Data.csv"))
  
  if(exists("log_output")) {
    log_output(file.path(fig_dir, "FigS4A_GSE65204_genewise_IgE_forest.pdf"), "figure", "Fig S4A", "Gene-wise IgE correlation forest plot (B=2000)", "18_supp_figS4_architecture.R")
    log_output(file.path(fig_dir, "FigS4B_cross_cohort_modular_consistency.pdf"), "figure", "Fig S4B", "Cross-cohort mean absolute correlation heatmap", "18_supp_figS4_architecture.R")
  }
}

cat("\n=== Figure S4 Generation Completed ===\n")
# sink() handled by on.exit()