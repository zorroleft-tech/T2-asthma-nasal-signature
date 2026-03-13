#!/usr/bin/env Rscript

# =========================================================================
# Script: 16_supp_figS2_landscape.R
# Purpose: Generate Figure S2B and Table S9 coverage matrix
#   - Read platform-wide TableS9 file for all cohorts except GSE152004
#   - Add GSE152004 as full coverage (11/11)
#   - Generate coverage tile plot and gene coverage matrix
# Author: Zuo
# Date: 2026-03-01
#
# Inputs:
#   - output/supplement/TableS9_platform_wide.csv  # Platform-wide coverage data
#
# Outputs:
#   - output/supplement/FigS2B_coverage_tile.pdf  # Coverage tile plot
#   - output/tables_main/TableS9_Gene_Coverage_Matrix.csv  # Gene coverage matrix
# =========================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(ggplot2)
  library(glue)
})

# ----------------------------
# 0) Paths
# ----------------------------
base_dir    <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
dir_in      <- file.path(base_dir, "output", "supplement")          # <= 实际路径
dir_out_fig <- file.path(base_dir, "output", "supplement")
dir_out_tab <- file.path(base_dir, "output", "tables_main")
dir_logs    <- file.path(base_dir, "output", "logs")

# ----------------------------
# 0.1) Load theme setup
# ----------------------------
theme_setup_file <- file.path(base_dir, "R", "00_theme_setup.R")
if (file.exists(theme_setup_file)) {
  source(theme_setup_file)
  if (exists("set_plot_theme")) set_plot_theme()
}

if (!dir.exists(dir_out_fig)) dir.create(dir_out_fig, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(dir_out_tab)) dir.create(dir_out_tab, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(dir_logs))    dir.create(dir_logs, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(dir_logs, "16_supp_figS2_landscape.log")
log_line <- function(...) {
  msg <- glue(...)
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

log_line("=== simplified 16_supp_figS2_landscape.R started at {format(Sys.time(), '%Y-%m-%d %H:%M:%S')} ===")

# ----------------------------
# 1) Locked genes (must match manuscript)
# ----------------------------
locked_genes <- c(
  "SERPINB2","CLCA1","CCR3","HDC","IL4","CPA3","MUC5B","MUC5AC","IL13","CCL24","POSTN"
)

# ----------------------------
# 2) Read TableS9 platform-wide file (10 datasets)
# ----------------------------
in_file <- file.path(dir_in, "TableS9_platform_wide.csv")  # <= 把文件放这里，或改路径
stopifnot(file.exists(in_file))

tab <- read_csv(in_file, show_col_types = FALSE)

# Expect columns: Dataset, Platform, Role, 11 genes..., Coverage_X_of_11, Missing_Genes
stopifnot(all(c("Dataset","Platform","Role") %in% names(tab)))
stopifnot(all(locked_genes %in% names(tab)))

# ----------------------------
# 3) Add GSE152004 (discovery): full coverage 11/11
# ----------------------------
gse152004_row <- tibble(
  Dataset = "GSE152004",
  Platform = "GPL11154",         # RNA-seq; platform字段只是展示用
  Role = "Discovery",
  SERPINB2 = 1, CLCA1 = 1, CCR3 = 1, HDC = 1, IL4 = 1, CPA3 = 1, MUC5B = 1, MUC5AC = 1, IL13 = 1, CCL24 = 1, POSTN = 1,
  Coverage_X_of_11 = 11,
  Missing_Genes = "None"
)

tab_all <- bind_rows(gse152004_row, tab)

# 强制按我们定义的 cohort 顺序（你也可按 Role 排序）
cohort_order <- c("GSE152004", tab$Dataset)
tab_all$Dataset <- factor(tab_all$Dataset, levels = unique(cohort_order))

# ----------------------------
# 4) Build long format for plotting
# ----------------------------
plot_df <- tab_all %>%
  select(Dataset, all_of(locked_genes)) %>%
  pivot_longer(cols = all_of(locked_genes), names_to = "Gene", values_to = "Present") %>%
  mutate(
    Gene = factor(Gene, levels = rev(locked_genes)),
    Status = if_else(Present == 1, "Present", "Missing"),
    Status = factor(Status, levels = c("Present", "Missing"))
  )

# ----------------------------
# 5) Write Table S9 matrix (gene x cohort; Present/Missing)
# ----------------------------
table_s9_out <- plot_df %>%
  select(Gene, Dataset, Status) %>%
  pivot_wider(names_from = Dataset, values_from = Status) %>%
  arrange(match(Gene, locked_genes))

out_table <- file.path(dir_out_tab, "TableS9_Gene_Coverage_Matrix.csv")
write_csv(table_s9_out, out_table)
log_line("Wrote: {out_table}")

# ----------------------------
# 6) Plot Fig S2B (tile)
# ----------------------------
status_colors <- c("Present" = "#2C7FB8", "Missing" = "grey75")

p <- ggplot(plot_df, aes(x = Dataset, y = Gene, fill = Status)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_manual(values = status_colors, drop = FALSE) +
  geom_text(
    data = plot_df %>% filter(Status == "Missing"),
    aes(label = "×"),
    color = "grey30",
    size = 5
  ) +
  labs(
    x = NULL, y = NULL, fill = NULL,
    title = "Fig S2B. Locked 11-gene coverage across 11 datasets (including discovery GSE152004)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 12)
  )

out_fig <- file.path(dir_out_fig, "FigS2B_coverage_tile.pdf")
ggsave(out_fig, p, width = 11, height = 5, useDingbats = FALSE)
log_line("Wrote: {out_fig}")

log_line("=== finished at {format(Sys.time(), '%Y-%m-%d %H:%M:%S')} ===")