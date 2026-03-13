# =========================================================================
# Script: 19_supp_figS5_robustness.R
# Purpose: Generate Supplementary Figure S5 (Robustness & Calibration)
#          A. Dose-response (Tertile boxplots via ntile)
#          B. Calibration curve (Logistic Calibration + Stats export)
#          C. Permutation test (Independent label shuffling)
# Author: Zuo
# Date: 2026-03-02
#
# Inputs:
#   - data/derived/locked_weights.csv  # Locked gene weights
#   - data/processed/expr_gse65204&GPL14550__FULL.rds  # GSE65204 expression data
#   - data/processed/phase3_outputs/GSE65204_pheno_raw.rds  # GSE65204 phenotype data
#   - R/03_index_logger.R  # Logging script
#
# Outputs:
#   - output/supplement/FigS5A_Dose_Response.pdf  # Dose-response boxplots
#   - output/supplement/FigS5B_Calibration_Curve.pdf  # Logistic calibration curve
#   - output/supplement/FigS5C_Permutation_Test.pdf  # Permutation test results
#   - output/supplement/FigS5B_calibration_stats.csv  # Calibration statistics
#   - output/supplement/FigS5_Source_Data.csv  # Source data for all panels
#   - output/logs/19_supp_figS5_robustness.log  # Log file
# =========================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(rms)
  library(cowplot)
  library(ggplot2)
})

base_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
processed_dir <- file.path(base_dir, "data", "processed")
derived_dir <- file.path(base_dir, "data", "derived")
supp_dir <- file.path(base_dir, "output", "supplement") 
logs_dir <- file.path(base_dir, "output", "logs")

dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

# 安全日志记录机制
log_file <- file.path(logs_dir, "19_supp_figS5_robustness.log")
sink(log_file, append = FALSE, split = TRUE)
on.exit({ if(sink.number() > 0) sink() }, add = TRUE)

cat("========================================================\n")
cat("=== Generating Figure S5 (Robustness & Calibration) ===\n")
cat("========================================================\n\n")

# Load theme setup
theme_setup_file <- file.path(base_dir, "R", "00_theme_setup.R")
if (file.exists(theme_setup_file)) {
  source(theme_setup_file)
  if (exists("set_plot_theme")) set_plot_theme()
}

source(file.path(base_dir, "R", "03_index_logger.R"))
set.seed(2026)

# Notion 建议：允许通过环境变量控制置换次数，方便 CI/CD 或快速调试
N_PERM <- as.integer(Sys.getenv("N_PERM", "10000"))

# ----------------------------------------------------------------------------
# 1. Load Data & Calculate Weighted Signature Score
# ----------------------------------------------------------------------------
cat("[1/4] Loading data and calculating platform-limited signature score...\n")

weights_file <- file.path(derived_dir, "locked_weights.csv")
expr_65204_path <- file.path(processed_dir, "expr_gse65204&GPL14550__FULL.rds")
pheno_65204_path <- file.path(processed_dir, "phase3_outputs", "GSE65204_pheno_raw.rds")

# Notion 建议 1：修复 all(file.exists) 的语法 Bug
if (!all(file.exists(c(weights_file, expr_65204_path, pheno_65204_path)))) {
  stop("Required input files missing! Please check data/processed and data/derived.")
}

weights_df <- read_csv(weights_file, show_col_types = FALSE)
expr_mat <- readRDS(expr_65204_path)
pheno_df <- readRDS(pheno_65204_path)

# 从raw_characteristics_all列中解析IgE值
pheno_df <- pheno_df %>%
  mutate(lnige = as.numeric(str_extract(raw_characteristics_all, "lnige: ([0-9.]+)") %>% str_replace("lnige: ", "")))

# Notion 建议 2：强制可用基因数断言与日志输出
avail_genes <- intersect(weights_df$Gene, rownames(expr_mat))
missing_genes <- setdiff(weights_df$Gene, avail_genes)
cat(sprintf("  GSE65204 available genes: %d/11\n", length(avail_genes)))
cat(sprintf("  Missing genes: %s\n", ifelse(length(missing_genes)==0, "None", paste(missing_genes, collapse = ", "))))

stopifnot("Error: Platform coverage too low (< 4 genes)." = length(avail_genes) >= 4)
stopifnot("Error: Gene mismatch in matrix." = all(avail_genes %in% rownames(expr_mat)))

# Notion 建议 3：安全、带维度的矩阵乘法 (抛弃 crossprod)
w_matched <- weights_df$Weight[match(avail_genes, weights_df$Gene)]
sig_scores <- as.vector(t(expr_mat[avail_genes, , drop=FALSE]) %*% w_matched)
names(sig_scores) <- colnames(expr_mat)

# ----------------------------------------------------------------------------
# 2. Process IgE & Define Tertiles
# ----------------------------------------------------------------------------
cat("[2/4] Processing IgE levels and defining tertiles...\n")

ige_col <- "lnige"
sample_id_col <- "std_sample_id"

if (is.na(ige_col)) stop("IgE column not found in GSE65204 pheno.")

valid_pheno <- pheno_df %>% filter(!is.na(.data[[ige_col]]))
common <- intersect(names(sig_scores), valid_pheno[[sample_id_col]])

score_sub <- sig_scores[common]
valid_pheno <- valid_pheno[match(common, valid_pheno[[sample_id_col]]), ]
ige_vals <- as.numeric(valid_pheno[[ige_col]])

# lnige already in log form, no need for transformation
cat("  Using pre-transformed lnIgE values.\n")

# Notion 建议 5：使用 dplyr::ntile 防御 ties 问题
df_main <- tibble(
  Sample_ID = common,
  Signature_Score = score_sub,
  lnIgE = ige_vals
) %>%
  mutate(IgE_tertile = factor(dplyr::ntile(lnIgE, 3), labels = c("Low", "Mid", "High")))

cat("  Tertile distribution:\n")
print(table(df_main$IgE_tertile))

# ----------------------------------------------------------------------------
# 3. Generating Panels A & B (Dose-response & Calibration)
# ----------------------------------------------------------------------------
cat("[3/4] Generating Figure S5A (Dose-response) & S5B (Calibration)...\n")

# --- S5A: Dose-response ---
p_s5a <- ggplot(df_main, aes(x = IgE_tertile, y = Signature_Score, fill = IgE_tertile)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 1.8, alpha = 0.6, color = "black") +
  scale_fill_manual(values = c("Low" = "#4DBBD5", "Mid" = "#E64B35", "High" = "#C20008")) +
  labs(title = "A. Dose-response relationship", x = "Serum IgE Tertile", y = "Platform-limited Signature Score") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none", plot.title = element_text(face="bold"))

ggsave(file.path(supp_dir, "FigS5A_Dose_Response.pdf"), p_s5a, width = 5, height = 5)

# --- S5B: Calibration Curve ---
df_cal <- df_main %>% 
  filter(IgE_tertile %in% c("Low", "High")) %>% 
  mutate(y = ifelse(IgE_tertile == "High", 1, 0))

dd <- datadist(df_cal); options(datadist='dd')
fit <- lrm(y ~ Signature_Score, data = df_cal, x=TRUE, y=TRUE)
pred_prob <- predict(fit, type="fitted")

pdf(file.path(supp_dir, "FigS5B_Calibration_Curve.pdf"), width = 5, height = 5, family = "Helvetica", pointsize = 14)
par(mar = c(5, 5, 4, 2), cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.1)
cal_stats <- val.prob(pred_prob, df_cal$y, statloc = FALSE, legendloc = "bottomright", 
                      xlab = "Predicted Probability (T2-High)", ylab = "Actual Proportion (T2-High)")
dev.off()

# Notion 建议 6：导出 Calibration 统计指标供审稿人查阅
write_csv(as.data.frame(t(cal_stats)), file.path(supp_dir, "FigS5B_calibration_stats.csv"))
cat("  Saved Calibration stats (Brier, slope, intercept) to FigS5B_calibration_stats.csv\n")

# ----------------------------------------------------------------------------
# 4. Generating Panel C (Permutations)
# ----------------------------------------------------------------------------
cat(sprintf("[4/4] Running %d Permutations for S5C...\n", N_PERM))

obs_rho <- cor(df_main$Signature_Score, df_main$lnIgE, method = "spearman")
obs_delta <- median(df_main$Signature_Score[df_main$IgE_tertile == "High"]) - 
  median(df_main$Signature_Score[df_main$IgE_tertile == "Low"])

scores <- df_main$Signature_Score
iges <- df_main$lnIgE
terts <- df_main$IgE_tertile

null_rho <- numeric(N_PERM)
null_delta <- numeric(N_PERM)

# Notion 建议 7：完全独立置换，保证 Null Distribution 逻辑纯粹性
for (i in seq_len(N_PERM)) {
  null_rho[i] <- cor(scores, sample(iges), method = "spearman")
  
  shuf_tert <- sample(terts)
  null_delta[i] <- median(scores[shuf_tert == "High"]) - median(scores[shuf_tert == "Low"])
}

p_rho <- (sum(abs(null_rho) >= abs(obs_rho)) + 1) / (N_PERM + 1)
p_delta <- (sum(abs(null_delta) >= abs(obs_delta)) + 1) / (N_PERM + 1)

# Plotting S5C
p_rho_plot <- ggplot(tibble(val = null_rho), aes(x = val)) +
  geom_histogram(bins = 60, fill = "grey80", color = "grey40") +
  geom_vline(xintercept = obs_rho, color = "#E64B35", linewidth = 1.2, linetype = "dashed") +
  annotate("text", x = 0, y = Inf, label = sprintf("Observed \u03C1 = %.3f\nP = %.1e", obs_rho, p_rho), 
           vjust = 1.5, hjust = 0.5, color = "black", fontface = "bold") +
  labs(title = "C. Permutation test validation", x = "Permuted Spearman \u03C1", y = "Frequency") +
  ylim(0, 700) +
  theme_classic(base_size = 12) + theme(plot.title = element_text(face="bold"))

p_del_plot <- ggplot(tibble(val = null_delta), aes(x = val)) +
  geom_histogram(bins = 60, fill = "grey80", color = "grey40") +
  geom_vline(xintercept = obs_delta, color = "#E64B35", linewidth = 1.2, linetype = "dashed") +
  annotate("text", x = 0, y = Inf, label = sprintf("Observed \u0394 = %.3f\nP = %.1e", obs_delta, p_delta), 
           vjust = 1.5, hjust = 0.5, color = "black", fontface = "bold") +
  labs(title = "", x = "Permuted Tertile Gradient (\u0394 High-Low)", y = "Frequency") +
  ylim(0, 700) +
  theme_classic(base_size = 12)

p_s5c <- plot_grid(p_rho_plot, p_del_plot, ncol = 2, align = "h")

# Notion 建议 8：所有输出转移到 output/supplement
ggsave(file.path(supp_dir, "FigS5C_Permutation_Test.pdf"), p_s5c, width = 10, height = 4.5)

# ----------------------------------------------------------------------------
# 5. Save Source Data & Log
# ----------------------------------------------------------------------------
write_csv(df_main, file.path(supp_dir, "FigS5_Source_Data.csv"))

if(exists("log_output")) {
  log_output(file.path(supp_dir, "FigS5A_Dose_Response.pdf"), "figure", "Fig S5A", "Dose-response boxplots", "19_supp_figS5_robustness.R")
  log_output(file.path(supp_dir, "FigS5B_Calibration_Curve.pdf"), "figure", "Fig S5B", "Logistic calibration curve", "19_supp_figS5_robustness.R")
  log_output(file.path(supp_dir, "FigS5C_Permutation_Test.pdf"), "figure", "Fig S5C", "10,000 permutations null distributions", "19_supp_figS5_robustness.R")
}

cat("\n=== Figure S5 Generation Completed ===\n")