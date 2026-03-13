# =========================================================================
# Script: 12_supp_tables_S5_to_S7.R
# Purpose: Generate Supplementary Tables S5 to S7 (Real Data Computation)
# Author: Zuo + Notion QA hardening
# Date: 2026-03-01 (finalized)
#
# Inputs:
#   - data/derived/locked_weights.csv
#   - data/processed/GSE65204_expr_z.rds
#   - data/processed/GSE65204_pheno.csv
#   - data/raw/drug_network/{categories.tsv,drugs.tsv,genes.tsv,interactions.tsv}  (preferred for S7)
#   - output/figures_main/**/gsea*.csv OR other gsea*.csv (for S6; auto-scan)
#   - R/03_index_logger.R
#
# Outputs:
#   - output/supplement/TableS5_NRI_IDI.csv
#   - output/supplement/TableS6_GSEA_full.csv
#   - output/supplement/TableS7_drug_interaction_full.csv
#   - output/logs/12_supp_tables_S5_to_S7.log
# =========================================================================

# ----------------------------
# 0) Project root and paths
# ----------------------------
base_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# Validate project root by checking for key directories
if (!dir.exists(file.path(base_dir, "data")) ||
    !dir.exists(file.path(base_dir, "analysis")) ||
    !dir.exists(file.path(base_dir, "data_preparation"))) {
  stop("Please run this script from the project root directory (the folder containing data/, analysis/, data_preparation/).",
       call. = FALSE)
}

processed_dir   <- file.path(base_dir, "data", "processed")
derived_dir     <- file.path(base_dir, "data", "derived")
raw_dir         <- file.path(base_dir, "data", "raw")
supplement_dir  <- file.path(base_dir, "output", "supplement")
logs_dir        <- file.path(base_dir, "output", "logs")

dir.create(supplement_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# 1) Logging (safe sink)
# ----------------------------
log_file <- file.path(logs_dir, "12_supp_tables_S5_to_S7.log")
sink(log_file, append = FALSE, split = TRUE)
on.exit({
  sink()
}, add = TRUE)

cat(sprintf("Base directory: %s\n", base_dir))
cat("Outputs:\n")
cat(sprintf("  - %s\n", file.path(supplement_dir, "TableS5_NRI_IDI.csv")))
cat(sprintf("  - %s\n", file.path(supplement_dir, "TableS6_drug_pathway_full.csv")))
cat(sprintf("  - %s\n", file.path(supplement_dir, "TableS6_GSEA_enrichment_full.csv")))
cat(sprintf("  - %s\n", file.path(supplement_dir, "TableS7_drug_interaction_full.csv")))
cat(sprintf("  - %s\n\n", log_file))

# ----------------------------
# 2) Packages (no auto-install)
# ----------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(stringr)
})

# ----------------------------
# 3) Helper function: Format Table S6 GSEA
# ----------------------------
format_table_s6_gsea <- function(df_raw) {
  # 1) 统一列名为小写+下划线，便于兼容不同 fgsea 输出
  nm0 <- names(df_raw)
  names(df_raw) <- nm0 |> 
    stringr::str_trim() |> 
    stringr::str_to_lower() |> 
    stringr::str_replace_all("[^a-z0-9]+", "_")

  # 2) 列名映射（兼容常见写法）
  if (!"pathway" %in% names(df_raw) && "term" %in% names(df_raw)) df_raw$pathway <- df_raw$term
  if (!"pathway" %in% names(df_raw) && "pathway_name" %in% names(df_raw)) df_raw$pathway <- df_raw$pathway_name
  if (!"pval" %in% names(df_raw) && "p_value" %in% names(df_raw)) df_raw$pval <- df_raw$p_value
  if (!"pval" %in% names(df_raw) && "pvalue" %in% names(df_raw)) df_raw$pval <- df_raw$pvalue
  if (!"padj" %in% names(df_raw) && "fdr" %in% names(df_raw)) df_raw$padj <- df_raw$fdr
  if (!"padj" %in% names(df_raw) && "qvalue" %in% names(df_raw)) df_raw$padj <- df_raw$qvalue
  if (!"padj" %in% names(df_raw) && "p_adjust" %in% names(df_raw)) df_raw$padj <- df_raw$p_adjust
  if (!"leadingedge" %in% names(df_raw) && "leading_edge" %in% names(df_raw)) df_raw$leadingedge <- df_raw$leading_edge
  if (!"leadingedge" %in% names(df_raw) && "leadingedge" %in% names(df_raw)) df_raw$leadingedge <- df_raw$leadingedge
  if (!"leadingedge" %in% names(df_raw) && "leadingEdge" %in% names(df_raw)) df_raw$leadingedge <- df_raw$leadingEdge
  if (!"nes" %in% names(df_raw) && "NES" %in% names(df_raw)) df_raw$nes <- df_raw$NES

  # 3) 强校验
  need <- c("pathway", "size", "nes", "pval", "padj", "leadingedge")
  miss <- setdiff(need, names(df_raw))
  if (length(miss) > 0) {
    stop("Table S6 input missing required columns: ", paste(miss, collapse = ", "))
  }

  # 4) Leading edge 截断到前5个
  # leadingedge 可能是 list-column 或字符串，这里统一成字符向量
  leading_edge_str <- vapply(df_raw$leadingedge, function(x) {
    genes <- NULL
    if (is.list(x)) {
      genes <- unlist(x)
    } else {
      # 常见格式： "GENE1,GENE2,GENE3" 或 "c(\"A\",\"B\")"
      s <- as.character(x)
      s <- stringr::str_replace_all(s, "[\\[\\]\\(\\)\\{\\}\\\"\\']", "")
      s <- stringr::str_replace_all(s, "\\s+", "")
      genes <- unlist(strsplit(s, ",", fixed = TRUE))
      genes <- genes[nzchar(genes)]
    }
    if (length(genes) == 0 || all(is.na(genes)) || all(genes == "")) {
      "Not available"
    } else if (length(genes) <= 5) {
      paste(genes, collapse = ", ")
    } else {
      paste0(paste(genes[1:5], collapse = ", "), ", ...")
    }
  }, character(1))

  # 5) Pathway 清洗：去掉 KEGG_ 前缀，下划线→空格，首字母大写（Title Case）
  pathway_clean <- df_raw$pathway |> 
    as.character() |> 
    stringr::str_replace("^kegg_", "") |> 
    stringr::str_replace("^KEGG_", "") |> 
    stringr::str_replace_all("_", " ") |> 
    stringr::str_to_lower() |> 
    stringr::str_to_title()

  out <- df_raw %>% 
    dplyr::mutate(
      `Collection` = "KEGG",
      `Pathway` = pathway_clean,
      `Size` = as.integer(.data$size),
      `NES` = as.numeric(.data$nes),
      `P -value` = as.numeric(.data$pval),
      `FDR` = as.numeric(.data$padj),
      `Pass network filter` = ifelse(.data$padj <= 0.10, "YES", "NO"),
      `Leading edge genes (truncated)` = leading_edge_str
    ) %>% 
    dplyr::select(
      `Collection`,
      `Pathway`,
      `Size`,
      `NES`,
      `P -value`,
      `FDR`,
      `Pass network filter`,
      `Leading edge genes (truncated)`
    )

  out
}

# ----------------------------
# 4) Helper function: Format Table S7 DGIdb
# ----------------------------
format_table_s7_dgidb <- function(df_raw) {
  nm0 <- names(df_raw)
  names(df_raw) <- nm0 |> 
    stringr::str_trim() |> 
    stringr::str_to_lower() |> 
    stringr::str_replace_all("[^a-z0-9]+", "_")

  # 识别列
  gene_col <- grep("gene|target", names(df_raw), value = TRUE)[1]
  drug_col <- grep("drug|compound", names(df_raw), value = TRUE)[1]
  type_col <- grep("interaction|type", names(df_raw), value = TRUE)[1]
  db_col   <- grep("database|source_db|source_databases|source", names(df_raw), value = TRUE)[1]

  if (is.na(gene_col) || is.na(drug_col)) {
    stop("Table S7 input missing gene/drug columns.")
  }

  out <- df_raw %>% 
    dplyr::mutate(
      `Target Gene` = toupper(as.character(.data[[gene_col]])),
      `Drug Name` = as.character(.data[[drug_col]]),
      `Interaction Type` = dplyr::if_else(
        !is.na(type_col) & !is.na(.data[[type_col]]) & nzchar(as.character(.data[[type_col]])),
        as.character(.data[[type_col]]),
        "Unknown/Unspecified"
      ),
      `Source Status` = "Unknown",
      `Source Evidence` = "none",
      `Sources Count` = 1L,
      `Source Databases` = dplyr::if_else(
        !is.na(db_col) & !is.na(.data[[db_col]]) & nzchar(as.character(.data[[db_col]])),
        as.character(.data[[db_col]]),
        "NCI"
      ),
      `Context Sources` = "NCI"
    ) %>% 
    dplyr::select(
      `Target Gene`,
      `Drug Name`,
      `Interaction Type`,
      `Source Status`,
      `Source Evidence`,
      `Sources Count`,
      `Source Databases`,
      `Context Sources`
    ) %>% 
    dplyr::distinct()

  out
}

# ----------------------------
# 3) Source utilities
# ----------------------------
idx_logger_path <- file.path(base_dir, "R", "03_index_logger.R")
if (file.exists(idx_logger_path)) {
  source(idx_logger_path)
} else {
  cat("WARNING: R/03_index_logger.R not found. OUTPUTS_INDEX logging will be skipped.\n")
}

# ----------------------------
# 4) Load locked weights
# ----------------------------
weights_file <- file.path(derived_dir, "locked_weights.csv")
if (!file.exists(weights_file)) stop("locked_weights.csv not found!", call. = FALSE)

locked_weights <- read_csv(weights_file, show_col_types = FALSE)
stopifnot(all(c("Gene", "Weight") %in% colnames(locked_weights)))

sig_genes <- locked_weights$Gene
cat(sprintf("Locked genes loaded: %d\n", length(sig_genes)))
cat(sprintf("Locked genes: %s\n\n", paste(sig_genes, collapse = ", ")))

# ----------------------------
# 5) Input resolution function
# ----------------------------
resolve_input <- function(base_dir, prefix, types) {
  # Define search directories in fixed priority order
  search_dirs <- c(
    file.path(base_dir, "data", "processed"),
    file.path(base_dir, "data", "processed_diet"),
    file.path(base_dir, "data", "processed_full"),
    file.path(base_dir, "data", "processed", "phase3_outputs")
  )
  
  results <- list()
  
  for (type in types) {
    type_results <- list()
    
    for (dir in search_dirs) {
      if (!dir.exists(dir)) next
      
      # Define pattern based on type
      if (type == "expr") {
        pattern <- "^GSE65204_.*expr.*\\.(rds|RDS)$"
      } else if (type == "pheno") {
        pattern <- "^GSE65204_.*pheno.*\\.(csv|CSV|rds|RDS)$"
      } else {
        next
      }
      
      # List files matching pattern
      files <- list.files(path = dir, pattern = pattern, full.names = TRUE)
      if (length(files) > 0) {
        # Sort files to ensure consistent selection
        files <- sort(files)
        
        # For expr, prioritize files containing "expr_z"
        if (type == "expr") {
          expr_z_files <- files[grepl("expr_z", files, ignore.case = TRUE)]
          if (length(expr_z_files) > 0) {
            files <- expr_z_files
          }
        }
        
        # For pheno, prioritize files named "pheno.csv"
        if (type == "pheno") {
          pheno_csv_files <- files[grepl("pheno\\.csv$", files, ignore.case = TRUE)]
          if (length(pheno_csv_files) > 0) {
            files <- pheno_csv_files
          } else {
            # Then try pheno_raw or pheno_std
            pheno_raw_std_files <- files[grepl("pheno_(raw|std)", files, ignore.case = TRUE)]
            if (length(pheno_raw_std_files) > 0) {
              files <- pheno_raw_std_files
            }
          }
        }
        
        type_results[[dir]] <- files
      }
    }
    
    # Flatten results and take the first one (highest priority)
    all_files <- unlist(type_results)
    if (length(all_files) > 0) {
      results[[type]] <- all_files[1]
    } else {
      results[[type]] <- NA
    }
  }
  
  return(results)
}

# =====================================================================
# Table S5: Incremental value metrics (continuous NRI and IDI)
# - LOOCV out-of-fold probabilities (logistic regression)
# - 2000 bootstrap replicates for 95% CI + empirical P
# - Full model is platform-limited (available genes in GSE65204 processed matrix)
# - Overlap module fixed and must remain intact (4 genes):
#     SERPINB2, MUC5AC, CPA3, HDC
# =====================================================================
cat("\n>>> Generating Table S5: NRI/IDI metrics (LOOCV + bootstrap)...\n")

table_s5_output   <- file.path(supplement_dir, "TableS5_NRI_IDI.csv")

# Resolve input files using the resolve_input function
input_files <- resolve_input(base_dir, "GSE65204", c("expr", "pheno"))
expr_65204_path <- input_files$expr
pheno_65204_path <- input_files$pheno

if (is.na(expr_65204_path) || is.na(pheno_65204_path)) {
  table_s5 <- tibble(
    Comparison = "Platform-limited Full vs Overlap (Core Epithelial–Mast Cell Module)",
    Outcome = NA_character_,
    N = NA_integer_,
    N_event = NA_integer_,
    N_nonevent = NA_integer_,
    Full_coverage = NA_character_,
    Overlap_coverage = NA_character_,
    Missing_genes_full = NA_character_,
    Missing_genes_overlap = NA_character_,
    NRI = NA_real_,
    NRI_CI_low = NA_real_,
    NRI_CI_high = NA_real_,
    NRI_P_emp = NA_real_,
    IDI = NA_real_,
    IDI_CI_low = NA_real_,
    IDI_CI_high = NA_real_,
    IDI_P_emp = NA_real_,
    Notes = "GSE65204 expr/pheno files missing."
  )
  write_csv(table_s5, table_s5_output)
  cat("  ⚠ GSE65204 data missing. Generated placeholder Table S5.\n")
  
} else {
  
  # Log the resolved files
  cat(sprintf("  Resolved expr file: %s\n", expr_65204_path))
  cat(sprintf("  Resolved pheno file: %s\n", pheno_65204_path))
  
  expr_data <- readRDS(expr_65204_path)
  
  # Read pheno file based on extension
  pheno_ext <- tolower(tools::file_ext(pheno_65204_path))
  if (pheno_ext %in% c("csv", "tsv")) {
    pheno_data <- read_csv(pheno_65204_path, show_col_types = FALSE)
  } else if (pheno_ext %in% c("rds")) {
    pheno_data <- readRDS(pheno_65204_path)
  } else {
    stop(sprintf("Unsupported pheno file format: %s", pheno_ext), call. = FALSE)
  }
  
  # --------
  # Harmonize expr matrix orientation
  # Expect: gene x sample
  # --------
  if (!is.null(dim(expr_data))) {
    cat(sprintf("Expression object dims: %d x %d\n", nrow(expr_data), ncol(expr_data)))
  } else {
    stop("Expression object is not a matrix-like object.", call. = FALSE)
  }
  
  # Heuristic: if rownames look like GSM and colnames look like genes, transpose
  # (Conservative: only transpose when many rownames start with "GSM" and many colnames are in sig_genes)
  looks_like_gsm <- function(x) mean(grepl("^GSM", x)) > 0.5
  looks_like_gene <- function(x) mean(x %in% sig_genes) > 0.3
  
  if (!is.null(rownames(expr_data)) && !is.null(colnames(expr_data))) {
    if (looks_like_gsm(rownames(expr_data)) && looks_like_gene(colnames(expr_data))) {
      cat("Detected expr matrix as sample x gene. Transposing to gene x sample.\n")
      expr_data <- t(expr_data)
    }
  }
  
  if (is.null(rownames(expr_data)) || is.null(colnames(expr_data))) {
    stop("Expression matrix must have rownames (genes) and colnames (samples).", call. = FALSE)
  }
  
  # --------
  # Outcome selection (binary) - prioritized by clinical anchor
  # --------
  outcome_col <- dplyr::case_when(
    "IgE_High_binary" %in% names(pheno_data) ~ "IgE_High_binary",
    "Case_Binary" %in% names(pheno_data) ~ "Case_Binary",
    "T2_Binary" %in% names(pheno_data) ~ "T2_Binary",
    "std_group_label" %in% names(pheno_data) ~ "std_group_label",
    "asthma_status" %in% names(pheno_data) ~ "asthma_status",
    "lnige" %in% names(pheno_data) ~ "lnige",
    TRUE ~ NA_character_
  )
  
  # Create binary outcome if needed
  if (outcome_col == "std_group_label") {
    # Create binary outcome from std_group_label
    pheno_data$Case_Binary <- ifelse(pheno_data$std_group_label == "case", 1, 
                                     ifelse(pheno_data$std_group_label == "control", 0, NA))
    outcome_col <- "Case_Binary"
  } else if (outcome_col == "asthma_status") {
    # Create binary outcome from asthma_status
    pheno_data$Case_Binary <- ifelse(pheno_data$asthma_status == "Asthma", 1, 
                                     ifelse(pheno_data$asthma_status == "Control", 0, NA))
    outcome_col <- "Case_Binary"
  } else if (outcome_col == "lnige") {
    # Create IgE_High_binary using median split
    median_ige <- median(as.numeric(pheno_data$lnige), na.rm = TRUE)
    pheno_data$IgE_High_binary <- ifelse(as.numeric(pheno_data$lnige) > median_ige, 1, 0)
    outcome_col <- "IgE_High_binary"
  } else if (is.na(outcome_col)) {
    stop("Binary outcome missing: IgE_High_binary, T2_Binary, Case_Binary, std_group_label, asthma_status, or lnige not found.", call. = FALSE)
  }
  
  # Check for IgE tertile fields and apply subset if available
  subset_rule <- "not_available"
  N_before <- nrow(pheno_data)
  N_after <- N_before
  
  if (any(c("IgE_tertile", "IgE_group") %in% names(pheno_data))) {
    # Use the first available tertile field
    tertile_col <- ifelse("IgE_tertile" %in% names(pheno_data), "IgE_tertile", "IgE_group")
    
    # Subset to high and low tertiles, exclude middle
    pheno_data <- pheno_data %>% 
      filter(.data[[tertile_col]] %in% c("high", "low", "High", "Low", 1, 3))
    
    N_after <- nrow(pheno_data)
    subset_rule <- "IgE tertile high vs low (middle excluded)"
    
    cat(sprintf("Applied IgE tertile subset: %s | N before=%d, N after=%d\n", subset_rule, N_before, N_after))
  }
  
  valid_pheno <- pheno_data %>% dplyr::filter(!is.na(.data[[outcome_col]]))
  
  # Try different sample ID column names
  sample_id_col <- dplyr::case_when(
    "Sample_ID" %in% names(valid_pheno) ~ "Sample_ID",
    "std_sample_id" %in% names(valid_pheno) ~ "std_sample_id",
    "sample_id" %in% names(valid_pheno) ~ "sample_id",
    "geo_accession" %in% names(valid_pheno) ~ "geo_accession",
    TRUE ~ NA_character_
  )
  
  if (is.na(sample_id_col)) {
    stop("pheno file must contain column 'Sample_ID', 'std_sample_id', 'sample_id', or 'geo_accession' matching expr colnames.", call. = FALSE)
  }
  
  common <- intersect(colnames(expr_data), valid_pheno[[sample_id_col]])
  if (length(common) < 10) stop("Too few matched samples between expr and pheno.", call. = FALSE)
  
  expr_subset <- expr_data[, common, drop = FALSE]
  valid_pheno <- valid_pheno[match(common, valid_pheno[[sample_id_col]]), , drop = FALSE]
  y <- as.integer(valid_pheno[[outcome_col]])
  if (!all(y %in% c(0L, 1L))) stop("Outcome must be binary 0/1.", call. = FALSE)
  
  n <- length(y)
  n_event <- sum(y == 1L)
  n_nonevent <- sum(y == 0L)
  cat(sprintf("Outcome: %s | N=%d (event=%d, nonevent=%d)\n", outcome_col, n, n_event, n_nonevent))
  
  # --------
  # Module definition (LOCKED)
  # Overlap module = Core Epithelial–Mast Cell Module (baseline)
  # --------
  overlap_genes <- c("SERPINB2", "MUC5AC", "CPA3", "HDC")
  
  # Full model = platform-limited implementation: available locked genes on this platform
  full_avail <- intersect(sig_genes, rownames(expr_subset))
  
  base_avail <- intersect(overlap_genes, rownames(expr_subset))
  missing_full <- setdiff(sig_genes, full_avail)
  missing_base <- setdiff(overlap_genes, base_avail)
  
  cat(sprintf("Full available genes: %d/11\n", length(full_avail)))
  cat(sprintf("Full missing genes: %s\n", ifelse(length(missing_full) == 0, "NA", paste(missing_full, collapse = ", "))))
  cat(sprintf("Overlap available genes: %d/4\n", length(base_avail)))
  cat(sprintf("Overlap missing genes: %s\n", ifelse(length(missing_base) == 0, "NA", paste(missing_base, collapse = ", "))))
  
  if (length(base_avail) != 4) {
    stop("Overlap module is expected to remain intact (4/4 available) in GSE65204. Please check expr gene symbols / mapping.", call. = FALSE)
  }
  if (length(full_avail) < 2) {
    stop("Too few available genes for platform-limited full model.", call. = FALSE)
  }
  
  # --------
  # Locked weighted score (weighted sum on available genes)
  # --------
  w_tbl <- locked_weights %>% dplyr::select(Gene, Weight)
  
  calc_locked_score <- function(expr_mat, weights_tbl, genes) {
    avail <- intersect(genes, rownames(expr_mat))
    if (length(avail) == 0) return(rep(NA_real_, ncol(expr_mat)))
    w2 <- weights_tbl %>% dplyr::filter(Gene %in% avail)
    w2 <- w2[match(avail, w2$Gene), , drop = FALSE]
    as.numeric(crossprod(w2$Weight, expr_mat[avail, , drop = FALSE]))
  }
  
  score_base <- calc_locked_score(expr_subset, w_tbl, overlap_genes)
  score_new  <- calc_locked_score(expr_subset, w_tbl, full_avail)
  
  # --------
  # LOOCV predicted probabilities (1D logistic regression)
  # --------
  loocv_prob <- function(score, y) {
    stopifnot(length(score) == length(y))
    n <- length(y)
    probs <- rep(NA_real_, n)

    for (i in seq_len(n)) {
      train <- data.frame(
        y = y[-i],
        score = score[-i]
      )
      test <- data.frame(score = score[i])

      if (anyNA(train$y) || anyNA(train$score) || anyNA(test$score)) {
        stop(sprintf("LOOCV fail: NA detected at i=%d", i))
      }

      m <- tryCatch(glm(y ~ score, family = binomial(), data = train),
                    error = function(e) stop(sprintf("LOOCV glm fail i=%d: %s", i, e$message)))

      p <- tryCatch(as.numeric(predict(m, newdata = test, type = "response")),
                    warning = function(w) stop(sprintf("LOOCV predict warning i=%d: %s", i, w$message)),
                    error = function(e) stop(sprintf("LOOCV predict error i=%d: %s", i, e$message)))

      if (length(p) != 1 || is.na(p)) stop(sprintf("LOOCV fail: invalid prob at i=%d", i))
      probs[i] <- p
    }

    stopifnot(!anyNA(probs))
    probs
  }
  
  p_base <- loocv_prob(score_base, y)
  p_new  <- loocv_prob(score_new,  y)
  
  # --------
  # Category-free continuous NRI + IDI
  # ties count as 0 (neither up nor down)
  # --------
  calc_nri_idi <- function(p_new, p_base, y) {
    ev <- y == 1L
    ne <- y == 0L
    
    idi <- (mean(p_new[ev]) - mean(p_base[ev])) - (mean(p_new[ne]) - mean(p_base[ne]))
    
    up_ev   <- mean(p_new[ev] > p_base[ev])
    down_ev <- mean(p_new[ev] < p_base[ev])
    nri_ev  <- up_ev - down_ev
    
    down_ne <- mean(p_new[ne] < p_base[ne])
    up_ne   <- mean(p_new[ne] > p_base[ne])
    nri_ne  <- down_ne - up_ne
    
    nri <- nri_ev + nri_ne
    c(NRI = nri, IDI = idi)
  }
  
  # Gate checks to ensure no NA values
  stopifnot(!anyNA(p_base), !anyNA(p_new))
  pt <- calc_nri_idi(p_new, p_base, y)
  stopifnot(!anyNA(pt))
  
  # --------
  # Bootstrap 2000: 95% CI + empirical two-sided P (vs 0)
  # --------
  set.seed(2026)
  B <- 2000
  boot_mat <- matrix(NA_real_, nrow = B, ncol = 2)
  colnames(boot_mat) <- c("NRI", "IDI")
  
  for (b in seq_len(B)) {
    idx <- sample.int(n, size = n, replace = TRUE)
    boot_mat[b, ] <- calc_nri_idi(p_new[idx], p_base[idx], y[idx])
  }
  
  ci <- apply(boot_mat, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
  p_emp <- sapply(seq_len(2), function(j) {
    mean(abs(boot_mat[, j]) >= abs(pt[j]), na.rm = TRUE)
  })
  
  # Generate audit notes
  audit_notes <- paste0(
    "LOOCV out-of-fold probabilities from logistic regression; 2000 bootstrap replicates for CI and empirical P.",
    " expr_file_used=", expr_65204_path,
    "; pheno_file_used=", pheno_65204_path,
    "; expr_dim=", paste(dim(expr_data), collapse = "x"),
    "; n_common_samples=", length(common),
    "; subset_rule=", subset_rule,
    "; N_before=", N_before,
    "; N_after=", N_after
  )
  
  table_s5 <- tibble(
    Comparison = "Platform-limited Full vs Overlap (Core Epithelial–Mast Cell Module)",
    Outcome = outcome_col,
    N = n,
    N_event = n_event,
    N_nonevent = n_nonevent,
    Full_coverage = paste0(length(full_avail), "/11"),
    Overlap_coverage = "4/4",
    Missing_genes_full = ifelse(length(missing_full) == 0, "NA", paste(missing_full, collapse = ";")),
    Missing_genes_overlap = "NA",
    NRI = unname(pt["NRI"]),
    NRI_CI_low = ci[1, "NRI"],
    NRI_CI_high = ci[2, "NRI"],
    NRI_P_emp = p_emp[1],
    IDI = unname(pt["IDI"]),
    IDI_CI_low = ci[1, "IDI"],
    IDI_CI_high = ci[2, "IDI"],
    IDI_P_emp = p_emp[2],
    Notes = audit_notes
  )
  
  write_csv(table_s5, table_s5_output)
  cat(paste0("Saved Table S5: ", basename(table_s5_output), "\n"))
}

# =====================================================================
# Table S6: Drug → Pathway mapping full table and GSEA enrichment
# - Drug-pathway: auto-scan
# - GSEA enrichment: only read Fig1C_GSEA_results.csv
# =====================================================================
cat("\n>>> Generating Table S6: Drug-Pathway mapping and GSEA enrichment tables...\n")

# Define output files
table_s6_drug_pathway_output <- file.path(supplement_dir, "TableS6_drug_pathway_full.csv")
table_s6_gsea_enrichment_output <- file.path(supplement_dir, "TableS6_GSEA_enrichment_full.csv")
table_s6_fig1c_output <- file.path(supplement_dir, "TableS6_GSEA_for_Fig1C.csv")

# ---------------------------------------------------------------------
# Process drug-pathway data (auto-scan)
# ---------------------------------------------------------------------
cat("\n[Step 1] Processing drug-pathway data...\n")

gsea_search_paths <- c(
  file.path(base_dir, "output", "figures_main"),
  file.path(base_dir, "output", "tables_main"),
  file.path(base_dir, "output", "supplement"),
  file.path(base_dir, "output"),
  file.path(base_dir, "analysis"),
  file.path(base_dir, "data")
)

gsea_files <- unlist(lapply(gsea_search_paths, function(p) {
  if (!dir.exists(p)) return(character(0))
  list.files(path = p, pattern = "(drug.*pathway|pathway.*drug).*\\.(csv|tsv)$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
}))

gsea_files <- unique(gsea_files)

# Exclude output files generated by this script
exclude <- c("TableS6_drug_pathway_full.csv",
             "TableS6_GSEA_enrichment_full.csv",
             "TableS6_GSEA_for_Fig1C.csv",
             "TableS6_GSEA_full.csv")
gsea_files <- gsea_files[!basename(gsea_files) %in% exclude]

# Function to extract cohort_id from filename
extract_cohort_id <- function(filename) {
  # Try to match GSE IDs
  gse_match <- stringr::str_extract(filename, "GSE[0-9]+")
  if (!is.na(gse_match)) {
    return(gse_match)
  }
  # Try to match other cohort identifiers
  cohort_patterns <- c("asthma", "control", "case", "healthy")
  for (pattern in cohort_patterns) {
    if (stringr::str_detect(filename, pattern, negate = FALSE)) {
      return(toupper(pattern))
    }
  }
  return("")
}

# Function to standardize drug-pathway data frame
standardize_drug_pathway_df <- function(df, source_file) {
  # Convert to tibble and ensure source_file is character
  df <- as_tibble(df)
  
  # Extract cohort_id from filename
  cohort_id <- extract_cohort_id(source_file)
  
  # Create standardized data frame with the same number of rows as df
  standardized_df <- tibble(
    source_file = rep(as.character(source_file), nrow(df)),
    cohort_id = rep(cohort_id, nrow(df)),
    DrugID = rep(NA_character_, nrow(df)),
    EntrezID = rep(NA_character_, nrow(df)),
    DrugName = rep(NA_character_, nrow(df)),
    PathwayID = rep(NA_character_, nrow(df)),
    PathwayName = rep(NA_character_, nrow(df))
  )
  
  # Map columns to standardized names
  col_map <- list(
    DrugID = c("DrugID", "drug_id", "DrugId", "drugid"),
    EntrezID = c("EntrezID", "entrez_id", "EntrezId", "entrezid", "GeneID", "gene_id"),
    DrugName = c("DrugName", "drug_name", "Drug", "drug"),
    PathwayID = c("PathwayID", "pathway_id", "PathwayId", "pathwayid", "ID", "id"),
    PathwayName = c("PathwayName", "pathway_name", "Pathway", "pathway", "Description", "description", "term")
  )
  
  # Fill in values based on column mapping
  for (target_col in names(col_map)) {
    for (source_col in col_map[[target_col]]) {
      if (source_col %in% names(df)) {
        # Safely convert to character
        standardized_df[[target_col]] <- tryCatch(
          as.character(df[[source_col]]),
          error = function(e) rep(NA_character_, nrow(df))
        )
        break
      }
    }
  }
  
  return(standardized_df)
}

# Process drug-pathway files
drug_pathway_list <- list()
for (f in gsea_files) {
  tryCatch({
    # Determine file type and read accordingly
    ext <- tolower(tools::file_ext(f))
    if (ext == "tsv") {
      df <- read_tsv(f, show_col_types = FALSE)
    } else {
      df <- read_csv(f, show_col_types = FALSE)
    }
    
    # Check if data frame is empty
    if (nrow(df) == 0 || ncol(df) == 0) {
      cat(sprintf("  ⚠ Skipping empty file: %s\n", f))
      next
    }
    
    # Use safer method to create relative path without regex
    source_file <- sub(paste0("^", base_dir), "", f)
    
    # Process as drug-pathway file
    standardized_df <- standardize_drug_pathway_df(df, source_file)
    drug_pathway_list[[f]] <- standardized_df
  }, error = function(e) {
    cat(sprintf("  ⚠ Failed to read: %s | %s\n", f, e$message))
  })
}

# Write drug-pathway table
if (length(drug_pathway_list) > 0) {
  table_s6_drug_pathway <- bind_rows(drug_pathway_list)
  
  # Deduplicate
  table_s6_drug_pathway <- table_s6_drug_pathway %>%
    distinct(cohort_id, DrugID, PathwayID, .keep_all = TRUE)
  
  # Calculate QC statistics
  rows_written <- nrow(table_s6_drug_pathway)
  unique_drug_ids <- length(unique(table_s6_drug_pathway$DrugID[!is.na(table_s6_drug_pathway$DrugID)]))
  unique_pathway_ids <- length(unique(table_s6_drug_pathway$PathwayID[!is.na(table_s6_drug_pathway$PathwayID)]))
  
  # Calculate top 3 pathways by count
  top_pathways <- table_s6_drug_pathway %>%
    filter(!is.na(PathwayID)) %>%
    group_by(PathwayID, PathwayName) %>%
    count() %>%
    arrange(desc(n)) %>%
    head(3)
  
  # Write to log
  cat(sprintf("S6 drug-pathway rows_written=%d\n", rows_written))
  cat(sprintf("unique DrugID=%d\n", unique_drug_ids))
  cat(sprintf("unique PathwayID=%d\n", unique_pathway_ids))
  
  if (nrow(top_pathways) > 0) {
    cat("top 3 pathways by count:\n")
    for (i in 1:nrow(top_pathways)) {
      cat(sprintf("  %s (%s): %d\n", top_pathways$PathwayName[i], top_pathways$PathwayID[i], top_pathways$n[i]))
    }
  }
  
  write_csv(table_s6_drug_pathway, table_s6_drug_pathway_output)
  cat(paste0("Saved Table S6 drug-pathway: ", basename(table_s6_drug_pathway_output), "\n"))
} else {
  # Create empty drug-pathway table
  empty_drug_pathway <- tibble(
    source_file = character(),
    cohort_id = character(),
    DrugID = character(),
    EntrezID = character(),
    DrugName = character(),
    PathwayID = character(),
    PathwayName = character()
  )
  write_csv(empty_drug_pathway, table_s6_drug_pathway_output)
  cat("  ⚠ No drug-pathway files found. Generated empty Table S6 drug-pathway table.\n")
  # Write QC stats for empty table
  cat("S6 drug-pathway rows_written=0\n")
  cat("unique DrugID=0\n")
  cat("unique PathwayID=0\n")
}

# ---------------------------------------------------------------------
# Process GSEA enrichment data (only read Fig1C_GSEA_results.csv)
# ---------------------------------------------------------------------
cat("\n[Step 2] Processing GSEA enrichment data...\n")

# Read Fig1C_GSEA_results.csv
fig1c_gsea_file <- file.path(base_dir, "output", "tables_main", "Fig1C_GSEA_results.csv")

if (!file.exists(fig1c_gsea_file)) {
  stop("Fig1C_GSEA_results.csv not found. Please run analysis/02_main_validation/07_fig1_discovery.R first to generate it.")
}

cat(paste0("Reading GSEA results from: ", fig1c_gsea_file, "\n"))

# Read the file
gsea_df <- read_csv(fig1c_gsea_file, show_col_types = FALSE)

if (nrow(gsea_df) == 0) {
  stop("Fig1C_GSEA_results.csv is empty. Please check the Fig1C script output.")
}

# Format Table S6 using the standardized function
table_s6_gsea <- format_table_s6_gsea(gsea_df)

# Verify column names
required_s6_cols <- c("Collection", "Pathway", "Size", "NES", "P -value", "FDR", "Pass network filter", "Leading edge genes (truncated)")
stopifnot(identical(names(table_s6_gsea), required_s6_cols))

# Assertion prints for Table S6
cat("=== Table S6 Assertion Prints ===\n")

# 1) Check column names
cat("1. Column names:\n")
print(names(table_s6_gsea))

# 2) Check FDR is numeric and has reasonable range
cat("\n2. FDR summary:\n")
table_s6_gsea$FDR_num <- as.numeric(table_s6_gsea$FDR)
print(summary(table_s6_gsea$FDR_num))

# 3) Top 10 by FDR (numeric)
cat("\n3. Top 10 by FDR (numeric):\n")
top10_fdr <- table_s6_gsea %>% 
  mutate(FDR_num = as.numeric(FDR)) %>% 
  arrange(FDR_num) %>% 
  select(Pathway, NES, FDR_num) %>% 
  head(10)
print(top10_fdr)

# 4) Re-counts (using numeric)
cat("\n4. Re-counts:\n")
cat("FDR <= 0.10:", sum(as.numeric(table_s6_gsea$FDR) <= 0.10, na.rm=TRUE), "\n")
cat("FDR <= 0.05:", sum(as.numeric(table_s6_gsea$FDR) <= 0.05, na.rm=TRUE), "\n")
cat("NES > 0:", sum(as.numeric(table_s6_gsea$NES) > 0, na.rm=TRUE), "\n")
cat("NES < 0:", sum(as.numeric(table_s6_gsea$NES) < 0, na.rm=TRUE), "\n")

# 5) Additional checks: nrow and NA values
cat("\n5. Additional checks:\n")
cat("Total rows:", nrow(table_s6_gsea), "\n")
cat("NA in NES:", sum(is.na(table_s6_gsea$NES)), "\n")
cat("NA in FDR:", sum(is.na(table_s6_gsea$FDR)), "\n")

cat("================================\n\n")

# Write GSEA enrichment table
write_csv(table_s6_gsea, table_s6_gsea_enrichment_output)
cat(paste0("Saved Table S6 GSEA enrichment: ", basename(table_s6_gsea_enrichment_output), "\n"))

# Generate Figure 1C ready table (using the formatted data)
table_s6_fig1c <- table_s6_gsea %>%
  filter(!is.na(`FDR`) & `FDR` <= 0.10) %>%
  arrange(`FDR`) %>%
  head(35)

if (nrow(table_s6_fig1c) > 0) {
  write_csv(table_s6_fig1c, table_s6_fig1c_output)
  cat(paste0("Saved Table S6 for Fig1C: ", basename(table_s6_fig1c_output), "\n"))
}

# Generate minimal CSV with top pathways
table_s6_top_pathways_output <- file.path(supplement_dir, "TableS6_top_pathways.csv")

# FDR top20
fdr_top20 <- table_s6_gsea %>%
  arrange(`FDR`) %>%
  head(20) %>%
  mutate(type = "FDR_top20")

# NES top20 (max)
nes_top20 <- table_s6_gsea %>%
  arrange(desc(`NES`)) %>%
  head(20) %>%
  mutate(type = "NES_top20")

# NES bottom20 (min)
nes_bottom20 <- table_s6_gsea %>%
  arrange(`NES`) %>%
  head(20) %>%
  mutate(type = "NES_bottom20")

# Combine and write
result <- bind_rows(fdr_top20, nes_top20, nes_bottom20)
write_csv(result, table_s6_top_pathways_output)
cat(paste0("Saved Table S6 top pathways: ", basename(table_s6_top_pathways_output), "\n"))

# =====================================================================
# Table S7: Drug–gene interaction full table
# Preferred input: data/raw/drug_network/*.tsv (DGIdb-derived network)
# Fallback: auto-scan csv/tsv with "dgidb/drug/interact" keywords
# =====================================================================
cat("\n>>> Generating Table S7: Drug interaction full table...\n")
table_s7_output <- file.path(supplement_dir, "TableS7_drug_interaction_full.csv")

drug_net_dir <- file.path(raw_dir, "drug_network")
genes_tsv <- file.path(drug_net_dir, "genes.tsv")
drugs_tsv <- file.path(drug_net_dir, "drugs.tsv")
interactions_tsv <- file.path(drug_net_dir, "interactions.tsv")
categories_tsv <- file.path(drug_net_dir, "categories.tsv")

sig_genes_upper <- toupper(sig_genes)

make_empty_s7 <- function(note) {
  tibble(
    Gene_Symbol = character(),
    Drug_Name = character(),
    Interaction_Type = character(),
    Source = character(),
    PMIDs = character(),
    Notes = note
  )
}

if (file.exists(genes_tsv) && file.exists(drugs_tsv) && file.exists(interactions_tsv)) {
  cat("Using preferred drug_network TSV inputs.\n")
  
  genes_df <- read_tsv(genes_tsv, show_col_types = FALSE, progress = FALSE)
  drugs_df <- read_tsv(drugs_tsv, show_col_types = FALSE, progress = FALSE)
  int_df <- read_tsv(interactions_tsv, show_col_types = FALSE, progress = FALSE)
  
  # Format Table S7 using the standardized function
  tryCatch({
    table_s7 <- format_table_s7_dgidb(int_df)
    
    # Filter to signature genes
    table_s7 <- table_s7 %>%
      filter(`Target Gene` %in% sig_genes_upper)
    
    if (nrow(table_s7) == 0) {
      # Create empty table with the correct structure
      table_s7 <- tibble(
        `Target Gene` = character(),
        `Drug Name` = character(),
        `Interaction Type` = character(),
        `Source Status` = character(),
        `Source Evidence` = character(),
        `Sources Count` = integer(),
        `Source Databases` = character(),
        `Context Sources` = character()
      )
    }
  }, error = function(e) {
    cat(sprintf("  ⚠ Error formatting Table S7: %s\n", e$message))
    # Create empty table with the correct structure
    table_s7 <- tibble(
      `Target Gene` = character(),
      `Drug Name` = character(),
      `Interaction Type` = character(),
      `Source Status` = character(),
      `Source Evidence` = character(),
      `Sources Count` = integer(),
      `Source Databases` = character(),
      `Context Sources` = character()
    )
  })
  
} else {
  cat("Preferred drug_network TSV inputs not found. Falling back to auto-scan.\n")
  
  # Fallback scan
  scan_paths <- c(file.path(base_dir, "data"), file.path(base_dir, "output"))
  dgi_files <- unlist(lapply(scan_paths, function(p) {
    if (!dir.exists(p)) return(character(0))
    list.files(path = p, recursive = TRUE, full.names = TRUE, ignore.case = TRUE,
               pattern = "(dgidb|drug.*interact|interaction).*\\.(csv|tsv)$")
  }))
  dgi_files <- unique(dgi_files)
  
  if (length(dgi_files) == 0) {
    # Create empty table with the correct structure
    table_s7 <- tibble(
      `Target Gene` = character(),
      `Drug Name` = character(),
      `Interaction Type` = character(),
      `Source Status` = character(),
      `Source Evidence` = character(),
      `Sources Count` = integer(),
      `Source Databases` = character(),
      `Context Sources` = character()
    )
    cat("  ⚠ No DGIdb files found. Generated empty Table S7.\n")
  } else {
    cat(sprintf("Found candidate drug interaction files: %d\n", length(dgi_files)))
    
    read_any <- function(f) {
      ext <- tolower(tools::file_ext(f))
      tryCatch({
        if (ext == "tsv") read_tsv(f, show_col_types = FALSE, progress = FALSE)
        else read_csv(f, show_col_types = FALSE, progress = FALSE)
      }, error = function(e) NULL)
    }
    
    dgi_list <- lapply(dgi_files, read_any)
    combined <- bind_rows(dgi_list)
    
    if (nrow(combined) == 0) {
      # Create empty table with the correct structure
      table_s7 <- tibble(
        `Target Gene` = character(),
        `Drug Name` = character(),
        `Interaction Type` = character(),
        `Source Status` = character(),
        `Source Evidence` = character(),
        `Sources Count` = integer(),
        `Source Databases` = character(),
        `Context Sources` = character()
      )
    } else {
      # Format Table S7 using the standardized function
      tryCatch({
        table_s7 <- format_table_s7_dgidb(combined)
        
        # Filter to signature genes
        table_s7 <- table_s7 %>%
          filter(`Target Gene` %in% sig_genes_upper)
        
        if (nrow(table_s7) == 0) {
          # Create empty table with the correct structure
          table_s7 <- tibble(
            `Target Gene` = character(),
            `Drug Name` = character(),
            `Interaction Type` = character(),
            `Source Status` = character(),
            `Source Evidence` = character(),
            `Sources Count` = integer(),
            `Source Databases` = character(),
            `Context Sources` = character()
          )
        }
      }, error = function(e) {
        cat(sprintf("  ⚠ Error formatting Table S7: %s\n", e$message))
        # Create empty table with the correct structure
        table_s7 <- tibble(
          `Target Gene` = character(),
          `Drug Name` = character(),
          `Interaction Type` = character(),
          `Source Status` = character(),
          `Source Evidence` = character(),
          `Sources Count` = integer(),
          `Source Databases` = character(),
          `Context Sources` = character()
        )
      })
    }
  }
}

# Verify column names for Table S7
required_s7_cols <- c("Target Gene", "Drug Name", "Interaction Type", "Source Status", "Source Evidence", "Sources Count", "Source Databases", "Context Sources")
stopifnot(identical(names(table_s7), required_s7_cols))

# Write Table S7
write_csv(table_s7, table_s7_output)
cat(paste0("Saved Table S7: ", basename(table_s7_output), "\n"))

# =====================================================================
# Outputs index logging
# =====================================================================
if (exists("log_output")) {
  tryCatch({
    log_output(file_path = table_s5_output, file_type = "table", figure_id = "Table S5",
               description = "Incremental value metrics (continuous NRI/IDI) in GSE65204 anchor (LOOCV + bootstrap)",
               script_name = "12_supp_tables_S5_to_S7.R")
    log_output(file_path = table_s6_drug_pathway_output, file_type = "table", figure_id = "Table S6 (drug–pathway mapping)",
               description = "Drug-Pathway mapping full table (auto-scan, standardized schema)",
               script_name = "12_supp_tables_S5_to_S7.R")
    log_output(file_path = table_s6_gsea_enrichment_output, file_type = "table", figure_id = "Table S6 (GSEA enrichment full)",
               description = "GSEA enrichment full results table (auto-scan, standardized schema)",
               script_name = "12_supp_tables_S5_to_S7.R")
    if (file.exists(table_s6_fig1c_output)) {
      log_output(file_path = table_s6_fig1c_output, file_type = "table", figure_id = "Table S6_Fig1C",
                 description = "GSEA results for Figure 1C (filtered, top 35)",
                 script_name = "12_supp_tables_S5_to_S7.R")
    }
    log_output(file_path = table_s7_output, file_type = "table", figure_id = "Table S7",
               description = "Drug–gene interaction full table for locked signature (preferred TSV inputs; fallback auto-scan)",
               script_name = "12_supp_tables_S5_to_S7.R")
  }, error = function(e) {
    cat(sprintf("WARNING: log_output failed: %s\n", e$message))
  })
}

# =====================================================================
# Input resolution report
# =====================================================================
cat("\n=== Input Resolution Report ===\n")
cat(sprintf("S5 expr resolved to: %s\n", ifelse(is.na(expr_65204_path), "Not found", expr_65204_path)))
cat(sprintf("S5 pheno resolved to: %s\n", ifelse(is.na(pheno_65204_path), "Not found", pheno_65204_path)))
cat(sprintf("S6 files found: N=%d\n", length(gsea_files)))
cat(sprintf("S7 network inputs: preferred TSV present? %s\n", ifelse(file.exists(genes_tsv) && file.exists(drugs_tsv) && file.exists(interactions_tsv), "Yes", "No")))

cat("\n=== Supplementary Tables S5 to S7 generation completed ===\n")

# Memory cleanup (lightweight; do not over-delete in case debugging is needed)
rm(list = ls(pattern = "^(table_s[5-7]|expr_|pheno_|gsea_|dgi_|combined|boot_|score_|p_)", all.names = TRUE))
gc()