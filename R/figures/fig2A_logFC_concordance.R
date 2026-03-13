# =========================================================================
# Function: make_fig2A_logfc_concordance
# Purpose: Generate Figure 2A - log2FC concordance plot
# Inputs:
#   - gse152004_expr_path: Path to GSE152004 expression file
#   - gse152004_pheno_path: Path to GSE152004 phenotype file
#   - gse40888_expr_path: Path to GSE40888 expression file
#   - gse40888_pheno_path: Path to GSE40888 phenotype file
#   - gse40888_probe_map_path: Path to GSE40888 probe to gene mapping file
#   - output_dir: Directory to save output files
#   - tables_dir: Directory to save statistics files
# Outputs:
#   - output_dir/Fig2A_log2FC.pdf
#   - output_dir/Fig2A_log2FC.png (optional)
#   - tables_dir/Fig2A_concordance_stats.csv
# =========================================================================

make_fig2A_logfc_concordance <- function(
  gse152004_expr_path,
  gse152004_pheno_path,
  gse40888_expr_path,
  gse40888_pheno_path,
  output_dir,
  tables_dir,
  base_dir
) {
  # Load required libraries
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  
  # Source theme setup for color palettes
  source(file.path(base_dir, "R", "00_theme_setup.R"))
  
  # Define 11-gene signature
  signature_genes <- c("CLCA1", "SERPINB2", "CPA3", "CCR3", "HDC", 
                      "MUC5AC", "IL4", "MUC5B", "IL13", "CCL24", "POSTN")
  
  # Step 1: Read GSE152004 data
  cat("Step 1: Reading GSE152004 data...\n")
  expr_152004 <- readRDS(gse152004_expr_path)
  
  # Read T2 labels for GSE152004
  gse152004_labels_path <- file.path(base_dir, "data", "derived", "GSE152004_T2_labels.csv")
  cat("Reading GSE152004 T2 labels...\n")
  gse152004_labels <- read.csv(gse152004_labels_path)
  
  # Step 2: Calculate logFC for GSE152004 (T2-high vs T2-low)
  cat("Step 2: Calculating logFC for GSE152004...\n")
  case_samples <- gse152004_labels$Sample_ID[gse152004_labels$T2_Status == "T2-high"]
  control_samples <- gse152004_labels$Sample_ID[gse152004_labels$T2_Status == "T2-low"]
  
  # Calculate mean expression for each gene in case and control
  case_means <- rowMeans(expr_152004[, case_samples])
  control_means <- rowMeans(expr_152004[, control_samples])
  
  # Calculate logFC
  logfc <- log2(case_means / control_means)
  
  # Create data frame and standardize gene symbols
  df_152004 <- data.frame(
    gene_symbol = toupper(trimws(names(logfc))),
    logFC_152004 = as.numeric(logfc)
  )
  
  # Remove NA gene symbols
  df_152004 <- df_152004[!is.na(df_152004$gene_symbol) & df_152004$gene_symbol != "", ]
  
  # Fail-fast assertions
  stopifnot(nrow(df_152004) > 5000)
  cat(paste("Calculated logFC for", nrow(df_152004), "genes from GSE152004\n"))
  
  # Step 3: Read GSE40888 expression and phenotype data
  cat("Step 3: Reading GSE40888 data...\n")
  expr_40888 <- readRDS(gse40888_expr_path)
  pheno_40888 <- readRDS(gse40888_pheno_path)
  
  # Step 4: Calculate logFC for GSE40888 (Atopic wheeze vs Healthy non-wheeze)
  cat("Step 4: Calculating logFC for GSE40888...\n")
  
  # Print group labels and counts
  cat("GSE40888 group labels and counts:\n")
  group_counts <- table(pheno_40888$std_group_label)
  print(group_counts)
  
  case_samples_40888 <- pheno_40888$std_sample_id[pheno_40888$std_group_label == "Atopic wheeze"]
  control_samples_40888 <- pheno_40888$std_sample_id[pheno_40888$std_group_label == "Healthy non-wheeze"]
  
  cat(paste("Case samples (Atopic wheeze):", length(case_samples_40888), "\n"))
  cat(paste("Control samples (Healthy non-wheeze):", length(control_samples_40888), "\n"))
  
  # Calculate mean expression for each gene in case and control
  case_means_40888 <- rowMeans(expr_40888[, case_samples_40888])
  control_means_40888 <- rowMeans(expr_40888[, control_samples_40888])
  
  # Calculate logFC
  logfc_40888 <- log2(case_means_40888 / control_means_40888)
  
  # Create data frame and standardize gene symbols
  df_40888 <- data.frame(
    gene_symbol = toupper(trimws(names(logfc_40888))),
    logFC_40888 = as.numeric(logfc_40888)
  )
  
  # Create means file with finite checks
  means_40888_file <- file.path(tables_dir, "Fig2A_means_40888_allgenes.csv")
  means_40888_df <- data.frame(
    gene_symbol = toupper(trimws(names(case_means_40888))),
    case_mean = as.numeric(case_means_40888),
    control_mean = as.numeric(control_means_40888),
    case_mean_is_finite = is.finite(case_means_40888),
    control_mean_is_finite = is.finite(control_means_40888)
  )
  write.csv(means_40888_df, means_40888_file, row.names = FALSE)
  cat(paste("GSE40888 means saved to", means_40888_file, "\n"))
  
  # Remove NA gene symbols
  df_40888 <- df_40888[!is.na(df_40888$gene_symbol) & df_40888$gene_symbol != "", ]
  
  cat(paste("Calculated logFC for", nrow(df_40888), "genes from GSE40888\n"))
  
  # Fail-fast assertions
  stopifnot(nrow(df_40888) > 5000)
  cat(paste("Mapped and aggregated", nrow(df_40888), "genes for GSE40888\n"))
  
  # Step 6: Merge data from both datasets
  cat("Step 6: Merging data...\n")
  merged_data <- merge(df_152004, df_40888, by = "gene_symbol")
  n_common_genes <- nrow(merged_data)
  
  # Fail-fast assertion
  stopifnot(n_common_genes > 3000)
  cat(paste("Merged data contains", n_common_genes, "common genes\n"))
  
  # Add is_signature column
  merged_data$is_signature <- ifelse(merged_data$gene_symbol %in% signature_genes, TRUE, FALSE)
  
  # Step 7: Save audit files
  cat("Step 7: Saving audit files...\n")
  
  # Save GSE152004 logFC
  logfc_152004_file <- file.path(tables_dir, "Fig2A_logFC_152004_allgenes.csv")
  write.csv(df_152004, logfc_152004_file, row.names = FALSE)
  cat(paste("GSE152004 logFC saved to", logfc_152004_file, "\n"))
  
  # Save GSE40888 logFC
  logfc_40888_file <- file.path(tables_dir, "Fig2A_logFC_40888_allgenes.csv")
  write.csv(df_40888, logfc_40888_file, row.names = FALSE)
  cat(paste("GSE40888 logFC saved to", logfc_40888_file, "\n"))
  
  # Save plot input
  plot_input_file <- file.path(tables_dir, "Fig2A_plot_input.csv")
  write.csv(merged_data, plot_input_file, row.names = FALSE)
  cat(paste("Plot input saved to", plot_input_file, "\n"))
  
  # Step 8: Classify genes for plotting
  merged_data$class <- ifelse(merged_data$gene_symbol %in% signature_genes, "signature", "background")
  
  # Find signature genes in data
  signature_in_data <- intersect(signature_genes, merged_data$gene_symbol)
  cat(paste("Found", length(signature_in_data), "signature genes in merged data\n"))
  
  # Create highlight data
  highlight_data <- merged_data[merged_data$class == "signature", ]
  
  # Check for missing MUC5AC
  missing_genes <- setdiff(signature_genes, signature_in_data)
  cat(paste("Missing signature genes:", paste(missing_genes, collapse = ", "), "\n"))
  
  # Step 9: Create plot
  cat("Step 9: Creating Figure 2A...\n")
  
  # Calculate dynamic X-axis limit based on signature genes
  if (nrow(highlight_data) > 0) {
    max_x <- ceiling(max(highlight_data$logFC_152004, na.rm = TRUE)) + 1
  } else {
    max_x <- 4  # Default value
  }
  cat(paste("Dynamic X-axis limit calculated:", max_x, "\n"))
  
  # Create ggplot
  plot <- ggplot() +
    # Background points (only common genes)
    geom_point(
      data = merged_data[merged_data$class == "background", ],
      aes(x = logFC_152004, y = logFC_40888),
      color = "grey85",
      alpha = 0.4,
      size = 0.8
    ) +
    # Density contours (only common genes)
    stat_density_2d(
      data = merged_data[merged_data$class == "background", ],
      aes(x = logFC_152004, y = logFC_40888),
      color = "grey50",
      bins = 12,
      linewidth = 0.6,
      alpha = 0.8
    ) +
    # Signature points
    geom_point(
      data = highlight_data,
      aes(x = logFC_152004, y = logFC_40888),
      color = "#E64B35",
      alpha = 1.0,
      size = 3.5
    ) +
    # Gene labels
    geom_text_repel(
      data = highlight_data,
      aes(x = logFC_152004, y = logFC_40888, label = gene_symbol),
      size = 4,
      fontface = "bold",
      color = "black",
      bg.color = "white",
      bg.r = 0.15,
      box.padding = 0.3,
      max.overlaps = 20,
      segment.color = "grey50",
      segment.size = 0.3
    ) +
    # Reference lines
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.6) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.6) +
    # Labels
    labs(
      x = "log₂FC in GSE152004 (Nasal: T2-high vs. T2-low)",
      y = "log₂FC in GSE40888 (PBMC: Atopic wheeze vs. Healthy non-wheeze)",
      title = "Cross-tissue Endotype Concordance",
      subtitle = ifelse("MUC5AC" %in% missing_genes, 
                       "Conservation of logFC directions from airway epithelium to peripheral blood (MUC5AC missing in GSE40888)",
                       "Conservation of logFC directions from airway epithelium to peripheral blood")
    ) +
    # Theme
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "grey60"),
      axis.title.x = element_text(size = 12, face = "bold"),
      axis.title.y = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, color = "grey80", linewidth = 0.5)
    ) +
    # Add coordinate limits to prevent outlier distortion
    coord_cartesian(xlim = c(-4, max_x), ylim = c(-4, 4))
  
  # Step 10: Save plot
  cat("Step 10: Saving Figure 2A...\n")
  
  # Create output directories if not exist
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Save as PDF
  output_file_pdf <- file.path(output_dir, "Fig2A_log2FC.pdf")
  ggsave(output_file_pdf, plot, width = 12, height = 10, dpi = 300)
  cat(paste("Figure 2A PDF saved to", output_file_pdf, "\n"))
  
  # Save as PNG
  output_file_png <- file.path(output_dir, "Fig2A_log2FC.png")
  ggsave(output_file_png, plot, width = 12, height = 10, dpi = 300)
  cat(paste("Figure 2A PNG saved to", output_file_png, "\n"))
  
  # Step 11: Calculate statistics
  cat("Step 11: Calculating concordance statistics...\n")
  if (nrow(merged_data) > 0) {
    # Calculate concordant considering NA values
    merged_data$concordant <- sign(merged_data$logFC_152004) == sign(merged_data$logFC_40888)
    total <- nrow(merged_data)
    # Exclude NA values when calculating concordance
    valid_data <- merged_data[!is.na(merged_data$concordant), ]
    valid_total <- nrow(valid_data)
    concordant <- sum(valid_data$concordant)
    rate <- (concordant / valid_total) * 100
    
    cat(paste("Total genes:", total, "\n"))
    cat(paste("Valid genes (non-NA):", valid_total, "\n"))
    cat(paste("Concordant genes:", concordant, "\n"))
    cat(paste("Concordance rate:", round(rate, 2), "%\n"))
    
    # Signature gene statistics
    if (nrow(highlight_data) > 0) {
      highlight_data$concordant <- sign(highlight_data$logFC_152004) == sign(highlight_data$logFC_40888)
      # Exclude NA values for signature genes
      valid_highlight <- highlight_data[!is.na(highlight_data$concordant), ]
      sig_total <- nrow(valid_highlight)
      sig_concordant <- sum(valid_highlight$concordant)
      sig_rate <- (sig_concordant / sig_total) * 100
      
      cat(paste("Signature genes total:", sig_total, "\n"))
      cat(paste("Signature genes concordant:", sig_concordant, "\n"))
      cat(paste("Signature concordance rate:", round(sig_rate, 2), "%\n"))
    }
    
    # Save statistics
    # Add group counts to stats
    group_counts_str <- paste(names(group_counts), group_counts, sep="=", collapse=", ")
    
    stats_df <- data.frame(
      metric = c("total_genes", "valid_genes", "concordant_genes", "total_concordance_rate",
                 "signature_genes_total", "signature_genes_concordant", "signature_concordance_rate",
                 "missing_genes", "gse40888_group_counts"),
      value = c(total, valid_total, concordant, rate,
                ifelse(exists("sig_total"), sig_total, 0),
                ifelse(exists("sig_concordant"), sig_concordant, 0),
                ifelse(exists("sig_rate"), sig_rate, 0),
                paste(missing_genes, collapse = ", "),
                group_counts_str)
    )
    
    stats_file <- file.path(tables_dir, "Fig2A_concordance_stats.csv")
    write.csv(stats_df, stats_file, row.names = FALSE)
    cat(paste("Statistics saved to", stats_file, "\n"))
    
  } else {
    cat("No common genes found between datasets\n")
  }
  
  # Add log output
  cat("Fig2A: n_152004=", nrow(df_152004), " n_40888=", nrow(df_40888), " n_common=", n_common_genes, "\n")
  
  # Sanity check for Fig2A_plot_input.csv
  cat("\nPerforming sanity check on Fig2A_plot_input.csv...\n")
  
  # Sample 20 rows
  sample_rows <- merged_data[sample(nrow(merged_data), min(20, nrow(merged_data))), ]
  cat("Sample 20 rows from Fig2A_plot_input.csv:\n")
  print(sample_rows)
  
  # Check signature genes
  signature_rows <- merged_data[merged_data$is_signature == TRUE, ]
  cat(paste("\nSignature genes found:", nrow(signature_rows), "\n"))
  print(signature_rows)
  
  # Check missing genes
  cat(paste("\nMissing genes:", paste(missing_genes, collapse = ", "), "\n"))
  
  cat("\nFigure 2A generation completed!\n")
  
  return(list(
    plot = plot,
    output_pdf = output_file_pdf,
    output_png = output_file_png,
    stats_file = stats_file,
    logfc_152004_file = logfc_152004_file,
    logfc_40888_file = logfc_40888_file,
    plot_input_file = plot_input_file
  ))
}
