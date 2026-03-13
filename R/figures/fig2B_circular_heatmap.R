# =========================================================================
# Function: make_fig2B_circular_heatmap
# Purpose: Generate Figure 2B - Circular heatmap of 11-gene signature
# Inputs:
#   - cohort_manifest: Data frame with cohort information
#   - signature_genes: Vector of 11 signature genes
#   - output_dir: Directory to save output files
#   - tables_dir: Directory to save statistics files
# Outputs:
#   - output_dir/Fig2B_heatmap.pdf
#   - output_dir/Fig2B_heatmap.png (optional)
#   - tables_dir/Fig2B_gene_coverage.csv
#   - tables_dir/Fig2B_sample_order.csv
# =========================================================================

make_fig2B_circular_heatmap <- function(
  cohort_manifest,
  signature_genes,
  output_dir,
  tables_dir,
  base_dir
) {
  # Load required libraries
  library(ComplexHeatmap)
  library(circlize)
  library(tidyverse)
  library(yaml)
  library(grid)
  library(viridisLite)
  library(dendextend)
  
  # Source theme setup for color palettes
  source(file.path(base_dir, "R", "00_theme_setup.R"))
  
  # Create output directories if not exist
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ==============================================================================
  # Patch: ensure anchor cohorts exist in manifest (minimal intervention)
  # ==============================================================================
  
  required_rows <- tibble::tribble(
    ~cohort_id,  ~tissue_group, ~is_anchor,
    "GSE65204",  "Nasal",       TRUE,
    "GSE201955", "Airway",      TRUE,
    "GSE45111",  "Airway",      TRUE
  )
  
  # If cohort_manifest doesn't have these cohort_ids, add them; if present but tissue_group/is_anchor are empty, fill them in
  cohort_manifest <- cohort_manifest %>% 
    dplyr::full_join(required_rows, by = "cohort_id", suffix = c("", ".req")) %>% 
    dplyr::mutate(
      tissue_group = dplyr::coalesce(.data$tissue_group, .data$tissue_group.req),
      is_anchor    = dplyr::coalesce(.data$is_anchor,    .data$is_anchor.req)
    ) %>% 
    dplyr::select(-dplyr::ends_with(".req"))
  
  # Use color palettes from theme setup
  trae_colors_tissue <- trae_colors_tissue
  trae_colors_cohort <- trae_colors_cohort
  trae_colors_heatmap_zscore <- trae_colors_heatmap_zscore
  
  # ==============================================================================
  # Step 1: Write read and align functions
  # ==============================================================================
  
  # Read and align function
  read_and_align <- function(file_path, signature_genes) {
    # Read expression matrix
    cat("[Debug] Attempting to read file: ", file_path, "\n")
    if (!file.exists(file_path)) {
      cat("[Error] File not found: ", file_path, "\n")
      # For non-existent files, return NA matrix
      aligned_mat <- matrix(NA, nrow = length(signature_genes), 
                            ncol = 0, 
                            dimnames = list(signature_genes, NULL))
      return(aligned_mat)
    }
    
    cat("[Debug] File exists, starting to read...\n")
    # Read RDS file
    mat <- readRDS(file_path)
    cat("[Debug] Read successful, matrix dimensions: ", nrow(mat), "×", ncol(mat), "\n")
    
    # Create strictly 11-row empty matrix, fill NA for undetected genes
    aligned_mat <- matrix(NA, nrow = length(signature_genes), 
                          ncol = ncol(mat), 
                          dimnames = list(signature_genes, colnames(mat)))
    
    # Fill in detected genes
    common_genes <- intersect(rownames(mat), signature_genes)
    cat("[Debug] Number of common genes: ", length(common_genes), "\n")
    if (length(common_genes) > 0) {
      aligned_mat[common_genes, ] <- mat[common_genes, ]
    }
    
    cat("[Debug] Aligned matrix dimensions: ", nrow(aligned_mat), "×", ncol(aligned_mat), "\n")
    return(aligned_mat)
  }
  
  # Read from diet expression file and calculate z-score
  read_diet_expr_and_calculate_zscore <- function(cohort_id, signature_genes, base_dir) {
    # Construct diet expression file path
    diet_expr_path <- file.path(base_dir, "data", "processed_diet", paste0(cohort_id, "_diet_expr.rds"))
    cat("[Debug] Attempting to read diet expression file: ", diet_expr_path, "\n")
    
    if (!file.exists(diet_expr_path)) {
      cat("[Error] Diet expression file not found: ", diet_expr_path, "\n")
      # For non-existent files, return NA matrix
      aligned_mat <- matrix(NA, nrow = length(signature_genes), 
                            ncol = 0, 
                            dimnames = list(signature_genes, NULL))
      return(aligned_mat)
    }
    
    cat("[Debug] Diet file exists, starting to read...\n")
    # Read RDS file
    mat <- readRDS(diet_expr_path)
    cat("[Debug] Read successful, matrix dimensions: ", nrow(mat), "×", ncol(mat), "\n")
    
    # Create strictly 11-row empty matrix, fill NA for undetected genes
    aligned_mat <- matrix(NA, nrow = length(signature_genes), 
                          ncol = ncol(mat), 
                          dimnames = list(signature_genes, colnames(mat)))
    
    # Fill in detected genes
    common_genes <- intersect(rownames(mat), signature_genes)
    cat("[Debug] Number of common genes: ", length(common_genes), "\n")
    if (length(common_genes) > 0) {
      aligned_mat[common_genes, ] <- mat[common_genes, ]
    }
    
    # Calculate z-score
    cat("[Debug] Calculating z-score...\n")
    scaled_mat <- safe_zscore(aligned_mat)
    
    cat("[Debug] Aligned and standardized matrix dimensions: ", nrow(scaled_mat), "×", ncol(scaled_mat), "\n")
    return(scaled_mat)
  }
  
  # Extremely rigorous safe Z-score standardization function
  safe_zscore <- function(mat) {
    scaled_mat <- mat
    for (i in 1:nrow(mat)) {
      row_data <- mat[i, ]
      # Find data with actual probes (avoid True Missing NA)
      non_na_idx <- !is.na(row_data)
      
      if (sum(non_na_idx) > 1) {
        valid_data <- row_data[non_na_idx]
        # Precisely check variance: if there are probes but variance is 0 (e.g., CST1 in blood), force Z-score to zero (representing population baseline white/light blue)
        if (var(valid_data) == 0) {
          scaled_mat[i, non_na_idx] <- 0  
        } else {
          # Normal Z-score calculation
          scaled_mat[i, non_na_idx] <- scale(valid_data)
        }
      } else if (sum(non_na_idx) == 1) {
        scaled_mat[i, non_na_idx] <- 0
      }
    }
    # Never batch process is.na() here, preserve real NA to draw gray missing blocks!
    return(scaled_mat)
  }
  
  # ==============================================================================
  # Step 2: Data reading, alignment, and standardization
  # ==============================================================================
  cat("\n=== Step 2: Batch processing all cohorts ===\n")
  
  processed_cohorts <- list()
  coverage_stats <- list()
  
  for (i in 1:nrow(cohort_manifest)) {
    cohort <- cohort_manifest[i, ]
    cat("[Processing] Cohort: ", cohort$cohort_id, " (", cohort$tissue_group, ")\n")
    
    # Read from diet expression file and calculate z-score
    mat <- read_diet_expr_and_calculate_zscore(cohort$cohort_id, signature_genes, base_dir)
    
    # Check if matrix is empty
    if (ncol(mat) == 0) {
      cat("[Warning] Cohort", cohort$cohort_id, "has no sample data, skipping processing\n")
      coverage_stats[[cohort$cohort_id]] <- 0
      # Store empty matrix
      processed_cohorts[[cohort$cohort_id]] <- list(
        mat = mat,
        tissue_group = cohort$tissue_group,
        is_anchor = cohort$is_anchor
      )
      next
    }
    
    # Calculate gene coverage
    coverage <- sum(!is.na(rowSums(mat)))
    coverage_stats[[cohort$cohort_id]] <- coverage
    cat("[Coverage] Detected genes: ", coverage, "/11\n")
    
    # Store processed data
    processed_cohorts[[cohort$cohort_id]] <- list(
      mat = mat,
      tissue_group = cohort$tissue_group,
      is_anchor = cohort$is_anchor
    )
    cat("[Debug] Cohort", cohort$cohort_id, "processing completed\n")
  }
  
  # ==============================================================================
  # Step 3: Generate Figure 2B (7 cohorts main figure)
  # ==============================================================================
  cat("\n=== Step 3: Generate Figure 2B ===\n")
  
  # Filter out non-anchor cohorts and exclude GSE152004
  manifest_7cohorts <- cohort_manifest %>% filter(!is_anchor, cohort_id != "GSE152004")
  cat("[Main figure] Number of cohorts: ", nrow(manifest_7cohorts), "\n")
  print(manifest_7cohorts)
  
  # Print main figure cohort count
  cat("[Main figure] Actual number of cohorts: ", nrow(manifest_7cohorts), "\n")
  
  # Check cohort data
  for (cohort_id in manifest_7cohorts$cohort_id) {
    mat <- processed_cohorts[[cohort_id]]$mat
    cat("[Main figure] Cohort", cohort_id, "matrix dimensions: ", nrow(mat), "×", ncol(mat), "\n")
  }
  
  # ==============================================================================
  # Step 4: Write main plotting function
  # ==============================================================================
  cat("\n=== Step 4: Write plotting function ===\n")
  
  plot_circos_heatmap <- function(manifest_subset, output_prefix, output_dir) {
    tryCatch({
      # Define NA color
      na_col <- "#D9D9D9"
      
      # Sort cohorts by tissue order
      cohort_order_subset <- c(
        manifest_subset %>% filter(tissue_group == "Nasal") %>% pull(cohort_id),
        manifest_subset %>% filter(tissue_group == "Airway") %>% pull(cohort_id),
        manifest_subset %>% filter(tissue_group == "Blood") %>% pull(cohort_id)
      )
      
      # Extract matrices in order
      ordered_mats_subset <- lapply(cohort_order_subset, function(cohort_id) {
        processed_cohorts[[cohort_id]]$mat
      })
      
      # Combine matrices
      combined_mat_subset <- do.call(cbind, ordered_mats_subset)
      
      # Check if matrix is empty
      if (ncol(combined_mat_subset) == 0) {
        cat("[Warning] Combined matrix is empty, skipping plotting: ", output_prefix, "\n")
        return(list(
          combined_mat = combined_mat_subset,
          cohort_labels = character(0),
          tissue_labels = character(0)
        ))
      }
      
      # Generate sample labels
      cohort_labels_subset <- unlist(lapply(cohort_order_subset, function(cohort_id) {
        mat <- processed_cohorts[[cohort_id]]$mat
        rep(cohort_id, ncol(mat))
      }))
      
      tissue_labels_subset <- unlist(lapply(cohort_order_subset, function(cohort_id) {
        tissue <- processed_cohorts[[cohort_id]]$tissue_group
        mat <- processed_cohorts[[cohort_id]]$mat
        rep(tissue, ncol(mat))
      }))
      
      # Check if label vectors are empty
      if (length(tissue_labels_subset) == 0) {
        cat("[Warning] Sample label vectors are empty, skipping plotting: ", output_prefix, "\n")
        return(list(
          combined_mat = combined_mat_subset,
          cohort_labels = cohort_labels_subset,
          tissue_labels = tissue_labels_subset
        ))
      }
      
      # 1. Force unique column names (resolve warnings from duplicate sample names)
      colnames(combined_mat_subset) <- paste0(cohort_labels_subset, "_", 1:ncol(combined_mat_subset))
      
      # Transpose matrix for circular heatmap
      mat_t <- t(combined_mat_subset)
      
      # [New optimization]: Reverse matrix columns!
      # Make longer gene names the first column (drawn on spacious outer circle), shorter ones the last column (drawn on narrow inner circle)
      mat_t <- mat_t[, rev(colnames(mat_t))]
      
      # Prepare tissue split factor
      tissue_levels <- c("Nasal", "Airway", "Blood")
      tissue_factor <- factor(tissue_labels_subset, levels = tissue_levels)
      
      # ==========================================================
      # Core trick: Pre-calculate and pre-color cluster trees
      # ==========================================================
      cat("[Debug] Starting to pre-calculate cluster trees for three major tissues...\n")
      dend_list <- lapply(tissue_levels, function(tissue) {
        sub_mat <- mat_t[tissue_labels_subset == tissue, , drop = FALSE]
        sub_mat_dist <- sub_mat
        sub_mat_dist[is.na(sub_mat_dist)] <- 0
        d <- dist(sub_mat_dist)
        d[is.na(d)] <- 0 
        hc <- hclust(d)
        dend <- as.dendrogram(hc)
        branch_col <- trae_colors_tissue[tissue]
        dend <- dendextend::color_branches(dend, k = 1, col = branch_col)
        return(dend)
      })
      names(dend_list) <- tissue_levels
      cat("[Debug] Cluster tree pre-calculation completed!\n")

      # ==========================================================
      # Define common plotting execution body
      # ==========================================================
      draw_circos_plot <- function() {
        circos.clear()
        circos.par(start.degree = 90, gap.degree = c(2, 2, 12), track.margin = c(0.01, 0.01))
        
        # 1. [First circle]: 7 major cohort color bands (Cohort)
        # ⚠️ Core trick: Use cluster=TRUE to let the system automatically calculate and determine the global sample order!
        circos.heatmap(
          cohort_labels_subset,
          col = trae_colors_cohort,
          split = tissue_factor,
          cluster = TRUE,            # <--- Use system default clustering to determine order
          dend.side = "none",        # <--- But hide it because we will draw it manually on the innermost circle!
          track.height = 0.03,
          rownames.side = "none"
        )
        
        # 2. [Second circle]: Three major tissue color bands (Tissue)
        circos.heatmap(
          tissue_labels_subset,
          col = trae_colors_tissue,
          track.height = 0.03,
          rownames.side = "none"
        )
        
        # Add Nasal/Airway/Blood large labels on the periphery of the tissue circle
        circos.track(track.index = 2, panel.fun = function(x, y) {
          circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + mm_y(1),
                      CELL_META$sector.index, facing = "bending.inside", font = 2, cex = 1.3)
        }, bg.border = NA)
        
        # 3. [Third circle]: Main heatmap (samples have been correctly ordered by the first circle)
        circos.heatmap(
          mat_t,
          col = trae_colors_heatmap_zscore,
          track.height = 0.35,
          rownames.side = "none"     # Absolutely no messy sample names on the outside
        )
        
        # 4. Precisely draw gene names at the 12-degree gap in the main heatmap
        gene_names <- colnames(mat_t)
        circos.track(track.index = 3, panel.fun = function(x, y) {
          if (CELL_META$sector.index == "Nasal") {
            n_genes <- length(gene_names)
            for (i in seq_along(gene_names)) {
              # [Matching coordinate fix]:
              # circos.heatmap Y-axis is from 0(inner) to 11(outer)
              # First column (long gene names) should be drawn on the outermost circle (y=10.5)
              # Eleventh column (short gene names) should be drawn on the innermost circle (y=0.5)
              y_pos <- n_genes - i + 0.5
              circos.text(CELL_META$cell.xlim[1], y_pos, gene_names[i],
                          facing = "inside", niceFacing = TRUE,
                          adj = c(1.1, 0.5), cex = 0.75, font = 2)
            }
          }
        }, bg.border = NA)

        # 5. [Fourth circle]: Manually draw cluster trees!
        # Get the height of the tallest tree to unify Y-axis scale
        max_height <- max(sapply(dend_list, function(x) attr(x, "height")))
        circos.track(ylim = c(0, max_height), panel.fun = function(x, y) {
          sector.index <- CELL_META$sector.index
          dend <- dend_list[[sector.index]]
          # facing="outside" means the tree leaves grow outward, perfectly fitting at the bottom of the third circle heatmap!
          circos.dendrogram(dend, max_height = max_height, facing = "outside")
        }, track.height = 0.16, bg.border = NA)
        
        # Dynamically get the truly existing cohorts in the current chart (7 for main figure)
        active_cohorts <- unique(cohort_labels_subset)
        
        # 6. Reconstruct centered legend
        # Adjust legend size based on device type
        if (dev.cur() == 1) {  # 1 represents default device (usually PNG)
          font_size <- 10
          gap_size <- unit(4, "mm")
        } else {  # Other devices (like PDF)
          font_size <- 8
          gap_size <- unit(3, "mm")
        }
        
        lg_zscore <- Legend(title = "Expression (Z-score)", col_fun = trae_colors_heatmap_zscore, direction = "horizontal", title_position = "topcenter", title_gp = gpar(fontface = "bold", fontsize = font_size), labels_gp = gpar(fontsize = font_size))
        
        lg_tissue <- Legend(title = "Tissue Group", labels = names(trae_colors_tissue), legend_gp = gpar(fill = trae_colors_tissue), ncol = 1, title_gp = gpar(fontface = "bold", fontsize = font_size), labels_gp = gpar(fontsize = font_size))
        
        # Only draw legends for existing active_cohorts
        lg_cohort <- Legend(title = "Cohort", labels = active_cohorts, legend_gp = gpar(fill = trae_colors_cohort[active_cohorts]), ncol = 2, title_gp = gpar(fontface = "bold", fontsize = font_size), labels_gp = gpar(fontsize = font_size))
        
        pd <- packLegend(lg_zscore, lg_tissue, lg_cohort, direction = "vertical", gap = gap_size)
        draw(pd, x = unit(0.5, "npc"), y = unit(0.5, "npc"))  # Absolutely centered at the center
        
        circos.clear()
      }

      # ==========================================================
      # Execute PDF and PNG output
      # ==========================================================
      
      # Draw PDF
      pdf_file <- file.path(output_dir, paste0(output_prefix, ".pdf"))
      tryCatch({
        # Adjust PDF output font size to match PNG
        pdf(pdf_file, width = 12, height = 12, onefile = TRUE, family = "Helvetica", paper = "a4", pointsize = 8)
        draw_circos_plot()
        dev.off()
        cat("[Save] PDF file: ", pdf_file, "\n")
      }, error = function(e) {
        cat("[Error] Failed to generate PDF: ", e$message, "\n")
        try(dev.off(), silent = TRUE)
      })

      # Draw PNG
      png_file <- file.path(output_dir, paste0(output_prefix, ".png"))
      tryCatch({
        png(png_file, width = 3600, height = 3600, res = 300)
        draw_circos_plot()
        dev.off()
        cat("[Save] PNG file: ", png_file, "\n")
      }, error = function(e) {
        cat("[Error] Failed to generate PNG: ", e$message, "\n")
        try(dev.off(), silent = TRUE)
      })

      return(list(
        combined_mat = combined_mat_subset,
        cohort_labels = cohort_labels_subset,
        tissue_labels = tissue_labels_subset
      ))
    }, error = function(e) {
      cat("[Error] Failed to generate heatmap: ", e$message, "\n")
      return(list(
        combined_mat = matrix(NA, nrow = length(signature_genes), ncol = 0),
        cohort_labels = character(0),
        tissue_labels = character(0)
      ))
    })
  }
  
  # Generate main figure
  cat("[Main figure] Starting to generate Figure 2B...\n")
  fig2b_result <- tryCatch({
    plot_circos_heatmap(
      manifest_subset = manifest_7cohorts,
      output_prefix = "Fig2B_heatmap",
      output_dir = output_dir
    )
  }, error = function(e) {
    cat("[Error] Failed to generate Figure 2B: ", e$message, "\n")
    return(list(
      combined_mat = matrix(NA, nrow = length(signature_genes), ncol = 0),
      cohort_labels = character(0),
      tissue_labels = character(0)
    ))
  })
  cat("[Main figure] Figure 2B generation completed\n")
  cat("[Main figure] Result matrix dimensions: ", nrow(fig2b_result$combined_mat), "×", ncol(fig2b_result$combined_mat), "\n")
  
  # ==============================================================================
  # Step 3b: Generate Figure S2A (7 replication-only + 3 anchors)
  # ==============================================================================
  
  cat("\n=== Step 3b: Generate Figure S2A (7 + 3 anchors) ===\n")
  
  # Check if cohort_manifest contains all necessary cohorts
  cat("[Supplementary] Checking cohort_manifest...\n")
  cat("[Supplementary] Number of cohorts in cohort_manifest: ", nrow(cohort_manifest), "\n")
  print(cohort_manifest)
  
  # Include anchor + non-anchor, but exclude GSE152004
  manifest_s2a <- cohort_manifest %>% filter(cohort_id != "GSE152004")
  cat("[Supplementary] Number of cohorts: ", nrow(manifest_s2a), "\n")
  print(manifest_s2a %>% dplyr::arrange(tissue_group, is_anchor, cohort_id))
  
  # Use existing supplement directory
  output_dir_supp <- file.path(base_dir, "output", "supplement")
  cat("[Supplementary] Output directory: ", output_dir_supp, "\n")
  dir.create(output_dir_supp, recursive = TRUE, showWarnings = FALSE)
  
  # Check if directory was created successfully
  if (dir.exists(output_dir_supp)) {
    cat("[Supplementary] Directory created successfully\n")
  } else {
    cat("[Supplementary] Failed to create directory\n")
  }
  
  # Generate FigS2A
  cat("[Supplementary] Starting to generate Figure S2A...\n")
  figs2a_result <- tryCatch({
    plot_circos_heatmap(
      manifest_subset = manifest_s2a,
      output_prefix   = "FigureS2A_heatmap",
      output_dir      = output_dir_supp
    )
  }, error = function(e) {
    cat("[Error] Failed to generate Figure S2A: ", e$message, "\n")
    return(list(
      combined_mat = matrix(NA, nrow = length(signature_genes), ncol = 0),
      cohort_labels = character(0),
      tissue_labels = character(0)
    ))
  })
  
  cat("[Supplementary] Figure S2A generation completed\n")
  cat("[Supplementary] Result matrix dimensions: ", nrow(figs2a_result$combined_mat), "×", ncol(figs2a_result$combined_mat), "\n")
  
  # ==============================================================================
  # Step 5: Export QC data
  # ==============================================================================
  cat("\n=== Step 5: Export QC data ===\n")
  
  # Generate coverage statistics
  cat("[Debug] Starting to generate coverage statistics...\n")
  tryCatch({
    coverage_df <- data.frame(
      cohort_id = names(coverage_stats),
      coverage = unlist(coverage_stats),
      tissue_group = cohort_manifest$tissue_group[match(names(coverage_stats), cohort_manifest$cohort_id)],
      is_anchor = cohort_manifest$is_anchor[match(names(coverage_stats), cohort_manifest$cohort_id)]
    )
    
    coverage_file <- file.path(tables_dir, "Fig2B_gene_coverage.csv")
    cat("[Debug] Coverage statistics file path: ", coverage_file, "\n")
    write.csv(coverage_df, coverage_file, row.names = FALSE)
    cat("[Save] Coverage statistics: ", coverage_file, "\n")
  }, error = function(e) {
    cat("[Error] Failed to generate coverage statistics: ", e$message, "\n")
  })
  
  # Export sample order
  cat("[Debug] Starting to generate sample order file...\n")
  tryCatch({
    sample_order_df <- data.frame(
      sample_id = colnames(fig2b_result$combined_mat),
      cohort_id = fig2b_result$cohort_labels,
      tissue_group = fig2b_result$tissue_labels
    )
    
    sample_order_file <- file.path(tables_dir, "Fig2B_sample_order.csv")
    write.csv(sample_order_df, sample_order_file, row.names = FALSE)
    cat("[Save] Sample order: ", sample_order_file, "\n")
  }, error = function(e) {
    cat("[Error] Failed to generate sample order file: ", e$message, "\n")
  })
  
  # ==============================================================================
  # Export Figure S2A QC data
  # ==============================================================================
  cat("\n=== Export Figure S2A QC data ===\n")
  
  # Export Figure S2A gene coverage
  coverage_file_s2a <- file.path(tables_dir, "FigureS2A_gene_coverage.csv")
  write.csv(coverage_df, coverage_file_s2a, row.names = FALSE)
  cat("[Save] Figure S2A coverage statistics: ", coverage_file_s2a, "\n")
  
  # Export Figure S2A sample order
  tryCatch({
    sample_order_df_s2a <- data.frame(
      sample_id = colnames(figs2a_result$combined_mat),
      cohort_id = figs2a_result$cohort_labels,
      tissue_group = figs2a_result$tissue_labels
    )
    
    sample_order_file_s2a <- file.path(tables_dir, "FigureS2A_sample_order.csv")
    write.csv(sample_order_df_s2a, sample_order_file_s2a, row.names = FALSE)
    cat("[Save] Figure S2A sample order: ", sample_order_file_s2a, "\n")
  }, error = function(e) {
    cat("[Error] Failed to generate Figure S2A sample order file: ", e$message, "\n")
  })
  
  cat("\n=== Figure 2B and Figure S2A generation completed ===\n")
  
  return(list(
    output_pdf = file.path(output_dir, "Fig2B_heatmap.pdf"),
    output_png = file.path(output_dir, "Fig2B_heatmap.png"),
    coverage_file = coverage_file,
    sample_order_file = sample_order_file,
    output_pdf_s2a = file.path(output_dir_supp, "FigureS2A_heatmap.pdf"),
    output_png_s2a = file.path(output_dir_supp, "FigureS2A_heatmap.png"),
    coverage_file_s2a = coverage_file_s2a,
    sample_order_file_s2a = sample_order_file_s2a
  ))
}
