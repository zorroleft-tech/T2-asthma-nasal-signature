# =========================================================================
# Function: make_fig2C_attenuation
# Purpose: Generate Figure 2C - Tissue gradient attenuation plot
# Inputs:
#   - output_dir: Directory to save output files
#   - tables_dir: Directory to save statistics files
# Outputs:
#   - output_dir/Fig2C_attenuation.pdf
#   - tables_dir/Fig2C_visual_guides.yaml
# =========================================================================

make_fig2C_attenuation <- function(
  output_dir,
  tables_dir,
  base_dir
) {
  # Load required libraries
  library(tidyverse)
  library(ggplot2)
  library(yaml)
  library(ggrepel)
  
  # Create output directories if not exist
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Source theme setup for color palettes
  source(file.path(base_dir, "R", "00_theme_setup.R"))
  
  # Use color palettes from theme setup
  trae_colors_cohort <- trae_colors_cohort
  
  # Define theme function
  theme_trae_publication <- function(base_size = 12, base_family = "sans") {
    theme_bw(base_size = base_size, base_family = base_family) +
      theme(
        # Title
        plot.title = element_text(size = base_size + 2, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = base_size, hjust = 0.5),
        
        # Axes
        axis.title = element_text(size = base_size, face = "bold"),
        axis.text = element_text(size = base_size - 2, color = "black"),
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.ticks = element_line(color = "black", linewidth = 0.5),
        
        # Legend
        legend.title = element_text(size = base_size, face = "bold"),
        legend.text = element_text(size = base_size - 2),
        legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3),
        legend.key = element_blank(),
        legend.position = "right",
        
        # Panel
        panel.background = element_rect(fill = "white", color = NA),
        panel.border = element_rect(color = "black", linewidth = 1),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        
        # Facet
        strip.background = element_rect(fill = "grey90", color = "black", linewidth = 0.5),
        strip.text = element_text(size = base_size, face = "bold"),
        
        # Margin
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
      )
  }
  
  # Step 1: Read Fig2C_cohort_effect_sizes.csv
  cat("Step 1: Reading Fig2C_cohort_effect_sizes.csv...\n")
  cohort_effect_file <- file.path(tables_dir, "Fig2C_cohort_effect_sizes.csv")
  
  if (!file.exists(cohort_effect_file)) {
    stop(paste("Fig2C_cohort_effect_sizes.csv not found:", cohort_effect_file))
  }
  
  cohort_stats <- read_csv(cohort_effect_file, show_col_types = FALSE)
  
  # Check if the file has the expected structure
  if ("cohens_d" %in% colnames(cohort_stats)) {
    # Use cohens_d column for y-axis
    cat("Using cohens_d column for y-axis\n")
  } else if ("d_value" %in% colnames(cohort_stats)) {
    # Old structure
    cohort_stats <- cohort_stats %>%
      rename(cohens_d = d_value)
    cat("Renamed d_value to cohens_d\n")
  } else {
    stop("Fig2C_cohort_effect_sizes.csv missing required column: cohens_d or d_value")
  }
  
  # Step 1.5: Unify tissue_group to 3 categories
  cat("Step 1.5: Unifying tissue_group to 3 categories...\n")
  
  # Define mapping function (vectorized with cohort_id)
  map_tissue_group <- function(tissue_raw, cohort_id) {
    mapply(function(t, cid) {
      if (is.na(t)) {
        return(NA)
      }
      tissue_lower <- tolower(t)
      if (tissue_lower %in% c("nasal", "nasal epithelium", "nasal swab", "upper airway")) {
        return("Nasal")
      } else if (tissue_lower %in% c("bronchial", "sputum", "tracheal", "airway", "lower airway", "lowerairway")) {
        return("Airway")
      } else if (tissue_lower %in% c("blood", "whole blood", "pbmc", "neutrophil", "bloodpbmc")) {
        return("Blood")
      } else {
        # For Other or unrecognized, map based on cohort_id
        if (cid %in% c("GSE118761", "GSE43696", "GSE103166")) {
          return("Airway")
        } else if (cid %in% c("GSE40888", "GSE230048", "GSE115770")) {
          return("Blood")
        } else if (cid %in% c("GSE123750")) {
          return("Nasal")
        } else {
          return(NA)
        }
      }
    }, tissue_raw, cohort_id)
  }
  
  # Apply mapping
  cohort_stats <- cohort_stats %>%
    mutate(tissue_group = map_tissue_group(tissue_group, cohort_id))
  
  # Hard gate validation
  if (any(is.na(cohort_stats$tissue_group)) || any(!cohort_stats$tissue_group %in% c("Nasal", "Airway", "Blood"))) {
    invalid_tissues <- cohort_stats %>%
      filter(is.na(tissue_group) | !tissue_group %in% c("Nasal", "Airway", "Blood"))
    stop(paste("Invalid tissue_group values found:", paste(invalid_tissues$cohort_id, collapse = ", ")))
  }
  
  cat(paste("Loaded", nrow(cohort_stats), "cohorts\n"))
  
  # Step 2: Calculate visual guides (median d per tissue group)
  cat("Step 2: Calculating visual guides...\n")
  
  visual_guides <- cohort_stats %>%
    group_by(tissue_group) %>%
    summarize(
      median_d = median(cohens_d, na.rm = TRUE),
      mean_d = mean(cohens_d, na.rm = TRUE),
      n_cohorts = n()
    ) %>%
    arrange(factor(tissue_group, levels = c("Nasal", "Airway", "Blood")))
  
  cat("Visual guides (median d):\n")
  print(visual_guides)
  
  # Step 3: Prepare data for plotting
  cat("Step 3: Preparing data for plotting...\n")
  
  # Ensure tissue_group factor order
  cohort_stats$tissue_group <- factor(cohort_stats$tissue_group, 
                                      levels = c("Nasal", "Airway", "Blood"))
  
  # Step 4: Create plot
  cat("Step 4: Creating Figure 2C...\n")
  
  p <- ggplot() +
    # Add horizontal line at y=0
    geom_hline(yintercept = 0, linetype = "solid", color = "grey50", linewidth = 0.5) +
    # Bottom layer: dashed line connecting medians
    geom_line(data = visual_guides, 
              aes(x = factor(tissue_group, levels = c("Nasal", "Airway", "Blood")), 
                  y = median_d, 
                  group = 1),
              linetype = "dashed", 
              color = "grey50", 
              linewidth = 1) +
    # Large points (summary): black diamonds representing median d
    geom_point(data = visual_guides, 
               aes(x = factor(tissue_group, levels = c("Nasal", "Airway", "Blood")), 
                   y = median_d),
               shape = 18, 
               size = 6, 
               color = "black") +
    # Scatter points (cohorts): jittered points for actual d values
    geom_jitter(data = cohort_stats, 
                aes(x = tissue_group, 
                    y = cohens_d, 
                    color = cohort_id),
                width = 0.1, 
                size = 6, 
                alpha = 0.9) +
    # Add cohort labels with repel
    geom_text_repel(data = cohort_stats, 
                    aes(x = tissue_group, 
                        y = cohens_d, 
                        label = cohort_id),
                    size = 6, 
                    color = "black", 
                    box.padding = 0.5, 
                    point.padding = 1.0, 
                    segment.color = "grey50", 
                    segment.size = 0.5) +
    # Color mapping
    scale_color_manual(values = trae_colors_cohort) +
    # Set y-axis limits based on actual data range
    coord_cartesian(ylim = c(-0.4, 0.4)) +
    # Add theme with larger base size
    theme_classic(base_size = 16) +
    # Remove legend
    theme(legend.position = "none") +
    # Add labels
    labs(
      x = "Tissue Group",
      y = "Cohen's d",
      title = "Tissue-Specific Attenuation of T2 Signature Effect Size",
      subtitle = "Visual guide based on median effect size (not pooled/meta-analysis)"
    )
  
  # Step 5: Save plot
  cat("Step 5: Saving Figure 2C...\n")
  
  # Save as PDF
  output_file_pdf <- file.path(output_dir, "Fig2C_attenuation.pdf")
  ggsave(output_file_pdf, p, width = 12, height = 10, dpi = 300)
  cat(paste("Figure 2C PDF saved to", output_file_pdf, "\n"))
  
  # Step 6: Save visual guides as YAML
  cat("Step 6: Saving visual guides...\n")
  
  # Convert to data frame first to avoid tidyverse split issues
  visual_guides_df <- as.data.frame(visual_guides)
  
  # Create list manually
  visual_guides_list <- list(
    visual_guides = list(
      "Nasal" = list(
        median_d = as.numeric(visual_guides_df$median_d[visual_guides_df$tissue_group == "Nasal"]),
        mean_d = as.numeric(visual_guides_df$mean_d[visual_guides_df$tissue_group == "Nasal"]),
        n_cohorts = as.integer(visual_guides_df$n_cohorts[visual_guides_df$tissue_group == "Nasal"])
      ),
      "Airway" = list(
        median_d = as.numeric(visual_guides_df$median_d[visual_guides_df$tissue_group == "Airway"]),
        mean_d = as.numeric(visual_guides_df$mean_d[visual_guides_df$tissue_group == "Airway"]),
        n_cohorts = as.integer(visual_guides_df$n_cohorts[visual_guides_df$tissue_group == "Airway"])
      ),
      "Blood" = list(
        median_d = as.numeric(visual_guides_df$median_d[visual_guides_df$tissue_group == "Blood"]),
        mean_d = as.numeric(visual_guides_df$mean_d[visual_guides_df$tissue_group == "Blood"]),
        n_cohorts = as.integer(visual_guides_df$n_cohorts[visual_guides_df$tissue_group == "Blood"])
      )
    ),
    note = "Visual guide based on median effect size (not pooled/meta-analysis)"
  )
  
  visual_guides_file <- file.path(tables_dir, "Fig2C_visual_guides.yaml")
  write_yaml(visual_guides_list, visual_guides_file)
  cat(paste("Visual guides saved to", visual_guides_file, "\n"))
  
  cat("\nFigure 2C generation completed!\n")
  
  return(list(
    plot = p,
    output_pdf = output_file_pdf,
    visual_guides_file = visual_guides_file
  ))
}
