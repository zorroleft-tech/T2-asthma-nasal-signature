# ==============================================================================
# Script: 00_theme_setup.R
# Purpose: Single Source of Truth for global themes, palettes, and mappings.
# Updated: 2026-02-24
# Note: ALL COMMENTS IN ENGLISH to prevent Encoding/UTF-8 crashes on reviewers' machines.
# Features: Optimized for Nature/Cell/Science (CNS) high-contrast visual standards.
# Architecture: 1 Discovery + 3 Clinical Anchors + 7 Replication Cohorts.
# ==============================================================================

library(ggplot2)
library(RColorBrewer)
library(scales)
# Note: Continuous heatmap palettes rely on the 'circlize' package (called directly in functions)

# ==============================================================================
# 1. Base Palette Libraries (Allergy Compliant Version)
# ==============================================================================

## 1.1 Main Figure Palette (Avoid Red/Green & Yellow/White clashes)
trae_colors_main <- c(
  "#E74C3C",  # Red (Case/High) - Kept for emphasis
  "#4CB5F5",  # Blue (Control/Low)
  "#f1a340",  # Orange (Intermediate) - Using Allergy Official Orange!
  "#542788",  # Purple (Normal/Train) - Using Allergy Official Purple! (Replaces Green)
  "#878787",  # Grey (Reference/Neutral)
  "#b35806"   # Dark Brown/Orange - Using Allergy Official Dark Orange!
)

## 1.2 Single-Cell Clustering Palette (35 colors, Nature/Cell style)
## Usage: Figure 3 (UMAP embeddings supporting up to 35 clusters)
trae_colors_celltype <- c(
  '#cb4936', '#ff0000', '#5ac1b3', '#1616c8', '#0D63A5',
  '#A12568', '#F499C1', '#F7C394', '#B2A157', '#ade87c',
  '#000000', '#03C4A1', '#bb4316', '#7382BC', '#F0E442',
  '#3B185F', '#0d6c0d', '#FEC260', '#FD7014', '#a67063',
  '#1B9B8A', '#D0EBE7', '#713045', '#F6E0EA', '#AD6D28',
  '#EAB67D', '#508CA4', '#EE4590', '#8B4513', '#00FFFF', 
  '#FF1493', '#7FFF00', '#4B0082', '#D2691E', '#B0C4DE'
)

## 1.3 Pathway & Network Palette (13 colors)
## Usage: Figure 1C (KEGG Network) & Figure 4D (Drug-Gene Network)
trae_colors_pathway <- c(
  '#007ABA', '#8B58A4', '#DF75AE', '#00B7CA',
  '#A8C7E9', '#E91E25', '#925047', '#F37121',
  '#FBB36E', '#F58D93', '#B4BE50', '#00A163', '#8CCA7C'
)

## 1.4 Diverging and Continuous Palettes
trae_colors_diverging <- colorRampPalette(c("#0000FF", "#FFFFFF", "#FF0000"))(100)
trae_colors_pvalue <- colorRampPalette(c("white", "#FF0000"))(100)


# ==============================================================================
# 2. Strict Semantic Mappings (Enforcing consistency across all figures)
# ==============================================================================

## 2.1 Tissue Compartment Palette (Colorblind safe: Red, Purple, Blue)
## Usage: Outer rings of circos plots, global tissue distinctions
trae_colors_tissue <- c(
  "Nasal"  = "#E91E25",  # Classic Red
  "Airway" = "#542788",  # Official Allergy Purple (Replaces the forbidden Green)
  "Blood"  = "#007ABA"   # Classic Blue
)

## 2.2 Comprehensive Cohort Palette (Updated to match colorblind safety)
trae_colors_cohort <- c(
  # Discovery (Dark Slate)
  "GSE152004" = "#2C3E50",
  
  # Clinical Biomarker Anchors (Vibrant, Colorblind Safe)
  "GSE65204"  = "#E74C3C",  # Anchor 1 (IgE) - Solid Red
  "GSE201955" = "#f1a340",  # Anchor 2 (FeNO) - Official Allergy Orange
  "GSE45111"  = "#007ABA",  # Anchor 3 (Sputum) - Solid Blue (Replaced Teal/Green)
  
  # Replication-only: Nasal (Soft Reds/Pinks)
  "GSE103166" = "#F1948A",  
  "GSE115770" = "#F5CBA7",  
  
  # Replication-only: Airway (Soft Purples - matching the new Airway color)
  "GSE118761" = "#d8daeb",  # Official Allergy light purple
  "GSE43696"  = "#998ec3",  # Official Allergy mid purple
  
  # Replication-only: Blood/PBMC (Soft Blues)
  "GSE123750" = "#AED6F1",  
  "GSE40888"  = "#5DADE2",  
  "GSE230048" = "#85C1E9"   
)

## 2.3 Heatmap Z-score Palette (Cell style)
## Usage: Figure 2B (Global expression landscape)
trae_colors_heatmap_zscore <- circlize::colorRamp2(
  c(-2, 0, 2), 
  c("#1f77b4", "white", "#d62728")
)

## 2.4 Phenotype and Clinical Labels Palette
trae_colors_group <- c(
  "Asthma"       = "#FF420E",
  "Control"      = "#4CB5F5",
  "T2-high"      = "#E91E25",
  "T2-low"       = "#007ABA",
  "IgE-high"     = "#E91E25",
  "IgE-low"      = "#007ABA",
  "FeNO-high"    = "#E91E25",
  "FeNO-low"     = "#007ABA",
  "Eos-high"     = "#E91E25",
  "Eos-low"      = "#007ABA"
)


# ==============================================================================
# 3. Global ggplot2 Themes
# ==============================================================================

## 3.1 Standard Publication Theme (Fully boxed, high contrast)
theme_trae_publication <- function(base_size = 12, base_family = "sans") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(size = base_size + 2, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = base_size, hjust = 0.5),
      
      axis.title = element_text(size = base_size, face = "bold"),
      axis.text = element_text(size = base_size - 2, color = "black"),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      
      legend.title = element_text(size = base_size, face = "bold"),
      legend.text = element_text(size = base_size - 2),
      legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3),
      legend.key = element_blank(),
      legend.position = "right",
      
      panel.background = element_rect(fill = "white", color = NA),
      panel.border = element_rect(color = "black", linewidth = 1),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      
      strip.background = element_rect(fill = "grey90", color = "black", linewidth = 0.5),
      strip.text = element_text(size = base_size, face = "bold"),
      
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
    )
}

# Global theme setting (backward compatibility)
theme_allergy <- function() {
  theme_trae_publication(base_size = 10)
}

# Publication-optimized theme (backward compatibility)
theme_pub <- function() {
  theme_trae_publication(base_size = 9)
}

# Color palette settings (backward compatibility)
color_palette <- trae_colors_main
color_gradient <- colorRampPalette(c("#4CB5F5", "#FFBB00"))
heatmap_palette <- colorRampPalette(c("#4CB5F5", "#FFFFFF", "#FF420E"))
cell_type_palette <- trae_colors_celltype[1:7]
names(cell_type_palette) <- c("T cell", "B cell", "Macrophage", "Dendritic cell", "Eosinophil", "Neutrophil", "Epithelial cell")
cohort_palette <- trae_colors_cohort


# ==============================================================================
# 4. Helper Functions for Export and Preview
# ==============================================================================

## 4.1 Universal figure saving function (Supports simultaneous PDF/PNG generation)
save_trae_figure <- function(plot, filename, width = 180, height = 120, 
                             dpi = 300, format = c("pdf", "png")) {
  output_dir <- dirname(filename)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  if ("pdf" %in% format) {
    ggsave(
      filename = paste0(tools::file_path_sans_ext(filename), ".pdf"),
      plot = plot, width = width, height = height, units = "mm", device = "pdf"
    )
  }
  
  if ("png" %in% format) {
    ggsave(
      filename = paste0(tools::file_path_sans_ext(filename), ".png"),
      plot = plot, width = width, height = height, units = "mm", dpi = dpi, device = "png"
    )
  }
  cat(sprintf("Success: Figure saved to %s\n", filename))
}

## 4.2 Palette Previewer
preview_colors <- function(colors, labels = NULL) {
  n <- length(colors)
  if (is.null(labels)) {
    labels <- if (!is.null(names(colors))) names(colors) else paste0("Color ", 1:n)
  }
  
  df <- data.frame(x = 1:n, y = 1, label = labels, color = unname(colors))
  
  p <- ggplot(df, aes(x = x, y = y, fill = color)) +
    geom_tile(color = "white", linewidth = 2) +
    geom_text(aes(label = label), angle = 90, hjust = 0.5, size = 3) +
    scale_fill_identity() +
    theme_void() +
    labs(title = paste0("TRAE Theme Palette Preview (n=", n, ")")) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
    )
  
  print(p)
  return(invisible(df))
}


# ==============================================================================
# 5. Theme Setting and Export Functions
# ==============================================================================

# Theme setting function
set_plot_theme <- function() {
  # Set global theme
  theme_set(theme_allergy())
  
  # Set global color palette
  options(ggplot2.discrete.colour = color_palette)
  options(ggplot2.discrete.fill = color_palette)
  
  cat("Global plotting theme and color palette set\n")
}

# Save themes and palettes for use in other scripts
export_theme_settings <- function() {
  list(
    theme_allergy = theme_allergy,
    theme_pub = theme_pub,
    theme_trae_publication = theme_trae_publication,
    color_palette = color_palette,
    color_gradient = color_gradient,
    heatmap_palette = heatmap_palette,
    cell_type_palette = cell_type_palette,
    cohort_palette = cohort_palette,
    trae_colors_main = trae_colors_main,
    trae_colors_celltype = trae_colors_celltype,
    trae_colors_pathway = trae_colors_pathway,
    trae_colors_diverging = trae_colors_diverging,
    trae_colors_pvalue = trae_colors_pvalue,
    trae_colors_tissue = trae_colors_tissue,
    trae_colors_cohort = trae_colors_cohort,
    trae_colors_group = trae_colors_group,
    save_trae_figure = save_trae_figure,
    preview_colors = preview_colors,
    set_plot_theme = set_plot_theme
  )
}


# ==============================================================================
# Initialization Console Message
# ==============================================================================
if (interactive()) {
  cat("\n=== TRAE Validation Stance Theme Loaded (English Version) ===\n")
  cat("[+] Enforced pure English comments to prevent encoding errors.\n")
  cat("[+] 11-Cohort Palette mapped: 1 Discovery + 3 Anchors (Vibrant) + 7 Replication (Morandi).\n")
  cat("[+] High-contrast tissue mapping and standard publication theme configured.\n")
  cat("=============================================================\n\n")
}


# Example usage
# source("R/00_theme_setup.R")
# set_plot_theme()
# ggplot(data, aes(x, y, color = group)) + geom_point() + theme_allergy()
# save_trae_figure(plot, "output/figures_main/figure1.pdf")
