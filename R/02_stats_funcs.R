#!/usr/bin/env Rscript

# Script 02: Statistical function encapsulation
# Purpose: Cohen's d / Wilcoxon calculation encapsulation, AUC calculation module physically isolated

# Strict policy enforcement: Intercept pROC package loading
if("pROC" %in% loadedNamespaces()) {
  warning("STRICT POLICY: pROC package is loaded. Ensure AUC is NOT calculated for replication-only cohorts.")
}

# Load necessary packages
library(tidyverse)

# Calculate Cohen's d
calculate_cohens_d <- function(group1, group2) {
  # Calculate means
  mean1 <- mean(group1, na.rm = TRUE)
  mean2 <- mean(group2, na.rm = TRUE)
  
  # Calculate standard deviations
  sd1 <- sd(group1, na.rm = TRUE)
  sd2 <- sd(group2, na.rm = TRUE)
  
  # Calculate sample sizes
  n1 <- sum(!is.na(group1))
  n2 <- sum(!is.na(group2))
  
  # Calculate pooled standard deviation
  pooled_sd <- sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))
  
  # Calculate Cohen's d
  cohens_d <- (mean1 - mean2) / pooled_sd
  
  return(cohens_d)
}

# Wilcoxon rank-sum test
calculate_wilcoxon <- function(group1, group2) {
  # Perform Wilcoxon rank-sum test
  test_result <- wilcox.test(group1, group2, paired = FALSE, alternative = "two.sided")
  
  # Extract p-value
  p_value <- test_result$p.value
  
  return(p_value)
}

# Paired Wilcoxon test
calculate_wilcoxon_paired <- function(group1, group2) {
  # Perform paired Wilcoxon test
  test_result <- wilcox.test(group1, group2, paired = TRUE, alternative = "two.sided")
  
  # Extract p-value
  p_value <- test_result$p.value
  
  return(p_value)
}

# Batch calculate statistics between two groups
batch_statistics <- function(data, group_col, value_col) {
  # Ensure data is a data frame
  data <- as.data.frame(data)
  
  # Get two groups
  groups <- unique(data[[group_col]])
  if (length(groups) != 2) {
    stop("Data must contain exactly two groups")
  }
  
  group1 <- data[data[[group_col]] == groups[1], ][[value_col]]
  group2 <- data[data[[group_col]] == groups[2], ][[value_col]]
  
  # Calculate statistics
  cohens_d <- calculate_cohens_d(group1, group2)
  wilcox_p <- calculate_wilcoxon(group1, group2)
  
  # Calculate means and standard deviations
  mean1 <- mean(group1, na.rm = TRUE)
  mean2 <- mean(group2, na.rm = TRUE)
  sd1 <- sd(group1, na.rm = TRUE)
  sd2 <- sd(group2, na.rm = TRUE)
  
  # Return results
  result <- list(
    group1 = groups[1],
    group2 = groups[2],
    mean1 = mean1,
    mean2 = mean2,
    sd1 = sd1,
    sd2 = sd2,
    cohens_d = cohens_d,
    wilcox_p = wilcox_p
  )
  
  return(result)
}

# Note: AUC calculation module is physically isolated in R/04_auc_anchor_only.R
# Only use in cohorts with independent clinical ground truth
# Strictly prohibited to use AUC calculation in non-anchor cohorts

# Example usage
# source("R/02_stats_funcs.R")
# group1 <- rnorm(50, mean = 0, sd = 1)
# group2 <- rnorm(50, mean = 1, sd = 1)
# cohens_d <- calculate_cohens_d(group1, group2)
# wilcox_p <- calculate_wilcoxon(group1, group2)

# For AUC calculation, use R/04_auc_anchor_only.R
# Only use in anchor cohorts with independent clinical truth values
