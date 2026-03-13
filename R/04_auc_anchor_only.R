# ============================================================================
# Script 04: AUC calculation functions - Anchor Only
# ============================================================================
# Purpose: Contains all functions that depend on pROC or ROCR for AUC calculation and ROC curve plotting
# WARNING: These functions should only be used with cohorts that have independent clinical truth values
# ============================================================================

# Load required packages
library(pROC)
library(ROCR)

# ============================================================================
# AUC calculation functions
# ============================================================================

# Calculate AUC for anchor cohorts
# Calculate AUC for anchor cohorts with independent clinical truth values
# 
# Parameters:
# predictions: Numeric vector of predicted scores
# labels: Binary vector of true labels (0/1 or factor)
# 
# Returns:
# AUC value
calculate_auc <- function(predictions, labels) {
  if (!requireNamespace("pROC", quietly = TRUE)) {
    stop("pROC package is required for AUC calculation")
  }
  
  roc_obj <- pROC::roc(response = labels, predictor = predictions)
  return(pROC::auc(roc_obj))
}

# Plot ROC curve for anchor cohorts
# Plot ROC curve for anchor cohorts with independent clinical truth values
# 
# Parameters:
# predictions: Numeric vector of predicted scores
# labels: Binary vector of true labels (0/1 or factor)
# title: Plot title
# 
# Returns:
# ggplot object
plot_roc_curve <- function(predictions, labels, title = "ROC Curve") {
  if (!requireNamespace("pROC", quietly = TRUE)) {
    stop("pROC package is required for ROC curve plotting")
  }
  
  roc_obj <- pROC::roc(response = labels, predictor = predictions)
  
  p <- pROC::ggroc(roc_obj, color = "#B31B1B", size = 1.5) +
    ggplot2::geom_abline(intercept = 1, slope = 1, linetype = "dashed", color = "gray") +
    ggplot2::labs(title = title,
                  x = "1 - Specificity",
                  y = "Sensitivity") +
    ggplot2::theme_minimal() +
    ggplot2::annotate("text", x = 0.75, y = 0.25,
                      label = paste("AUC =", round(pROC::auc(roc_obj), 3)),
                      size = 5, color = "#B31B1B")
  
  return(p)
}

# Calculate sensitivity and specificity at specific threshold
# Calculate sensitivity and specificity at a specific threshold
# 
# Parameters:
# predictions: Numeric vector of predicted scores
# labels: Binary vector of true labels (0/1 or factor)
# threshold: Threshold value
# 
# Returns:
# List with sensitivity and specificity
calculate_sensitivity_specificity <- function(predictions, labels, threshold) {
  if (!requireNamespace("pROC", quietly = TRUE)) {
    stop("pROC package is required for sensitivity/specificity calculation")
  }
  
  roc_obj <- pROC::roc(response = labels, predictor = predictions)
  coords <- pROC::coords(roc_obj, x = threshold, input = "threshold",
                         ret = c("sensitivity", "specificity"))
  
  return(list(
    sensitivity = coords[["sensitivity"]],
    specificity = coords[["specificity"]]
  ))
}

# ============================================================================
# ROCR-based functions (alternative implementation)
# ============================================================================

# Calculate AUC using ROCR
# Calculate AUC using ROCR package
# 
# Parameters:
# predictions: Numeric vector of predicted scores
# labels: Binary vector of true labels (0/1)
# 
# Returns:
# AUC value
calculate_auc_rocr <- function(predictions, labels) {
  if (!requireNamespace("ROCR", quietly = TRUE)) {
    stop("ROCR package is required for this function")
  }
  
  pred <- ROCR::prediction(predictions, labels)
  perf <- ROCR::performance(pred, "auc")
  return(as.numeric(perf@y.values[[1]]))
}

# Plot ROC curve using ROCR
# Plot ROC curve using ROCR package
# 
# Parameters:
# predictions: Numeric vector of predicted scores
# labels: Binary vector of true labels (0/1)
# title: Plot title
# 
# Returns:
# Plot object
plot_roc_curve_rocr <- function(predictions, labels, title = "ROC Curve") {
  if (!requireNamespace("ROCR", quietly = TRUE)) {
    stop("ROCR package is required for this function")
  }
  
  pred <- ROCR::prediction(predictions, labels)
  perf <- ROCR::performance(pred, "tpr", "fpr")
  
  plot(perf, main = title, col = "#B31B1B", lwd = 2)
  abline(a = 0, b = 1, lty = 2, col = "gray")
  
  # Add AUC to plot
  auc <- calculate_auc_rocr(predictions, labels)
  text(0.7, 0.2, paste("AUC =", round(auc, 3)), col = "#B31B1B", cex = 1.2)
  
  return(perf)
}

# ============================================================================
# Strict policy enforcement
# ============================================================================

# Check if AUC calculation is allowed
# Check if AUC calculation is allowed for the given cohort type
# 
# Parameters:
# cohort_type: Character string indicating cohort type ("anchor" or "replication")
# 
# Returns:
# TRUE if AUC is allowed, FALSE otherwise
check_auc_allowed <- function(cohort_type) {
  if (tolower(cohort_type) != "anchor") {
    warning("STRICT POLICY: AUC calculation is only allowed for anchor cohorts with independent clinical truth values")
    return(FALSE)
  }
  return(TRUE)
}
