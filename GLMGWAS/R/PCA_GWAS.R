#' @title GLM-based GWAS Analysis
#' @description Performs genome-wide association studies (GWAS) using a generalized linear model (GLM).
#' @param y Phenotype data (n x 1).
#' @param X Genotype data (n x m).
#' @param C Covariate data (n x t).
#' @return A numeric vector of p-values (1 x m).
#' @export
PCA_GWAS <- function(X, C, num_PCs = 5) {
  pca_res <- prcomp(X, scale. = TRUE)
  PCs <- as.data.frame(pca_res$x[, 1:num_PCs])
  
  # Check for linear dependence
  C_combined <- cbind(C, PCs)
  lm_test <- lm(as.matrix(PCs) ~ as.matrix(C))
  dependent_PCs <- which(summary(lm_test)$r.squared > 0.95)  # Less strict threshold
  
  if (length(dependent_PCs) > 0) {
    PCs <- PCs[, -dependent_PCs, drop = FALSE]
  }
  
  # Ensure at least one PC is retained
  if (ncol(PCs) == 0) {
    warning("All PCs are collinear with covariates. Retaining first PC to avoid empty matrix.")
    PCs <- as.data.frame(pca_res$x[, 1, drop = FALSE])
  }
  
  return(PCs)
}
