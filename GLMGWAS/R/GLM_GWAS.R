#' @title GLM-based GWAS Analysis
#' @description Performs genome-wide association studies (GWAS) using a generalized linear model (GLM).
#' @param y Phenotype data (n x 1).
#' @param X Genotype data (n x m).
#' @param C Covariate data (n x t).
#' @return A numeric vector of p-values (1 x m).
#' @export
GLM_GWAS <- function(y, X, C) {
  if (!is.data.frame(y) | !is.matrix(X) | !is.data.frame(C)) {
    stop("y must be a data frame, X must be a matrix, and C must be a data frame.")
  }
  
  n <- nrow(X)
  if (nrow(y) != n | nrow(C) != n) stop("Mismatch in number of individuals (rows).")
  
  p_values <- numeric(ncol(X))  # Store p-values
  
  for (i in 1:ncol(X)) {
    df <- data.frame(y = y[[1]], SNP = X[, i], C)
    model <- glm(y ~ SNP + ., data = df, family = gaussian())
    p_values[i] <- summary(model)$coefficients["SNP", 4]
  }
  
  names(p_values) <- colnames(X)
  return(p_values)
}
