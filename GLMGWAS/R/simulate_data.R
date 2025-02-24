#' @title GLM-based GWAS Analysis
#' @description Performs genome-wide association studies (GWAS) using a generalized linear model (GLM).
#' @param y Phenotype data (n x 1).
#' @param X Genotype data (n x m).
#' @param C Covariate data (n x t).
#' @return A numeric vector of p-values (1 x m).
#' @export
simulate_data <- function(n, m, causal_snp = 10) {
  X <- matrix(sample(0:2, n * m, replace = TRUE), n, m)
  beta <- rep(0, m)
  beta[sample(1:m, causal_snp)] <- rnorm(causal_snp)
  y <- X %*% beta + rnorm(n)
  return(list(y = as.data.frame(y), X = X))
}