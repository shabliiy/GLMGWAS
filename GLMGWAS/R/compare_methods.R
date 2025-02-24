#' @title GLM-based GWAS Analysis
#' @description Performs genome-wide association studies (GWAS) using a generalized linear model (GLM).
#' @param y Phenotype data (n x 1).
#' @param X Genotype data (n x m).
#' @param C Covariate data (n x t).
#' @return A numeric vector of p-values (1 x m).
#' @export
compare_methods <- function(n = 100, m = 1000, replicates = 30) {
  results <- data.frame(GLM_p = numeric(replicates), GWASbyCor_p = numeric(replicates))
  
  for (i in 1:replicates) {
    data <- simulate_data(n, m)
    y <- data$y
    X <- data$X
    C <- matrix(rnorm(n), n, 1)
    
    p_glm <- GLM_GWAS(y, X, C)
    p_gwasbycor <- cor.test(y[[1]], X[, 1])$p.value  # Placeholder for GWASbyCor
    
    results[i, ] <- c(mean(p_glm), p_gwasbycor)
  }
  
  return(results)
}