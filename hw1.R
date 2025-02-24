# Load required libraries
library(stats)
library(data.table)
library(ggplot2)
library(MASS)
library(Matrix)
library(parallel)
library(qqman)


# Load data files
genotype_data <- fread("mdp_numeric.txt", header = TRUE)
snp_info <- fread("mdp_SNP_information.txt", header = TRUE)
phenotype_data <- fread("CROP545_Phenotype.txt", header = TRUE)
covariates_data <- fread("CROP545_Covariates.txt", header = TRUE)

# Function to perform GLM GWAS
glm_gwas <- function(y, X, C = NULL) {
  n <- length(y)
  m <- ncol(X)
  p_values <- numeric(m)
  
  for (i in 1:m) {
    model <- glm(y ~ X[, i] + C, family = gaussian())
    p_values[i] <- summary(model)$coefficients[2, 4]  # Extract p-value for SNP
  }
  
  return(data.frame(SNP = colnames(X), p_value = p_values))
}

# Function to compute PCA for population structure
compute_pcs <- function(X, num_pcs = 10) {
  pca_res <- prcomp(X, scale. = TRUE)
  pcs <- pca_res$x[, 1:num_pcs]
  return(pcs)
}

# Function to filter linearly dependent PCs
filter_independent_pcs <- function(pcs, C) {
  independent_pcs <- pcs[, !apply(pcs, 2, function(pc) any(abs(cor(pc, C)) > 0.9))]
  return(independent_pcs)
}

# Function for multiple testing correction
multiple_testing_correction <- function(p_values, method = "bonferroni") {
  adj_p <- p.adjust(p_values, method = method)
  return(adj_p)
}

# Function to simulate GWAS data
simulate_gwas_data <- function(n, m, t) {
  y <- rnorm(n)
  X <- matrix(sample(0:2, n * m, replace = TRUE), nrow = n, ncol = m)
  C <- matrix(rnorm(n * t), nrow = n, ncol = t)
  return(list(y = y, X = X, C = C))
}

# Function to compare GWAS results with GWASbyCor
simulate_and_compare <- function(n_reps = 30) {
  pvals_glm <- matrix(NA, n_reps, 1000)
  pvals_gwasbycor <- matrix(NA, n_reps, 1000)
  
  for (i in 1:n_reps) {
    sim_data <- simulate_gwas_data(200, 1000, 3)
    pvals_glm[i, ] <- glm_gwas(sim_data$y, sim_data$X, sim_data$C)$p_value
    pvals_gwasbycor[i, ] <- run_gwasbycor(sim_data$y, sim_data$X)  # Placeholder function
  }
  
  return(list(glm = pvals_glm, gwasbycor = pvals_gwasbycor))
}

# Function to visualize Manhattan plot
plot_manhattan <- function(gwas_results) {
  gwas_results$SNP <- as.factor(gwas_results$SNP)
  gwas_results$logP <- -log10(gwas_results$p_value)
  
  ggplot(gwas_results, aes(x = SNP, y = logP)) +
    geom_point() +
    theme_minimal() +
    labs(title = "Manhattan Plot", x = "SNP", y = "-log10(p-value)")
}

# Example Usage
# results <- glm_gwas(phenotype_data$Phenotype, genotype_data, covariates_data)
# results$p_value_adj <- multiple_testing_correction(results$p_value)
# plot_manhattan(results)

