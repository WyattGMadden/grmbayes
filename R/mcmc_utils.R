

# Function to compute the log-density of a multivariate normal distribution
# using univariate density and lower diagonal cholesky decomposition of covariance.
# Sum of univariate normal densities is multiplied by jacobian (1/2) * determinant(Sigma)
d_mvn_chol_uvn <- function(x, chol_inv, det_sigma) {
    x_transformed <- chol_inv %*% x
    univ_density <- sum(stats::dnorm(x_transformed, log = T))
    jacobian <- 0.5 * det_sigma
    full_density <- univ_density - jacobian
    return(full_density)
}
