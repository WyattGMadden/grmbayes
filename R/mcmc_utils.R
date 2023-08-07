

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


# Function that calculates the jump for discrete theta metropolis-hastings and corresponding probability
# returns either -1 or 1, and accounts for edge cases that have double probability
get_discrete_theta_mh_jump <- function(which_theta_curr, discrete.theta.values) {
    if (which_theta_curr == 1) { 
        jump <- 1
        lik_jump_curr_to_prop <- 1
        lik_jump_prop_to_curr <- 0.5
    } else if (which_theta_curr == length(discrete.theta.values)) {
        jump <- -1
        lik_jump_curr_to_prop <- 1
        lik_jump_prop_to_curr <- 0.6
    } else {
        jump <- sample(c(-1, 1), 1)
        lik_jump_curr_to_prop <- 0.5
        if ((jump + which_theta_curr) %in% c(1, length(discrete.theta.values))) {

            lik_jump_prop_to_curr <- 1
        } else {

            lik_jump_prop_to_curr <- 0.5
        }
    }
    return(list(jump = jump, 
                lik_jump_curr_to_prop = lik_jump_curr_to_prop, 
                lik_jump_prop_to_curr = lik_jump_prop_to_curr))
}
