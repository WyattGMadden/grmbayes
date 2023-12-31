

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

#precalculate lists of kernals for discrete theta nngp 
# used for theta alpha and theta beta in grm
get_discrete_nngp_kernals <- function(discrete.values, dist_matrices, cov_kern) {

    kernals <- lapply(discrete.values, 
                      function(x) {
                          lapply(dist_matrices, 
                                 function(y) {
                                     cov_kern(distance = y, 
                                              theta = x)
                                 })
                      })

    kernals_inv <- lapply(kernals, 
                          function(x) {
                              lapply(x, 
                                     function(y) {
                                         solve(y)
                                     })
                          })

    kernals_partial_inv <- lapply(kernals, 
                                  function(x) {
                                      lapply(x[1:(length(x) - 1)],
                                             function(y) {
                                                 solve(y[-1, -1, drop = FALSE])
                                             })
                                  })
    return(list(kernals = kernals, 
                kernals_inv = kernals_inv, 
                kernals_partial_inv = kernals_partial_inv))
}

#precalculate lists of kernals for discrete theta gp 
# used for theta alpha and theta beta in grm
get_discrete_gp_kernals <- function(discrete.values, dist.space.mat, cov_kern) {

    kernals <- lapply(discrete.values,
                           function(x) cov_kern(distance = dist.space.mat, 
                                                theta = x))
    kernals_inv <- lapply(kernals, solve)
    kernals_chol <- lapply(kernals, 
                           function(x) t(chol(x)))
    kernals_chol_inv <- lapply(kernals_chol, 
                               function(x) solve(x))
    kernals_det <- lapply(kernals, 
                          function(x) as.numeric(determinant(x)$modulus))

    return(list(kernals = kernals, 
                kernals_inv = kernals_inv, 
                kernals_chol = kernals_chol,
                kernals_chol_inv = kernals_chol_inv,
                kernals_det = kernals_det))
}

init_discrete_theta <- function(discrete.theta.values) {

    which_theta_curr <- ceiling(length(discrete.theta.values) / 2)
    theta = discrete.theta.values[which_theta_curr]
    return(list(which_theta_curr = which_theta_curr, 
                theta = theta))
}
