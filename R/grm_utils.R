
#standardize covariate matrices in beginning of grm function
covariate_matrix_standardize <- function(x) {
  x <- as.matrix(x)
  x.mean <- apply(x, 2, mean)
  x.sd <- apply(x, 2, stats::sd)
  x[, x.sd != 0] <- sweep(sweep(x[, x.sd != 0, drop = F], 
                                2, 
                                x.mean[x.sd != 0, drop = F], 
                                "-"),
                          2,
                          x.sd[x.sd != 0],
                          "/") 
  return(list(x = x, x.mean = x.mean, x.sd = x.sd))
}

#density of inverse gamma
dinvgamma <- function(x, alpha, beta, log = F) {
    if (log) {
        return(stats::dgamma(1 / x, alpha, beta, log = T) - 2 * log(x))
    } else {
        return(stats::dgamma(1 / x, alpha, beta) * (1 / x)^2)
    }
}
