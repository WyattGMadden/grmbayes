
#standardize covariate matrices in beginning of grm function
covariate_matrix_standardize <- function(x) {
  x <- as.matrix(x)
  x.mean <- apply(x, 2, mean)
  x.sd <- apply(x, 2, stats::sd)
  x[, x.sd != 0] <- sweep(sweep(x[, x.sd != 0], 
                                2, 
                                x.mean[x.sd != 0], 
                                "-"),
                          2,
                          x.sd[x.sd != 0],
                          "/") 
  return(list(x = x, x.mean = x.mean, x.sd = x.sd))
}
