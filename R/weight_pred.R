
#' Generate weights from ensemble spatial fit
#'
#' This function generates weights using the output from the ensemble_spatial function
#'
#' @inheritParams grm_pred
#' @param q Matrix of samples output from ensemble_spatial function
#' @param theta Vector of theta samples output from ensemble_spatial function
#' @param tau2 Vector of tau2 samples output from ensemble_spatial function
#'
#' @return A matrix containing weights
#'
#' @examples
#' # grm_pred()
#' 
#' 
#' @export

weight_pred = function(q, 
                       theta, 
                       tau2, 
                       locations.Y, 
                       locations.pred, 
                       verbose = TRUE) {

    n.iter = length(theta)
  
  
    #Some data processing
    q.mat = as.matrix(q[, -1])
  
    N.mon = nrow(locations.Y)
    N.cell = nrow(locations.pred)
  
    XY = rbind(locations.Y, 
               locations.pred)[, -1]
    D22 = as.matrix(stats::dist(locations.Y, 
                                diag = TRUE, 
                                upper = TRUE))
    D12 = as.matrix(stats::dist(XY, 
                                diag = TRUE, 
                                upper = TRUE))[c(1:N.mon), -c(1:N.mon)]
  
    q.pred = matrix(NA, N.cell, n.iter)
  
    if (verbose == TRUE) {
        print("Imputing Spatial Weights") 
    }
  
    for (m in 1:n.iter) {
   
        if (verbose == TRUE & m%%100 == 0) {
            print(paste0("Imputing ", m, " of ", n.iter))
        } 
    
        q.m = q.mat[, m]
        tau2.m = tau2[m]
        theta.m = theta[m]
        Sigma11.m = tau2.m
        Sigma12.m = tau2.m * exp(-1 / theta.m * D12)
        Sigma22.m = tau2.m * exp(-1 / theta.m * D22)
        InvSigma22.m = solve(Sigma22.m)
    
        q.mu.m = t(Sigma12.m) %*% InvSigma22.m %*% q.m
        q.cov.m = Sigma11.m - diag(t(Sigma12.m) 
                                   %*% InvSigma22.m 
                                   %*% Sigma12.m)
        q.m.post = stats::rnorm(N.cell, q.mu.m, sqrt(q.cov.m))
    
        q.pred[, m] = q.m.post
    }
  
    return(q.pred)
}
  
