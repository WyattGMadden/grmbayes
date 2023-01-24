#' Fit the Spatial Ensemble model
#'
#' This function fits the spatial ensemble, using output from the grm_cv function applied to multiple data sets
#'
#' @inheritParams grm
#' @param d1 Vector of densities from first dataset grm_cv() output
#' @param d2 Vector of densities from second dataset grm_cv() output
#' @param X.pred Standardized primary independent variable vector for all predictions (n_pred)
#' @param L.pred Standardized spatial covariate matrix for all predictions (n_pred, p1)
#' @param M.pred Standardized spatio-temporal covariate matrix for all predictions (n_pred, p2)
#' @param locations.Y Matrix of primary independent variable x, y coodinates (n_obs, 2)
#' @param locations.pred Matrix of prediction x, y coodinates (n_pred, 2)
#' @param n.iter Number of iterations used in predictions. Must be <= than post-thined and burned iterations from grm.fit
#'
#' @return A data frame containing cross validation predictions
#'
#' @examples
#' # grm_pred()
#' 
#' 
#' @export

library (mvtnorm)
library (BayesLogit)
library (MASS)

ensemble_spatial = function(d1, 
                            d2, 
                            dist.space.mat, 
                            space.id, 
                            n.iter = 25000, 
                            burn = 5000, 
                            thin = 4,
                            a1 = 0.001, 
                            a2 = 0.001, 
                            theta.tune = 0.2, 
                            theta.a = 5, 
                            theta.b = 0.05) {
  
  S = max(space.id)
  n = length (d1)
  
  ##Results to save
  #Total number of samples in the end
  K = (n.iter - burn)/thin
  q.save = matrix (NA, ncol = S, nrow = K)
  theta.save = rep (NA, K)
  tau2.save = rep (NA, K)
  dev.save = rep (NA, K)
  theta.acc = 0
  
  ##Initial values
  #w = rep (0.5, S)
  q = rep(0,S)
  z = rbinom (n, 1, 0.5)
  tau2 = 0.01
  theta = 100
  Sigma = tau2*exp(-1/theta*dist.space.mat)

  ##Update W
  for (i in 1:n.iter){
    
    if ( (i %% 5000) == 0  ){ print (paste("     Iteration", i, "of", n.iter)) }
 
    ##Update log-odds q
    m = q[space.id]
    omega = rpg.devroye(n, 1, m)
    
    O = diag (tapply (omega, space.id, sum))
    VVV = solve ( O+solve(Sigma))
    MMM = VVV %*% (tapply(z-0.5, space.id, sum))  
    q = mvrnorm (1, MMM, VVV)  
    
    #update tau
    SSS = t(q)%*% solve(exp(-dist.space.mat/theta))%*%q
    tau2 = 1/rgamma (1, S/2+a1, SSS/2+a1 )
    Sigma = tau2*exp(-1/theta*dist.space.mat)
    
    #Update theta
    theta.prop = rlnorm(1, log(theta), theta.tune)
    SSS.curr = Sigma
    SSS.prop = tau2*exp(-dist.space.mat/theta.prop)
    
    lik.curr =  dmvnorm ( q, rep(0,S), SSS.curr, log = T)
    lik.prop =  dmvnorm ( q, rep(0,S), SSS.prop, log = T)
    
    ratio = lik.prop + dgamma (theta.prop, theta.a, theta.b, log = T) + log (theta.prop) -
      lik.curr - dgamma (theta, theta.a, theta.b, log = T) - log(theta)
    
    if ( log(runif(1)) < ratio){
      theta = theta.prop;
      theta.acc = theta.acc + 1
    }
    
    ##Update indicator z
    w = 1/(1+exp(-q))[space.id]
    p = w*d1/(w*d1 + (1-w)*d2)
    z = rbinom (n, 1, p)
    
    if (i > burn & i %%thin == 0){
      k = (i - burn)/thin
      q.save[k,] = q
      tau2.save[k] = tau2
      theta.save[k] = theta
      dev.save[k] = -2*sum(log(w*d1+(1-w)*d2))
    }
    
  }
  
  q.save = data.frame (space.id = 1:S, t(q.save))
  names (q.save) = c("space.id", paste0("Sample",1:K))
  
  other.save = data.frame (tau2 = tau2.save, theta = theta.save, dev = dev.save)
  
  list (q = q.save, other= other.save, theta.acc = theta.acc/n.iter)

}

