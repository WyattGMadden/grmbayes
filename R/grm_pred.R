#' Fit the Geostatistical Regression Model (GRM) With Cross Validation
#'
#' This function fits Bayesian Hierarchical Model (BHM) in the form of Y ~ beta X + gamma L + delta M with cross-validation
#'
#' @inheritParams grm
#' @param grm.fit Fit object created with grm()
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

grm_pred = function(grm.fit,
                    X.pred, 
                    L.pred, 
                    M.pred, 
                    locations.Y, 
                    locations.pred, 
                    space.id, 
                    time.id, 
                    spacetime.id, 
                    include.additive.annual.resid = T,
                    include.multiplicative.annual.resid = T,
                    n.iter = 500,
                    verbose = TRUE) {
  
  ###Print some information
  print("Preparing for Prediction")
  
  N = length(X.pred) #Total AOD observations
  N.space = max(space.id) #Total number of prediction cells
  N.time = max(time.id) #Maximum number of time interval (weeks) in prediction
  N.time.obs = length(unique(time.id)) #Number of observed time interval (weeks)
  N.spacetime = max(spacetime.id) #Time points where spatial trends vary by (year)
  
  N.Lmax = ncol(L.pred) #Number of spatial predictors to use
  N.Mmax = ncol(M.pred) #Number of spatial-temporal predictors to use
  
  if (verbose == TRUE) {
      print("#################################### ")
      print("######## Preparing for MCMC ######## ")
      print("#################################### ")
      print(paste0("     ", 
                   names(X.pred)[4], 
                   " GRM Model"))
      print(paste0("     Total number of prediction time points: ", 
                   N.time.obs, 
                   " (out of ", 
                   N.time,")") )
      print(paste0("     Total number of prediction cells: ", 
                   N.space))
      print(paste0("     Total number of spatial covariates: ", 
                   N.Lmax))
      print(paste0("     Total number of spatial-temporal covariates: ", 
                   N.Mmax))
  }
  
  N.mon = nrow(locations.Y)
  N.cell = nrow(locations.pred)
  
  ####Predict alpha and beta at grid cells
  XY = rbind(locations.Y, locations.pred)
  D22 = as.matrix(stats::dist(locations.Y, diag = TRUE, upper = TRUE))
  D12 = as.matrix(stats::dist(XY, diag = TRUE, upper = TRUE))[c(1:N.mon), -c(1:N.mon)]
  
  alpha_space_pred  = data.frame(expand.grid(1:N.space, 
                                             1:N.spacetime))
  beta_space_pred = data.frame(expand.grid(1:N.space, 
                                           1:N.spacetime))

  names(alpha_space_pred) = c("space.id", "spacetime.id")
  names(beta_space_pred) = c("space.id", "spacetime.id")

  alpha_space_pred[paste0("Sample", 1:n.iter)] = 0
  beta_space_pred[paste0("Sample", 1:n.iter)] = 0
  
  #For alpha's
  if (include.additive.annual.resid) {

    if (verbose == TRUE) {

      print("Imputing Spatial Alphas") 

    }

    for (m in 1:n.iter) {
        tau.m = grm.fit$others$tau.alpha[m]
        theta.m = grm.fit$others$theta.alpha[m]
        Sigma11.m = tau.m
        Sigma12.m = tau.m * exp(-1 / theta.m * D12)
        Sigma22.m = tau.m * exp(-1 / theta.m * D22)
        InvSigma22.m = solve(Sigma22.m)
  
        for (j in 1:N.spacetime) {
            alpha.m = grm.fit$alpha.space[grm.fit$alpha.space$spacetime.id == j, 
                                          paste0("Sample", m)]
            alpha.mu.m = t(Sigma12.m) %*% InvSigma22.m %*% alpha.m
            alpha.cov.m = Sigma11.m - diag(t(Sigma12.m) %*% InvSigma22.m %*% Sigma12.m)
            alpha.m.post = stats::rnorm(N.cell, alpha.mu.m, sqrt(alpha.cov.m))
            alpha_space_pred[alpha_space_pred$spacetime.id == j, 
                             paste0("Sample",m)] = alpha.m.post
        } #End of locations
      
        if (verbose == TRUE) {
            if (m %% (n.iter / 10) == 0) {
                print(paste("     Iteration", m, "of", n.iter))
            }
        }
    }
  }
   
  #For betas's
  if (include.multiplicative.annual.resid) {
      if (verbose == TRUE) {
          print("Imputing Spatial Betas") 
      }

      for (m in 1:n.iter) {
          tau.m = grm.fit$others$tau.beta[m]
          theta.m = grm.fit$others$theta.beta[m]
          Sigma11.m = tau.m
          Sigma12.m = tau.m * exp(-1 / theta.m * D12)
          Sigma22.m = tau.m * exp(-1 / theta.m * D22)
          InvSigma22.m = solve(Sigma22.m)
    
          for (j in 1:N.spacetime) {
              beta.m = grm.fit$beta.space[grm.fit$beta.space$spacetime.id == j, 
                                  paste0("Sample", m)]
              beta.mu.m = t(Sigma12.m) %*% InvSigma22.m %*% beta.m
              beta.cov.m = Sigma11.m - diag(t(Sigma12.m) %*% 
                                            InvSigma22.m %*% 
                                            Sigma12.m)
              beta.m.post = stats::rnorm(N.cell, beta.mu.m, sqrt(beta.cov.m))
              beta_space_pred[beta_space_pred$spacetime.id == j, 
                              paste0("Sample", m)] = beta.m.post
          } #End of locations
          if (verbose == TRUE) {
              if (m %% (n.iter / 10) == 0) {
                  print(paste("     Iteration", m, "of", n.iter))
              }
          }
  
      }#End of iterations
  }
  
  
  ####Make Predictions
  results = data.frame(time.id, space.id, spacetime.id)
  results$estimate = 0
  results$sd = 0
  
  id.temp.pred = paste0(alpha_space_pred$space.id, 
                        "_", 
                        alpha_space_pred$spacetime.id)

  id.temp = paste0(space.id, "_", spacetime.id)
  
  delta = as.matrix(grm.fit$delta)
  gamma = as.matrix(grm.fit$gamma)
  
  print("Wrapping up Predictions") 
  
  for (m in 1:n.iter) {
      intercept = grm.fit$others$alpha0[m] + 
          grm.fit$alpha.time[time.id, m + 1] + 
          alpha_space_pred[match(id.temp, id.temp.pred), m + 2]
      slope = grm.fit$others$beta0[m] + 
          grm.fit$beta.time[time.id, m + 1] + 
          beta_space_pred[match(id.temp, id.temp.pred), m + 2]
      fix.L = L.pred %*% gamma[m,]
      fix.M = M.pred %*% delta[m,]
    
      pred.mu = intercept + slope * X.pred + fix.L + fix.M 
      pred.mu = pred.mu + stats::rnorm(length(pred.mu), 0, sqrt(grm.fit$others$sigma2[m])) 
      results$estimate = results$estimate + pred.mu / n.iter
      results$sd = results$sd + pred.mu^2 / n.iter
  }

  results$sd = sqrt((results$sd - results$estimate^2))


  return(results)

}
