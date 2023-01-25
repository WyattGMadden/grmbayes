#' Fit the Spatial Ensemble model
#'
#' This function fits the spatial ensemble, using output from the grm_cv function applied to multiple data sets
#'
#' @inheritParams grm
#
#' @param d1 Vector of densities from first dataset grm_cv() output
#' @param d2 Vector of densities from second dataset grm_cv() output
#'
#' @return A data frame containing cross validation predictions
#'
#' @examples
#' # ensemble_spatial()
#' 
#' 
#' @export


ensemble_spatial = function(d1, 
                            d2, 
                            dist.space.mat, 
                            space.id, 
                            n.iter = 25000, 
                            burn = 5000, 
                            thin = 4,
                            tau.a = 0.001, 
                            tau.b = 0.001, 
                            theta.tune = 0.2, 
                            theta.a = 5, 
                            theta.b = 0.05) {
  
    S = max(space.id)
    n = length(d1)
  
    ##Results to save
    #Total number of samples in the end
    K = (n.iter - burn) / thin
    q.save = matrix(NA, ncol = S, nrow = K)
    theta.save = rep(NA, K)
    tau2.save = rep(NA, K)
    dev.save = rep(NA, K)
    theta.acc = 0
  
    ##Initial values
    #w = rep (0.5, S)
    q = rep(0, S)
    z = stats::rbinom(n, 1, 0.5)
    tau2 = 0.01
    theta = 100
    Sigma = tau2 * exp(-1 / theta * dist.space.mat)

    ##Update W
    for (i in 1:n.iter){

        if ((i %% 5000) == 0) {
            cat(paste("  Iteration", i, "of", n.iter, "\n"))
        }
 
        ##Update log-odds q
        m = q[space.id]
        omega = BayesLogit::rpg.devroye(n, 1, m)
    
        O = diag(tapply(omega, space.id, sum))
        VVV = solve(O + solve(Sigma))
        MMM = VVV %*% (tapply(z - 0.5, space.id, sum))  
        q = MASS::mvrnorm(1, MMM, VVV)
    
        #update tau
        SSS = t(q) %*% solve(exp(-dist.space.mat / theta)) %*% q
        tau2 = 1 / stats::rgamma(1, S / 2 + tau.a, SSS / 2 + tau.b)
        Sigma = tau2 * exp(-1 / theta * dist.space.mat)
    
        #Update theta
        theta.prop = stats::rlnorm(1, log(theta), theta.tune)
        SSS.curr = Sigma
        SSS.prop = tau2 * exp(-dist.space.mat / theta.prop)
    
        lik.curr = mvtnorm::dmvnorm(q, rep(0, S), SSS.curr, log = T)
        lik.prop = mvtnorm::dmvnorm(q, rep(0, S), SSS.prop, log = T)
    
        ratio = lik.prop + 
            stats::dgamma(theta.prop, theta.a, theta.b, log = T) + 
            log(theta.prop) -
            lik.curr - 
            stats::dgamma(theta, theta.a, theta.b, log = T) - 
            log(theta)
    
        if (log(stats::runif(1)) < ratio) {
            theta = theta.prop
            theta.acc = theta.acc + 1
        }
    
        ##Update indicator z
        w = 1 / (1 + exp(-q))[space.id]
        p = w * d1 / (w * d1 + (1 - w) * d2)
        z = stats::rbinom(n, 1, p)
    
        if (i > burn & i %% thin == 0) {
            k = (i - burn) / thin
            q.save[k, ] = q
            tau2.save[k] = tau2
            theta.save[k] = theta
            dev.save[k] = -2 * sum(log(w * d1 + (1 - w) * d2))
        }
    }
  
    q.save = data.frame(space.id = 1:S, t(q.save))
    names(q.save) = c("space.id", paste0("Sample",1:K))
  
    other.save = data.frame(tau2 = tau2.save, 
                            theta = theta.save, 
                            dev = dev.save)
  
    list(q = q.save, 
         other = other.save, 
         theta.acc = theta.acc / n.iter)

}

