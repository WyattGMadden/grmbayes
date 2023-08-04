
#' Fit the Geostatistical Regression Model (GRM)
#'
#' This function fits Bayesian Hierarchical Model (BHM) in the form of Y ~ beta X + gamma L + delta M
#'
#' @param Y Dependent variable vector (n)
#' @param X Unstandardized primary independent variable vector (n)
#' @param L Unstandardized spatial covariate matrix (n, p1)
#' @param M Unstandardized spatio-temporal covariate matrix (n, p2)
#' @param coords Matrix of x y coordinates, with colnames(coords) == c("x", "y"), (n, 2)
#' @param space.id Spatial location ID vector (n)
#' @param time.id Temporal location ID vector (n)
#' @param spacetime.id ID vector of time points where spatial trends vary (n)
#' @param include.additive.weekly.resid Include additive temporal (weekly) residual bias not explained by other inputes
#' @param include.additive.annual.resid Include additive spatial (annual) residual bias not explained by other inputes
#' @param include.multiplicative.weekly.resid Include multiplicative temporal (weekly) residual bias not explained by other inputes
#' @param include.multiplicative.annual.resid Include multiplicative spatial (weekly) residual bias not explained by other inputes
#' @param nngp Use nearest neighbor Gaussian process (NNGP) in place of Gaussian process
#' @param num_neighbors Number of nearest neighbors to use in NNGP
#' @param discrete.theta.alpha.values Values of theta (GP intercept range parameter) to use in discrete uniform prior. If NULL (default) continuous theta prior is used
#' @param discrete.theta.beta.values Values of theta (GP slope range parameter) to use in discrete uniform prior. If NULL (default) continuous theta prior is used
#' @param discrete.theta.gibbs If TRUE, use Gibbs sampling to sample theta.alpha and theta.beta. If FALSE, use Metropolis-Hastings (jump's to nearest discrete values)
#' @param n.iter Number of iterations used in predictions. 
#' @param burn Number of pre-covergence simulations
#' @param thin Save every thin'th simulation
#' @param covariance Specify covariance function (from "exponential", "matern", "custom")
#' @param covariance.kernal Specify a custom covariance function if covariance = "custom". Must be a function with "distance" and "theta" parameters.
#' @param matern.nu Specify nu parameter for Matern covariance function if used (from 0.5, 1.5, and 2.5)
#' @param tau.alpha.tune Tau alpha Metropolis-Hastings proposal tuning parameter, only used if nngp = T
#' @param tau.alpha.a First tau alpha prior hyperparameter
#' @param tau.alpha.b Second tau alpha prior hyperparameter
#' @param tau.beta.tune Tau beta Metropolis-Hastings proposal tuning parameter, only used if nngp = T
#' @param tau.beta.a First tau beta prior hyperparameter
#' @param tau.beta.b Second tau beta prior hyperparameter
#' @param omega.a First omega prior hyperparameter
#' @param omega.b Second omega prior hyperparameter
#' @param theta.alpha.tune Theta alpha Metropolis-Hastings proposal tuning parameter
#' @param theta.alpha.a First theta alpha prior hyperparameter
#' @param theta.alpha.b Second theta alpha prior hyperparameter
#' @param theta.alpha.init Initial value for theta alpha
#' @param theta.beta.tune Theta beta Metropolis-Hastings proposal tuning parameter
#' @param theta.beta.a First theta beta prior hyperparameter
#' @param theta.beta.b Second theta beta prior hyperparameter
#' @param theta.beta.init Initial value for theta beta
#' @param rho.alpha.init Initial value for rho alpha
#' @param rho.beta.init Initial value for rho beta
#' @param sigma.a First sigma prior hyperparameter
#' @param sigma.b Second sigma prior hyperparameter
#' @param verbose Print MCMC output
#' @param verbose.iter print MCMC output step number each 'verbose.iter' iterations
#'
#' @return A list containing MCMC output 
#'
#' @examples
#' # grm()
#' 
#' 
#' @export
grm = function(Y, 
               X, 
               L = NULL, 
               M = NULL, 
               coords, 
               space.id, 
               time.id, 
               spacetime.id,
               include.additive.weekly.resid = T,
               include.additive.annual.resid = T,
               include.multiplicative.weekly.resid = T,
               include.multiplicative.annual.resid = T,
               nngp = F,
               num_neighbors = 10,
               discrete.theta.alpha.values = NULL,
               discrete.theta.beta.values = NULL,
               discrete.theta.gibbs = T,
               n.iter = 25000,
               burn = 5000,
               thin = 4,
               covariance = "exponential",
               covariance.kernal = NULL,
               matern.nu = 1.5,
               tau.alpha.tune = 0.2, 
               tau.alpha.a = 0.5,
               tau.alpha.b = 0.005,
               tau.beta.tune = 0.2, 
               tau.beta.a = 0.5,
               tau.beta.b = 0.005,
               omega.a = 0.5,
               omega.b = 0.005,
               theta.alpha.tune = 0.2, 
               theta.alpha.a = 5, 
               theta.alpha.b = 0.05,
               theta.alpha.init = 100,
               theta.beta.tune = 0.2, 
               theta.beta.a = 5, 
               theta.beta.b = 0.05,
               theta.beta.init = 100,
               rho.alpha.init = 0.9999,
               rho.beta.init = 0.9999,
               sigma.a = 0.001, 
               sigma.b = 0.001,
               verbose = TRUE,
               verbose.iter = 1000) {
 
    ####################################
    ### Get some summary statistics ####
    ####################################
  
    N = length(Y) #Total AOD-PM linked pairs
    N.space = max(space.id) #Total number of monitors in training
    N.time = max(time.id) #Maximum number of time interval (weeks) in training
    N.time.obs = length(unique(time.id)) #Number of observed time interval (weeks)
    N.spacetime = max(spacetime.id) #Time points where spatial trends vary by (year)
    
    Space.labels = 1:N.space 
    Day.labels = 1:N.time
    
    N.Lmax = ncol(L) #Number of spatial predictors to use
    N.Mmax = ncol(M) #Number of spatial-temporal predictors to use
    
    if (verbose == TRUE) {
      print("#################################### ")
      print("######## Preparing for MCMC ######## ")
      print("#################################### ")
      print(paste0("     Total number of available days: ", 
                    N.time.obs, 
                    " (out of ", N.time,")")
      )
      print(paste0("     Total number of monitors: ", 
                   N.space)
      )
      print(paste0("     Total number of spatial covariates: ", 
                   ncol(L))
      )
      print(paste0("     Total number of spatial-temporal covariates: ", 
                   ncol(M))
      )
    }

    #################################
    ### Specify Covariance Kernal ###
    #################################

    if (covariance == "exponential") {

      cov_kern <- exponential_kernal

    } else if (covariance == "matern") {

      cov_kern <- function(distance, theta) {
          matern_kernal(distance = distance, 
                        theta = theta, 
                        nu = matern.nu)
      }
    } else if (covariance == "custom") {

      cov_kern <- covariance.kernal
      
    } else {

      stop("Covariance function not recognized")

    }

    ##############################
    ### Standardize X, L and M ###
    ##############################
    L = as.matrix(L)
    M = as.matrix(M)
    
    ### Calculate means and standard deviation
    X.mean = mean(X)
    if (!is.null(M)) M.mean = apply(M, 2, mean)
    if (!is.null(L)) L.mean = apply(L, 2, mean)
    
    X.sd = stats::sd(X)
    if (!is.null(M)) M.sd = apply(M, 2, stats::sd)
    if (!is.null(L)) L.sd = apply(L, 2, stats::sd)
    
    ### Standardize
    X = (X - mean(X)) / X.sd
    if (!is.null(M) & !any(M.sd == 0)) { 
        M = sweep(sweep(M, 2, M.mean, "-"), 2, M.sd, "/") 
    }

    if (!is.null(L) & !any(L.sd == 0)) {
        L = sweep(sweep(L, 2, L.mean, "-"), 2, L.sd, "/") 
    }

    ######################################
    ### Create Spatial Distance Matrix ###
    ######################################

    if (!nngp) {
        dist.space.mat <- unique(cbind(space.id, coords))
        dist.space.mat <- dist.space.mat[order(dist.space.mat$space.id), ]
        dist.space.mat <- as.matrix(stats::dist(dist.space.mat[, c("x", "y")], 
                                                diag = TRUE, 
                                                upper = TRUE))
    }

    
    #######################################
    ### Create Temporal Distance Matrix ###
    #######################################
    
    ### Create full temporal distance matrix
    A = spam::as.spam(matrix(0, ncol=N.time, N.time)) #Adjacency matrix
    for (i in 1:(N.time - 1)) {
        A[i, i + 1] = 1
    } 

    for (i in 2:(N.time)) { 
        A[i, i - 1] = 1 
    }

    Time.obs = 1:N.time %in% time.id ## T/F for which days are observed
    
    Q = spam::diag.spam(apply(A, 2, sum)) #Row sums of A 

    ##########################################################################
    ### Calculate Nearest Neighbor Gaussian Process Neighbors and Ordering ###
    ##########################################################################
    
    ### get coordinate ordering (from upper right to lower left of spatial grid)
    if (nngp) {
        unique_coords <- unique(cbind(space.id, coords))
        unique_coords <- unique_coords[order(unique_coords$space.id), ]
        coord_list <- order_coords(coords = unique_coords[, c("x", "y")],
                                   space_id = unique_coords$space.id)
        ordered_coords <- coord_list[["ordered_coords"]]
        coord_ordering <- coord_list[["coord_ordering"]]
        coord_reverse_ordering <- coord_list[["coord_reverse_ordering"]]

        ### get list of neighbors & dist matrices w.r.t. coord ordering
        neighbors <- get_neighbors(ordered_coords, m = num_neighbors)
        dist_matrices <- get_dist_matrices(ordered_coords, neighbors)

    }

    ##################################################################
    ### Pre-compute quantities for discrete theta (gp range parameter)
    ##################################################################

    if (!is.null(discrete.theta.alpha.values)) {

        kernals_alpha <- lapply(discrete.theta.alpha.values,
                               function(x) cov_kern(distance = dist.space.mat, 
                                                    theta = x))
        kernals_inv_alpha <- lapply(kernals_alpha, solve)
        kernals_chol_alpha <- lapply(kernals_alpha, 
                               function(x) t(chol(x)))
        kernals_chol_inv_alpha <- lapply(kernals_chol_alpha, 
                                   function(x) solve(x))
        kernals_det_alpha <- lapply(kernals_alpha, 
                              function(x) as.numeric(determinant(x)$modulus))

        #initialize at middle of range
        which_theta_alpha_curr <- ceiling(length(discrete.theta.alpha.values) / 2)
        theta_alpha = discrete.theta.alpha.values[which_theta_alpha_curr]


    }

    if (!is.null(discrete.theta.beta.values)) {

        kernals_beta <- lapply(discrete.theta.beta.values,
                               function(x) cov_kern(distance = dist.space.mat, 
                                                    theta = x))
        kernals_inv_beta <- lapply(kernals_beta, solve)
        kernals_chol_beta <- lapply(kernals_beta, 
                               function(x) t(chol(x)))
        kernals_chol_inv_beta <- lapply(kernals_chol_beta, 
                                   function(x) solve(x))
        kernals_det_beta <- lapply(kernals_beta, 
                              function(x) as.numeric(determinant(x)$modulus))

        #initialize at middle of range
        which_theta_beta_curr <- ceiling(length(discrete.theta.beta.values) / 2)
        theta_beta = discrete.theta.beta.values[which_theta_beta_curr]

    }



    
    ####################################################
    ### Pre-calculate quantities to save computation ###
    ####################################################

    ### space assignment to spacetime vector
    space_to_spacetime_assign <- rep(unique(spacetime.id), each = N.space)
    
    ###Indicator Matrix for Time
    Gamma_time = matrix(0, nrow = N, ncol = N.time)
    for (i in 1:N){
        Gamma_time[i, time.id[i]] <- 1
    }
    Gamma_time = spam::as.spam(Gamma_time)
    GtG_time = spam::diag.spam(t(Gamma_time) %*% Gamma_time)
    
    ###Indicator Matrix for SpaceTime
    Gamma_space = matrix(0, nrow = N, ncol = N.spacetime * N.space)
    for (i in 1:N) {
        Gamma_space[i, (spacetime.id[i]-1) * N.space + space.id[i]] <- 1
    }
    Gamma_space = spam::as.spam(Gamma_space)
    GtG_space = spam::diag.spam(t(Gamma_space) %*% Gamma_space)
    Z_ID = (spacetime.id-1) * N.space + space.id
    
    
    ### Design matrix for temporal random effects
    X_W = NULL
    #GammatX_W = NULL
    for (time.i in (1:N.time)) {
        X.i = X[time.id == time.i]
      
        if (length (X.i)>0){
            X_W = c(X_W, t(X.i) %*% (X.i))
        } else {
            X_W = c(X_W, 0)
        }
    }
    
    ### Design matrix for spatial random effects
    X_S = NULL
    GammatX_S = NULL
    Space.perm = NULL
    SpaceTime.avail = NULL
    
    for (spacetime.i in 1:N.spacetime) {
        for (mon.i in 1:N.space) {
            use = which(space.id == mon.i & spacetime.id == spacetime.i) 
            SpaceTime.avail = rbind(SpaceTime.avail, 
                                    c(spacetime.i, mon.i, length(use)))
          
            if (length(use) != 0) {
                Space.perm = c(Space.perm, use)
                X.i = as.matrix(X[use])
                X_S = c(X_S, t(X.i) %*% X.i)
            } else {
                X_S = c(X_S, 0)
            }
        }
    }
    
    ### Spatial and Spatial-temporal coefficients
    if (!is.null(M)) MtM = t(M) %*% M
    if (!is.null(L)) LtL = t(L) %*% L
    XtX = t(cbind(1, X)) %*% cbind(1, X)
    XtX_inv = solve(XtX)


    ### Discretize CAR parameters
    nrho = 2000
    jun = 1 / sqrt(Q)
    dt = eigen(jun %*% A %*% jun)$values
    canrho = detpart = stats::qbeta(seq(0, 0.9999, length = nrho), 1, 1)
    for (j in 1:nrho) {
        detpart[j] = 1 / 2 * sum(log(1 - canrho[j] * dt)) 
    }
    
    
    #####################
    ### Initial Value ###
    #####################
    
    ###Initialize alpha0, beta0, gamma, and delta
    ###Initialize lambda if needed
    if (!is.null(L) & !is.null(M)) {
        fit = stats::lm(Y ~ cbind(X, L, M))
        alpha0 = stats::coef(fit)[1]
        beta0 = stats::coef(fit)[2]
        gamma = stats::coef(fit)[3:(ncol(L) + 2)]
        delta = stats::coef(fit)[(3 + ncol(L)):length(stats::coef(fit))]
        lambda_gamma = stats::var(gamma)
        lambda_delta = stats::var(delta)
        mu = alpha0 + beta0 * X + L %*% gamma + M %*% delta
    } 
    
    if (is.null(L) & !is.null(M)) {
        fit = stats::lm (Y ~ cbind(X, M))
        alpha0 = stats::coef(fit)[1]
        beta0 = stats::coef(fit)[2]
        delta = stats::coef(fit)[3:length(stats::coef(fit))]
        lambda_delta = stats::var(delta)
        mu = alpha0 + beta0 * X + M %*% delta
    }
    
    if (!is.null(L) & is.null(M)) {
        fit = stats::lm (Y ~ cbind(X, L))
        alpha0 = stats::coef(fit)[1]
        beta0 = stats::coef(fit)[2]
        gamma = stats::coef(fit)[3:(ncol(L)+2)]
        lambda_gamma = stats::var(gamma)
        mu = alpha0 + beta0 * X + L %*% gamma 
    }
    
    if (is.null(L) & is.null(M)) {
        fit = stats::lm (Y ~ X)
        alpha0 = stats::coef(fit)[1]
        beta0 = stats::coef(fit)[2]
        mu = alpha0 + beta0 * X 
    } 
    
    
    ###Initialize alpha_time and its CAR variance omega_alpha
    alpha_time = rep(0, N.time)
    omega_alpha = 0
    if (include.additive.weekly.resid) { 
        alpha_time = (1 / GtG_time) * t(Gamma_time) %*% (Y - mu)
        omega_alpha = as.numeric(stats::var(alpha_time, na.rm = T))
    }
    
    ###Initialize beta_time and its CAR variance omega_beta
    beta_time = rep(0, N.time)
    omega_beta = 0
    if (include.multiplicative.weekly.resid) { 
        RRR = Y - mu - alpha_time[time.id]
        beta_time=  1 / X_W * t(Gamma_time) %*% (X * RRR)
        omega_beta = as.numeric(stats::var(beta_time, na.rm = T))
    }
    
    ###Initialize alpha_spacetime and its spatial variance tau_alpha
    alpha_space = rep(0, N.spacetime * N.space)
    tau_alpha = 0
    if (include.additive.annual.resid) {
        RRR = Y - mu - alpha_time[time.id] - beta_time[time.id] * X
        alpha_space = as.numeric((1 / GtG_space) * t(Gamma_space) %*% RRR)
        tau_alpha = as.numeric(stats::var(alpha_space, na.rm = T))
    }
    
    ###Initialize beta_spacetime and its spatial variance tau_beta
    beta_space = rep(0, N.spacetime * N.space)
    tau_beta = 0
    if (include.multiplicative.annual.resid) {
        RRR = Y - mu - alpha_time[time.id] - beta_time[time.id] * X - alpha_space[Z_ID]
        beta_space = as.numeric(1 / X_S * t(Gamma_space) %*% (X * RRR))
        tau_beta = as.numeric(stats::var(beta_space, na.rm = T))
    }
    
    ###Initialize mean
    MMM = alpha0 + beta0 * X + 
        L %*% gamma + 
        M %*% delta + 
        alpha_time[time.id] + 
        beta_time[time.id] * X + 
        alpha_space[Z_ID] + 
        beta_space[Z_ID] * X
    
    ###Initializae sigma2
    sigma2 = stats::var(as.numeric(Y - MMM), na.rm = T)
    if (is.null(discrete.theta.alpha.values)) {
        theta_alpha = theta.alpha.init
        theta_beta = theta.beta.init
    }
    rho_alpha = rho.alpha.init
    rho_beta = rho.beta.init
    
    #############################
    ### MCMC saved parameters ###
    #############################
    ###Declare  matrices and vectors to save results
    
    ###Total number of samples in the end
    K = (n.iter - burn) / thin
    
    alpha0.save = rep(NA, K)
    beta0.save = rep(NA, K)
    
    gamma.save = matrix(NA, nrow = K, ncol = length(gamma))
    delta.save = matrix(NA, nrow = K, ncol = length(delta))
    
    lambda_gamma.save = rep(NA, K)
    lambda_delta.save = rep(NA, K)
    
    alpha_time.save = matrix(NA, nrow = K, ncol = N.time)
    beta_time.save = matrix(NA, nrow = K, ncol = N.time)
    
    alpha_space.save = matrix(NA, nrow = K, ncol = N.spacetime * N.space)
    beta_space.save = matrix(NA, nrow = K, ncol = N.spacetime * N.space)
    
    sigma2.save = rep(NA, K)
    theta_alpha.save = theta_beta.save = rep(NA, K)
    rho_alpha.save = rho_beta.save = rep(NA, K)
    
    omega_alpha.save = omega_beta.save = rep(NA, K)
    lambda_alpha.save = lambda_beta.save = rep(NA, K)
    tau_alpha.save = tau_beta.save = rep(NA, K)
    
    LL.save = rep(NA, K)
    
    Y.hat = rep(0, length(Y))
    
    tau.acc = c(0,0)
    theta.acc = c(0,0)

    if (!is.null(discrete.theta.alpha.values)) {
        which.theta.alpha.discrete <- rep(0, K)
    }

    if (!is.null(discrete.theta.beta.values)) {
        which.theta.beta.discrete <- rep(0, K)
    }
    
    
    ###########################
    ### BEGIN MCMC SAMPLING ###
    ###########################
    if (verbose == TRUE) {
      print ("#################################### ")
      print ("########    MCMC Begins!    ######## ")
      print ("#################################### ")
    }

    for (i in 1:n.iter) {

        if (verbose == TRUE) {
            if ((i %% verbose.iter) == 0) print(paste("     Iteration", i, "of", n.iter))
        }
    
        ###Update alpha0 and beta0
        MMM =  MMM - alpha0 - beta0 * X 
        RRR = Y - MMM
        XXX = t(cbind(1,X)) %*% RRR
        VVV = XtX_inv
        alpha0_beta0 = matrix(mvnfast::rmvn(1, VVV %*% XXX, sigma2 * VVV), ncol = 1)
        alpha0 = alpha0_beta0[1]
        beta0 = alpha0_beta0[2]
        MMM = MMM + alpha0 + beta0 * X
      
        #Update gamma
        if (!is.null(L)) {
            MMM = MMM - L %*% gamma
            RRR =  Y - MMM
            XXX = 1 / sigma2 * t(L) %*% RRR
            VVV = 1 / sigma2 * LtL + 1 / lambda_gamma * diag(N.Lmax)
            VVV = solve(VVV)
            gamma = matrix(mvnfast::rmvn(1, VVV %*% XXX, VVV), ncol = 1)
            MMM = MMM + L %*% gamma

            ##Update lambda_gamma
            lambda_gamma = 1 / stats::rgamma(1,  
                                             length(gamma) / 2 + sigma.a, 
                                             sum(gamma^2)/2 + sigma.b)
        }
      
        #Update delta
        if (!is.null(M)) {
            MMM = MMM - M %*% delta
            RRR =  Y - MMM
            XXX = 1 / sigma2 * t(M) %*% RRR
            VVV = 1 / sigma2 * MtM + 1 / lambda_delta * diag(N.Mmax)
            VVV = solve(VVV)
            delta = matrix(mvnfast::rmvn(1, VVV %*% XXX, VVV), ncol = 1)
            MMM = MMM + M %*% delta
        
            ##Update lambda_delta
            lambda_delta = 1 / stats::rgamma(1, 
                                             length(delta) / 2 + sigma.a, 
                                             sum(delta ^ 2) / 2 + sigma.b)
        }
      
        #Update residual error sigma2
        RRR = Y - MMM
        sigma2 = 1 / stats::rgamma(1, length(RRR) / 2 + sigma.a, sum(RRR ^ 2) / 2 + sigma.b)
       
        #Update spatial intercepts and parameters if GP
        if (include.additive.annual.resid & !nngp & is.null(discrete.theta.alpha.values)) {
            MMM = MMM - alpha_space[Z_ID]
            RRR = Y - MMM
            XXX = 1 / sigma2 * t(Gamma_space) %*% RRR
            kern = cov_kern(distance = dist.space.mat, theta = theta_alpha)
            kern_inv = solve(kern)
            SSS = tau_alpha * kern
            SSS_inv = (1 / tau_alpha) * kern_inv
            for (st in unique(spacetime.id)) {
                GtG_space_st = GtG_space[space_to_spacetime_assign == st]
                XXX_st = XXX[space_to_spacetime_assign == st]
                VVV_st = diag(1 / sigma2 * GtG_space_st) + SSS_inv
                VVV_inv_st = solve(VVV_st)
                alpha_space_st = as.vector(mvnfast::rmvn(1, VVV_inv_st %*% XXX_st, 
                                                         VVV_inv_st))
                alpha_space[space_to_spacetime_assign == st] = alpha_space_st
            }
            MMM = MMM + alpha_space[Z_ID]
      
            #update tau_alpha
            SSS = 0
            for (st in unique(spacetime.id)) {
                alpha_space_st = alpha_space[space_to_spacetime_assign == st]
                SSS_st = t(alpha_space_st) %*% kern_inv %*% alpha_space_st
                SSS = SSS + SSS_st
            }
            tau_alpha = 1 / stats::rgamma(1, 
                                          N.space * N.spacetime / 2 + tau.alpha.a, 
                                          SSS / 2 + tau.alpha.b)
      
            #Update theta_alpha
            theta.prop = stats::rlnorm(1, 
                                       log(theta_alpha), 
                                       theta.alpha.tune)
            SSS.curr = tau_alpha * kern
            SSS.prop = tau_alpha * cov_kern(distance = dist.space.mat, 
                                            theta = theta.prop)
      
            lik.curr = 0
            lik.prop = 0

            for (st in unique(spacetime.id)) {
                alpha_space_st = alpha_space[space_to_spacetime_assign == st]
                lik.prop_st = mvtnorm::dmvnorm(alpha_space_st, 
                                               rep(0, N.space), 
                                               SSS.prop, 
                                               log = T)
                lik.curr_st = mvtnorm::dmvnorm(alpha_space_st,
                                               rep(0, N.space), 
                                               SSS.curr, 
                                               log = T)
                lik.prop = lik.prop + lik.prop_st
                lik.curr = lik.curr + lik.curr_st
            }
      
            ratio = lik.prop + 
                stats::dgamma(theta.prop, 
                              theta.alpha.a, 
                              theta.alpha.b, 
                              log = T) + 
                    log(theta.prop) -
                    lik.curr - 
                    stats::dgamma(theta_alpha, 
                                  theta.alpha.a, 
                                  theta.alpha.b, 
                                  log = T) - 
                    log(theta_alpha)
            if (log(stats::runif(1)) < ratio) {
                theta_alpha = theta.prop
                theta.acc[1] = theta.acc[1] + 1
            }
        }

        #Update spatial intercepts and parameters if NNGP
        if (include.additive.annual.resid & nngp & is.null(discrete.theta.alpha.values)) {
            MMM = MMM - alpha_space[Z_ID]
            RRR = Y - MMM
            XXX = 1 / sigma2 * t(Gamma_space) %*% RRR

            alpha_space <- rep(0, N.space * N.spacetime)

            #calculate separately for each spacetime
            for (st in unique(spacetime.id)) {
                XXX_st <- XXX[space_to_spacetime_assign == st, , drop = FALSE]
                XXX_st_ord <- XXX_st[coord_ordering, , drop = FALSE]
                GtG_space_st <- GtG_space[space_to_spacetime_assign == st]
                GtG_space_st_ord <- GtG_space_st[coord_ordering]
                coords_st_ord <- ordered_coords

                alpha_space_st <- alpha_space[space_to_spacetime_assign == st]


                #iterate through ordered conditionals
                SSS_n <- tau_alpha * cov_kern(distance = 0, theta = theta_alpha)
                VVV_n <- (1 / sigma2 * GtG_space_st_ord[length(alpha_space_st)]) + (1 / SSS_n)
                joint_cov_n <- 1 / VVV_n
                joint_mean_n <- joint_cov_n * XXX_st_ord[length(alpha_space_st)]
                alpha_space_st[length(alpha_space_st)] <- stats::rnorm(1, joint_mean_n, sqrt(joint_cov_n))

                for (j in (length(alpha_space_st) - 1):1) {
                    j_and_neighbors <- c(j, neighbors[[j]])
                    dist_space_mat_j <- dist_matrices[[j]]
                    SSS_j <- tau_alpha * cov_kern(distance = dist_space_mat_j, 
                                                  theta = theta_alpha)
                    VVV_j <- diag(1 / sigma2 * GtG_space_st_ord[j_and_neighbors]) + solve(SSS_j)
                    joint_cov <- solve(VVV_j)
                    joint_mean <- joint_cov %*% XXX_st_ord[j_and_neighbors]
                    sigma_12_inv_sigma_22 <- joint_cov[1, -1] %*% solve(joint_cov[-1, -1])
                    cond_mean <- joint_mean[1] - 
                        sigma_12_inv_sigma_22 %*% 
                        (alpha_space_st[neighbors[[j]]] - 
                         joint_mean[-1])
                    cond_cov <- joint_cov[1, 1] - 
                        sigma_12_inv_sigma_22 %*% 
                        joint_cov[-1, 1]
                    alpha_space_st[j] <- stats::rnorm(1, cond_mean, sqrt(cond_cov))
                }
                #reverse order back to original
                alpha_space_st <- alpha_space_st[coord_reverse_ordering]
                alpha_space[space_to_spacetime_assign == st] <- alpha_space_st
            }

            MMM = MMM + alpha_space[Z_ID]
      
      
            #Update tau_alpha
            tau_prop = stats::rlnorm(1, 
                                     log(tau_alpha), 
                                     tau.alpha.tune)
            lik_prop <- 0
            lik_cur <- 0
            for (st in unique(spacetime.id)) {
                alpha_space_st <- alpha_space[space_to_spacetime_assign == st]
                alpha_space_st <- alpha_space_st[coord_ordering]
                lik_prop <- lik_prop +
                    sum(dnngp(y = alpha_space_st,
                              ordered_coords = ordered_coords,
                              neighbors = neighbors,
                              dist_matrices = dist_matrices,
                              phi = tau_prop,
                              r = theta_alpha,
                              cov_kern = cov_kern,
                              log = T))
                lik_cur <- lik_cur +
                    sum(dnngp(y = alpha_space_st,
                              ordered_coords = ordered_coords,
                              neighbors = neighbors,
                              dist_matrices = dist_matrices,
                              phi = tau_alpha,
                              r = theta_alpha,
                              cov_kern = cov_kern,
                              log = T))
            }


      
            ratio = lik_prop + 
                stats::dgamma(tau_prop, 
                              tau.alpha.a, 
                              tau.alpha.b, 
                              log = T) + 
                log(tau_prop) -
                lik_cur - 
                stats::dgamma(tau_alpha, 
                              tau.alpha.a, 
                              tau.alpha.b, 
                              log = T) - 
                log(tau_alpha)

            if (log(stats::runif(1)) < ratio) {
                tau_alpha = tau_prop
                tau.acc[1] = tau.acc[1] + 1
            }

            #Update theta_alpha
            theta_prop = stats::rlnorm(1, 
                                       log(theta_alpha), 
                                       theta.alpha.tune)
            lik_prop <- 0
            lik_cur <- 0
            for (st in unique(spacetime.id)) {
                alpha_space_st <- alpha_space[space_to_spacetime_assign == st]
                alpha_space_st <- alpha_space_st[coord_ordering]
                lik_prop <- lik_prop +
                    sum(dnngp(y = alpha_space_st,
                              ordered_coords = ordered_coords,
                              neighbors = neighbors,
                              dist_matrices = dist_matrices,
                              phi = tau_alpha,
                              r = theta_prop,
                              cov_kern = cov_kern,
                              log = T))
                lik_cur <- lik_cur +
                    sum(dnngp(y = alpha_space_st,
                              ordered_coords = ordered_coords,
                              neighbors = neighbors,
                              dist_matrices = dist_matrices,
                              phi = tau_alpha,
                              r = theta_alpha,
                              cov_kern = cov_kern,
                              log = T))
            }

            ratio = lik_prop + 
                stats::dgamma(theta_prop, 
                              theta.alpha.a, 
                              theta.alpha.b, 
                              log = T) + 
                log(theta_prop) -
                lik_cur - 
                stats::dgamma(theta_alpha, 
                              theta.alpha.a, 
                              theta.alpha.b, 
                              log = T) - 
                log(theta_alpha)

            if (log(stats::runif(1)) < ratio) {
                theta_alpha = theta_prop
                theta.acc[1] = theta.acc[1] + 1
            }
        }

        #Update spatial intercept if discrete theta (GP range parameter)
        if (include.additive.annual.resid & !nngp & !is.null(discrete.theta.alpha.values)) {

            kern_curr = kernals_alpha[[which_theta_alpha_curr]]
            kern_inv_curr = kernals_inv_alpha[[which_theta_alpha_curr]]
            kern_chol_curr = kernals_chol_alpha[[which_theta_alpha_curr]]
            kern_chol_inv_curr = kernals_chol_inv_alpha[[which_theta_alpha_curr]]
            kern_det_curr = kernals_det_alpha[[which_theta_alpha_curr]]


            MMM = MMM - alpha_space[Z_ID] 
            RRR = Y - MMM
            XXX = 1 / sigma2 * t(Gamma_space) %*% RRR
            SSS = tau_alpha * kern_curr
            SSS_inv = (1 / tau_alpha) * kern_inv_curr
            for (st in unique(spacetime.id)) {
                GtG_space_st = GtG_space[space_to_spacetime_assign == st]
                XXX_st = XXX[space_to_spacetime_assign == st]
                VVV_st = diag(1 / sigma2 * GtG_space_st) + SSS_inv
                VVV_inv_st = solve(VVV_st)
                alpha_space_st = as.vector(mvnfast::rmvn(1, VVV_inv_st %*% XXX_st, 

                                                         VVV_inv_st))
                alpha_space[space_to_spacetime_assign == st] <- alpha_space_st
            }
            MMM = MMM + alpha_space[Z_ID] 

      
            #update tau_alpha
            SSS = 0
            for (st in unique(spacetime.id)) {
                alpha_space_st = alpha_space[space_to_spacetime_assign == st]
                SSS_st = t(alpha_space_st) %*% kern_inv_curr %*% alpha_space_st
                SSS = SSS + SSS_st
            }
            tau_alpha = 1 / stats::rgamma(1, N.space * N.spacetime /2 + tau.alpha.a, SSS / 2 + tau.alpha.b)
      
            #update theta_alpha

            #Update theta_alpha - mh jump 
            #jump is which direction to jump in the discrete theta value indices
            #adjustment is the mh proposal likelihood adjustment
            if (!discrete.theta.gibbs) {

                if (which_theta_alpha_curr == 1) { 
                    jump <- 1
                    lik_jump_alpha_curr_to_prop <- 1
                    lik_jump_alpha_prop_to_curr <- 0.5
                } else if (which_theta_alpha_curr == length(discrete.theta.alpha.values)) {
                    jump <- -1
                    lik_jump_alpha_curr_to_prop <- 1
                    lik_jump_alpha_prop_to_curr <- 0.5
                } else {
                    jump <- sample(c(-1, 1), 1)
                    lik_jump_alpha_curr_to_prop <- 0.5
                    if ((jump + which_theta_alpha_curr) %in% c(1, length(discrete.theta.alpha.values))) {
                        lik_jump_alpha_prop_to_curr <- 1
                    } else {
                        lik_jump_alpha_prop_to_curr <- 0.5
                    }
                }

                which_theta_alpha_prop <- which_theta_alpha_curr + jump
                theta_alpha_prop = discrete.theta.alpha.values[which_theta_alpha_prop]
                kern_prop = kernals_alpha[[which_theta_alpha_prop]]
                kern_inv_prop = kernals_inv_alpha[[which_theta_alpha_prop]]
                kern_chol_prop = kernals_chol_alpha[[which_theta_alpha_prop]]
                kern_chol_inv_prop = kernals_chol_inv_alpha[[which_theta_alpha_prop]]
                kern_det_prop = kernals_det_alpha[[which_theta_alpha_prop]]

                SSS_chol_curr = sqrt(tau_alpha) * kern_chol_curr
                SSS_det_curr = ncol(kern_curr) * log(tau_alpha) + kern_det_curr
                SSS_chol_inv_curr = (1 / sqrt(tau_alpha)) * kern_chol_inv_curr

                SSS_chol_prop = sqrt(tau_alpha) * kern_chol_prop
                SSS_det_prop = ncol(kern_prop) * log(tau_alpha) + kern_det_prop
                SSS_chol_inv_prop = (1 / sqrt(tau_alpha)) * kern_chol_inv_prop

                lik.curr = 0
                lik.prop = 0

                for (st in unique(spacetime.id)) {
                    alpha_space_st = alpha_space[space_to_spacetime_assign == st]
                    lik.curr_st = d_mvn_chol_uvn(alpha_space_st, 
                                                 SSS_chol_inv_curr, 
                                                 SSS_det_curr)
                    lik.prop_st = d_mvn_chol_uvn(alpha_space_st, 
                                                 SSS_chol_inv_prop, 
                                                 SSS_det_prop)
                    lik.curr = lik.curr + lik.curr_st
                    lik.prop = lik.prop + lik.prop_st
                }


                ratio = lik.prop + 
                    stats::dgamma(theta_alpha_prop, 
                                  theta.alpha.a, 
                                  theta.alpha.b, 
                                  log = T) + 
                    log(lik_jump_alpha_prop_to_curr) -
                    lik.curr - 
                    stats::dgamma(theta_alpha, 
                                  theta.alpha.a, 
                                  theta.alpha.b, 
                                  log = T) -
                    log(lik_jump_alpha_curr_to_prop)


                if(log(stats::runif(1)) < ratio) {
                    which_theta_alpha_curr <- which_theta_alpha_prop
                    theta_alpha = discrete.theta.alpha.values[which_theta_alpha_curr]
                    theta.acc[1] = theta.acc[1] + 1
                }
            }

            #Update theta_alpha - gibbs
            if (discrete.theta.gibbs) {
                #jump is which direction to jump in the discrete theta value indices
                #adjustment is the mh proposal likelihood adjustment


                SSS_chol = lapply(kernals_chol_alpha, 
                                       function(x) sqrt(tau_alpha) * x)
                SSS_det = lapply(kernals_det_alpha,
                                 function(x) N.space * log(tau_alpha) + x)
                SSS_chol_inv = lapply(kernals_chol_inv_alpha,
                                      function(x) (1 / sqrt(tau_alpha)) * x)


                lik = rep(0, length(discrete.theta.alpha.values))



                for (st in unique(spacetime.id)) {
                    alpha_space_st = alpha_space[space_to_spacetime_assign == st]
                    lik <- lik + mapply(function(x, y) d_mvn_chol_uvn(alpha_space_st, x, y),
                                        SSS_chol_inv,
                                        SSS_det)
                }
                lik <- lik + stats::dgamma(discrete.theta.alpha.values, 
                                           theta.alpha.a, 
                                           theta.alpha.b, 
                                           log = T)
                theta_alpha <- sample(x = discrete.theta.alpha.values, 
                                     size = 1, 
                                     prob = exp(lik - max(lik)))
                which_theta_alpha_curr <- which(discrete.theta.alpha.values == theta_alpha)
            }

        }
      
        #Update spatial coefficent for AOD if GP
        if (include.multiplicative.annual.resid & !nngp & is.null(discrete.theta.beta.values)) {

            MMM = MMM - beta_space[Z_ID] * X
            RRR = Y - MMM
            XXX = 1 / sigma2 * t(Gamma_space) %*% (X * RRR)
            kern = cov_kern(distance = dist.space.mat, 
                            theta = theta_beta)
            kern_inv = solve(kern)
            SSS = tau_beta * kern
            SSS_inv = (1 / tau_beta) * kern_inv
            for (st in unique(spacetime.id)) {
                X_S_st = X_S[space_to_spacetime_assign == st]
                XXX_st = XXX[space_to_spacetime_assign == st]
                VVV_st = diag(1 / sigma2 * X_S_st)
                VVV_st = VVV_st + SSS_inv
                VVV_inv_st = solve(VVV_st)
                beta_space_st = as.vector(mvnfast::rmvn(1, VVV_inv_st %*% XXX_st, 
                                                        VVV_inv_st))
                beta_space[space_to_spacetime_assign == st] <- beta_space_st
            }
            MMM = MMM + beta_space[Z_ID] * X

      
            #update tau_beta
            SSS = 0
            for (st in unique(spacetime.id)) {
                beta_space_st = beta_space[space_to_spacetime_assign == st]
                SSS_st = t(beta_space_st) %*% 
                    kern_inv %*% 
                    beta_space_st
                SSS = SSS + SSS_st
            }
            tau_beta = 1 / stats::rgamma(1, N.space * N.spacetime /2 + tau.beta.a, SSS / 2 + tau.beta.b)
      
            #Update theta_beta
            theta.prop = stats::rlnorm(1, log(theta_beta), theta.beta.tune)
            SSS.curr = tau_beta * kern
            SSS.prop = tau_beta * cov_kern(distance = dist.space.mat, 
                                           theta = theta.prop)

            lik.curr = 0
            lik.prop = 0
            for (st in unique(spacetime.id)) {
                beta_space_st = beta_space[space_to_spacetime_assign == st]
                lik.prop_st = mvtnorm::dmvnorm(beta_space_st, 
                                               rep(0, N.space), 
                                               SSS.prop, 
                                               log = T)
                lik.curr_st = mvtnorm::dmvnorm(beta_space_st,
                                               rep(0, N.space), 
                                               SSS.curr, 
                                               log = T)
                lik.prop = lik.prop + lik.prop_st
                lik.curr = lik.curr + lik.curr_st
            }
      
            ratio = lik.prop + 
                stats::dgamma(theta.prop,  
                              theta.beta.a, 
                              theta.beta.b, 
                              log = T) + 
                log(theta.prop) -
                lik.curr - 
                stats::dgamma(theta_beta, 
                              theta.beta.a, 
                              theta.beta.b, 
                              log = T) - 
                log(theta_beta)

            if(log(stats::runif(1)) < ratio) {
                theta_beta = theta.prop
                theta.acc[2] = theta.acc[2] + 1
            }
        }

        #Update spatial coefficent for AOD if NNGP
        if (include.multiplicative.annual.resid & nngp & is.null(discrete.theta.beta.values)) {
            MMM = MMM - beta_space[Z_ID] * X
            RRR = Y - MMM
            XXX = 1 / sigma2 * t(Gamma_space) %*% (X * RRR)

            beta_space <- rep(0, N.space * N.spacetime)
            #calculate separately for each spacetime
            for (st in unique(spacetime.id)) {
                XXX_st <- XXX[space_to_spacetime_assign == st, , drop = FALSE]
                XXX_st_ord <- XXX_st[coord_ordering, , drop = FALSE]
                X_S_st <- X_S[space_to_spacetime_assign == st]
                X_S_st_ord <- X_S_st[coord_ordering]
                coords_st_ord <- ordered_coords

                beta_space_st <- beta_space[space_to_spacetime_assign == st]
                SSS_n <- tau_beta * cov_kern(distance = 0, 
                                             theta = theta_beta)
                VVV_n <- (1 / sigma2 * X_S_st_ord[length(beta_space_st)]) + (1 / SSS_n)
                # iterate through ordered conditionals
                joint_cov_n <- 1 / VVV_n
                joint_mean_n <- joint_cov_n * XXX_st_ord[length(beta_space_st)]
                beta_space_st[length(beta_space_st)] <- stats::rnorm(1, joint_mean_n, sqrt(joint_cov_n))

                for (j in (length(beta_space_st) - 1):1) {
                    j_and_neighbors <- c(j, neighbors[[j]])
                    dist_space_mat_j <- dist_matrices[[j]]
                    SSS_j <- tau_beta * cov_kern(distance = dist_space_mat_j, 
                                                 theta = theta_beta)
                    VVV_j = diag(1 / sigma2 * X_S_st_ord[j_and_neighbors]) + solve(SSS_j)
                    joint_cov = solve(VVV_j)
                    joint_mean <- joint_cov %*% XXX_st_ord[j_and_neighbors]
                    sigma_12_inv_sigma_22 <- joint_cov[1, -1] %*% solve(joint_cov[-1, -1])
                    cond_mean <- joint_mean[1] - 
                        sigma_12_inv_sigma_22 %*% 
                        (beta_space_st[neighbors[[j]]] - 
                         joint_mean[-1])
                    cond_cov <- joint_cov[1, 1] - 
                        sigma_12_inv_sigma_22 %*% 
                        joint_cov[-1, 1]
                    beta_space_st[j] <- stats::rnorm(1, cond_mean, sqrt(cond_cov))
                }
                #reverse order back to original
                beta_space_st <- beta_space_st[coord_reverse_ordering]
                beta_space[space_to_spacetime_assign == st] <- beta_space_st
            }

            MMM = MMM + beta_space[Z_ID] * X


            #Update tau_beta
            tau_prop = stats::rlnorm(1, 
                                     log(tau_beta),
                                     tau.beta.tune)
            lik_prop <- 0
            lik_cur <- 0
            for (st in unique(spacetime.id)) {
                beta_space_st <- beta_space[space_to_spacetime_assign == st]
                beta_space_st <- beta_space_st[coord_ordering]
                lik_prop <- lik_prop +
                    sum(dnngp(y = beta_space_st,
                              ordered_coords = ordered_coords,
                              neighbors = neighbors,
                              dist_matrices = dist_matrices,
                              phi = tau_prop,
                              r = theta_beta,
                              cov_kern = cov_kern, 
                              log = T))
                lik_cur <- lik_cur +
                    sum(dnngp(y = beta_space_st,
                              ordered_coords = ordered_coords,
                              neighbors = neighbors,
                              dist_matrices = dist_matrices,
                              phi = tau_beta,
                              r = theta_beta,
                              cov_kern = cov_kern,
                              log = T))
            }


            ratio = lik_prop + 
                stats::dgamma(tau_prop, 
                              tau.beta.a, 
                              tau.beta.b, 
                              log = T) + 
                log(tau_prop) -
                lik_cur - 
                stats::dgamma(tau_beta, 
                              tau.beta.a, 
                              tau.beta.b, 
                              log = T) - 
                log(tau_beta)

            if (log(stats::runif(1)) < ratio) {
                tau_beta = tau_prop
                tau.acc[2] = tau.acc[2] + 1
            }

            #Update theta_beta
            theta_prop = stats::rlnorm(1, 
                                       log(theta_beta), 
                                       theta.beta.tune)
            lik_prop <- 0
            lik_cur <- 0
            for (st in unique(spacetime.id)) {
                beta_space_st <- beta_space[space_to_spacetime_assign == st]
                beta_space_st <- beta_space_st[coord_ordering]
                lik_prop <- lik_prop +
                    sum(dnngp(y = beta_space_st,
                              ordered_coords = ordered_coords,
                              neighbors = neighbors,
                              dist_matrices = dist_matrices,
                              phi = tau_beta,
                              r = theta_prop,
                              cov_kern = cov_kern,
                              log = T))
                lik_cur <- lik_cur +
                    sum(dnngp(y = beta_space_st,
                              ordered_coords = ordered_coords,
                              neighbors = neighbors,
                              dist_matrices = dist_matrices,
                              phi = tau_beta,
                              r = theta_beta,
                              cov_kern = cov_kern,
                              log = T))
            }


        
            ratio = lik_prop + 
                stats::dgamma(theta_prop, 
                              theta.beta.a, 
                              theta.beta.b, 
                              log = T) + 
                log(theta_prop) - 
                lik_cur - 
                stats::dgamma(theta_beta, 
                              theta.beta.a, 
                              theta.beta.b, 
                              log = T) - 
                log(theta_beta)

            if (log(stats::runif(1)) < ratio) {
                theta_beta = theta_prop
                theta.acc[2] = theta.acc[2] + 1
            }
        }

        #Update spatial coefficent for AOD if discrete theta (GP range parameter)
        if (include.multiplicative.annual.resid & !nngp & !is.null(discrete.theta.beta.values)) {

            kern_curr = kernals_beta[[which_theta_beta_curr]]
            kern_inv_curr = kernals_inv_beta[[which_theta_beta_curr]]
            kern_chol_curr = kernals_chol_beta[[which_theta_beta_curr]]
            kern_chol_inv_curr = kernals_chol_inv_beta[[which_theta_beta_curr]]
            kern_det_curr = kernals_det_beta[[which_theta_beta_curr]]


            MMM = MMM - beta_space[Z_ID] * X
            RRR = Y - MMM
            XXX = 1 / sigma2 * t(Gamma_space) %*% (X * RRR)
            SSS = tau_beta * kern_curr
            SSS_inv = (1 / tau_beta) * kern_inv_curr
            for (st in unique(spacetime.id)) {
                X_S_st = X_S[space_to_spacetime_assign == st]
                XXX_st = XXX[space_to_spacetime_assign == st]
                VVV_st = diag(1 / sigma2 * X_S_st)
                VVV_st = VVV_st + (1 / tau_beta) * kern_inv_curr
                VVV_inv_st = solve(VVV_st)
                beta_space_st = as.vector(mvnfast::rmvn(1, VVV_inv_st %*% XXX_st, 
                                                        VVV_inv_st))
                beta_space[space_to_spacetime_assign == st] <- beta_space_st
            }
            MMM = MMM + beta_space[Z_ID] * X

      
            #update tau_beta
            SSS = 0
            for (st in unique(spacetime.id)) {
                beta_space_st = beta_space[space_to_spacetime_assign == st]
                SSS_st = t(beta_space_st) %*% 
                    kern_inv_curr %*% 
                    beta_space_st
                SSS = SSS + SSS_st
            }
            tau_beta = 1 / stats::rgamma(1, N.space * N.spacetime /2 + tau.beta.a, SSS / 2 + tau.beta.b)
      
            #update theta_beta

            #Update theta_beta - mh jump 
            #jump is which direction to jump in the discrete theta value indices
            #adjustment is the mh proposal likelihood adjustment
            if (!discrete.theta.gibbs) {

                if (which_theta_beta_curr == 1) { 
                    jump <- 1
                    lik_jump_beta_curr_to_prop <- 1
                    lik_jump_beta_prop_to_curr <- 0.5
                } else if (which_theta_beta_curr == length(discrete.theta.beta.values)) {
                    jump <- -1
                    lik_jump_beta_curr_to_prop <- 1
                    lik_jump_beta_prop_to_curr <- 0.5
                } else {
                    jump <- sample(c(-1, 1), 1)
                    lik_jump_beta_curr_to_prop <- 0.5
                    if ((jump + which_theta_beta_curr) %in% c(1, length(discrete.theta.beta.values))) {
                        lik_jump_beta_prop_to_curr <- 1
                    } else {
                        lik_jump_beta_prop_to_curr <- 0.5
                    }
                }

                which_theta_beta_prop <- which_theta_beta_curr + jump
                theta_beta_prop = discrete.theta.beta.values[which_theta_beta_prop]
                kern_prop = kernals_beta[[which_theta_beta_prop]]
                kern_inv_prop = kernals_inv_beta[[which_theta_beta_prop]]
                kern_chol_prop = kernals_chol_beta[[which_theta_beta_prop]]
                kern_chol_inv_prop = kernals_chol_inv_beta[[which_theta_beta_prop]]
                kern_det_prop = kernals_det_beta[[which_theta_beta_prop]]

                SSS_chol_curr = sqrt(tau_beta) * kern_chol_curr
                SSS_det_curr = ncol(kern_curr) * log(tau_beta) + kern_det_curr
                SSS_chol_inv_curr = (1 / sqrt(tau_beta)) * kern_chol_inv_curr

                SSS_chol_prop = sqrt(tau_beta) * kern_chol_prop
                SSS_det_prop = ncol(kern_prop) * log(tau_beta) + kern_det_prop
                SSS_chol_inv_prop = (1 / sqrt(tau_beta)) * kern_chol_inv_prop

                lik.curr = 0
                lik.prop = 0

                for (st in unique(spacetime.id)) {
                    beta_space_st = beta_space[space_to_spacetime_assign == st]
                    lik.curr_st = d_mvn_chol_uvn(beta_space_st, 
                                                 SSS_chol_inv_curr, 
                                                 SSS_det_curr)
                    lik.prop_st = d_mvn_chol_uvn(beta_space_st, 
                                                 SSS_chol_inv_prop, 
                                                 SSS_det_prop)
                    lik.curr = lik.curr + lik.curr_st
                    lik.prop = lik.prop + lik.prop_st
                }

                ratio = lik.prop + 
                    stats::dgamma(theta_beta_prop, 
                                  theta.beta.a, 
                                  theta.beta.b, 
                                  log = T) + 
                    log(lik_jump_beta_prop_to_curr) -
                    lik.curr - 
                    stats::dgamma(theta_beta, 
                                  theta.beta.a, 
                                  theta.beta.b, 
                                  log = T) -
                    log(lik_jump_beta_curr_to_prop)


                if(log(stats::runif(1)) < ratio) {
                    which_theta_beta_curr <- which_theta_beta_prop
                    theta_beta = discrete.theta.beta.values[which_theta_beta_curr]
                    theta.acc[2] = theta.acc[2] + 1
                }
            }

            #Update theta_beta - gibbs
            if (discrete.theta.gibbs) {
                #jump is which direction to jump in the discrete theta value indices
                #adjustment is the mh proposal likelihood adjustment


                SSS_chol = lapply(kernals_chol_beta, 
                                       function(x) sqrt(tau_beta) * x)
                SSS_det = lapply(kernals_det_beta,
                                 function(x) N.space * log(tau_beta) + x)
                SSS_chol_inv = lapply(kernals_chol_inv_beta,
                                      function(x) (1 / sqrt(tau_beta)) * x)


                lik = rep(0, length(discrete.theta.beta.values))



                for (st in unique(spacetime.id)) {
                    beta_space_st = beta_space[space_to_spacetime_assign == st]
                    lik <- lik + mapply(function(x, y) d_mvn_chol_uvn(beta_space_st, x, y),
                                        SSS_chol_inv,
                                        SSS_det)
                }

                lik <- lik + stats::dgamma(discrete.theta.beta.values, 
                                           theta.beta.a, 
                                           theta.beta.b, 
                                           log = T)

                theta_beta <- sample(x = discrete.theta.beta.values, 
                                     size = 1, 
                                     prob = exp(lik - max(lik)))
                which_theta_beta_curr <- which(discrete.theta.beta.values == theta_beta)
            }

        }
      
        #Update temporal intercepts and its parameters
        if (include.additive.weekly.resid) {
            MMM = MMM - alpha_time[time.id]
            RRR = Y - MMM
            XXX = 1 / sigma2 * t(Gamma_time) %*% RRR
            VVV = diag(1 / sigma2 * GtG_time) + 1 / omega_alpha * (Q - rho_alpha * A)
            VVV = solve(VVV)
            alpha_time = mvnfast::rmvn(1, VVV %*% XXX, VVV)[1,]
            MMM = MMM + alpha_time[time.id]
      
            #Update omega_alpha
            SS1 = alpha_time %*% Q %*% alpha_time
            SS2 = alpha_time %*% A %*% alpha_time
            omega_alpha = 1 / stats::rgamma(1, N.time / 2 + omega.a, (SS1 - rho_alpha * SS2) / 2 + omega.b)
            
            #Update rho_alpha
            R = detpart + 0.5 * omega_alpha * canrho * c(SS2)
            rho_alpha = sample(canrho, 1, prob = exp(R - max(R)))
        }

        #Update temporal AOD coefficients
        if (include.multiplicative.weekly.resid) {
            MMM = MMM - beta_time[time.id] * X
            RRR = Y - MMM
            XXX = 1 / sigma2 * t(Gamma_time) %*% (X * RRR)
            VVV = diag(1 / sigma2 * X_W) + 1 / omega_beta * (Q - rho_beta * A)
            VVV = solve(VVV)
            beta_time = mvnfast::rmvn(1, VVV %*% XXX, VVV)[1,]
            MMM = MMM + beta_time[time.id] * X
        
            #Update omega_beta
            SS1 = beta_time %*% Q %*% beta_time
            SS2 = beta_time %*% A %*% beta_time
            omega_beta = 1 / stats::rgamma(1, N.time / 2 + omega.a, (SS1 - rho_beta * SS2) / 2 + omega.b)
            
            #Update rho_beta
            R = detpart + 0.5 * omega_beta * canrho * c(SS2)
            rho_beta = sample(canrho, 1, prob = exp(R - max(R)))
        }
    
     
        ###Save Samples 
     
        if (i > burn & i %% thin == 0) {
            k = (i - burn) / thin

            #Save statistics
            LL.save[k] = sum(-2 * stats::dnorm(Y, MMM, sqrt(sigma2), log = T))
  
            alpha0.save[k] = alpha0
            beta0.save[k] = beta0
            gamma.save[k, ] = gamma
            delta.save[k, ] = delta
       
            lambda_gamma.save[k] = lambda_gamma
            lambda_delta.save[k] = lambda_delta
       
            alpha_time.save[k,] = alpha_time
            beta_time.save[k,] = beta_time
       
            alpha_space.save[k,] = alpha_space
            beta_space.save[k,] = beta_space
       
            sigma2.save[k] = sigma2
            theta_alpha.save[k] = theta_alpha
            theta_beta.save[k] = theta_beta
            rho_alpha.save[k] = rho_alpha
            rho_beta.save[k] = rho_beta
       
            omega_alpha.save[k] = omega_alpha
            omega_beta.save[k] = omega_beta
            tau_alpha.save[k] = tau_alpha
            tau_beta.save[k] = tau_beta
       
            if (!is.null(discrete.theta.alpha.values)) {
                which.theta.alpha.discrete[k] = which_theta_alpha_curr
            }

            if (!is.null(discrete.theta.beta.values)) {
                which.theta.beta.discrete[k] = which_theta_beta_curr
            }

            Y.hat = Y.hat + MMM / K

        }

    } ## End of MCMC iterations
    
    if (verbose == TRUE) {
        print("#################################### ")
    }

    ##Process for output
    delta.save = data.frame(delta.save)
    names(delta.save) = colnames(M)
    gamma.save = data.frame(gamma.save);
    names(gamma.save) = colnames(L)
    
    alpha_time.save = data.frame(time.id = 1:N.time, t(alpha_time.save))
    names(alpha_time.save) = c("time.id", 
                               paste0("Sample", 1:K)) 
    row.names(alpha_time.save) = NULL
    beta_time.save = data.frame(time.id = 1:N.time, 
                                t(beta_time.save))
    names(beta_time.save) = c("time.id", 
                              paste0("Sample", 1:K)) 
    row.names(alpha_time.save) = NULL
    
    alpha_space.save = data.frame(space.id = rep(1:N.space, N.spacetime),
                                  spacetime.id = rep(1:N.spacetime, each = N.space),  
                                  t(alpha_space.save))
    names(alpha_space.save) = c("space.id", 
                                "spacetime.id", 
                                paste0("Sample", 1:K)) 
    row.names(alpha_space.save) = NULL
    beta_space.save = data.frame(space.id = rep(1:N.space, N.spacetime),
                                 spacetime.id = rep(1:N.spacetime, each = N.space),  
                                t(beta_space.save))
    names(beta_space.save) = c("space.id", "spacetime.id", paste0("Sample", 1:K))
    row.names(beta_space.save) = NULL
    
    other.save = data.frame(alpha0 = alpha0.save, 
                            beta0 = beta0.save,
                            sigma2 = sigma2.save,
                            lambda.delta = lambda_delta.save, 
                            lambda.gamma = lambda_gamma.save,
                            theta.alpha = theta_alpha.save, 
                            theta.beta = theta_beta.save, 
                            tau.alpha = tau_alpha.save, 
                            tau.beta = tau_beta.save, 
                            rho.alpha = rho_alpha.save, 
                            rho.beta = rho_beta.save, 
                            omega.alpha = omega_alpha.save, 
                            omega.beta = omega_beta.save,
                            dev = LL.save)
    
    standardize.param = rbind(data.frame(Type = "X", 
                                         Name = "AOD/CTM", 
                                         Mean = X.mean, 
                                         SD = X.sd),
                               data.frame(Type = "L", 
                                          Name = colnames(L), 
                                          Mean = L.mean, 
                                          SD = L.sd), 
                               data.frame(Type = "M", 
                                          Name = colnames(M), 
                                          Mean = M.mean, 
                                          SD = M.sd))
    row.names(standardize.param) = NULL

    nngp.info.save <- NULL
    if (nngp) {
      nngp.info.save <- list(ordered.coords = ordered_coords,
                             coord.ordering = coord_ordering,
                             coord.reverse.ordering = coord_reverse_ordering,
                             neighbors = neighbors,
                             dist.matrices = dist_matrices,
                             space.to.spacetime.assign = space_to_spacetime_assign,
                             num_neighbors = num_neighbors)
    }

    discrete.theta.alpha.info.save <- NULL
    if (!is.null(discrete.theta.alpha.values)) {
        discrete.theta.alpha.info.save <- list(which.theta.alpha.discrete = which.theta.alpha.discrete,
                                               kernals.inv.alpha = kernals_inv_alpha)

    }

    discrete.theta.beta.info.save <- NULL
    if (!is.null(discrete.theta.beta.values)) {
        discrete.theta.beta.info.save <- list(which.theta.beta.discrete = which.theta.beta.discrete,
                                               kernals.inv.beta = kernals_inv_beta)
    }
    
    list(delta = delta.save, 
         gamma = gamma.save,
         alpha.time = alpha_time.save, 
         beta.time = beta_time.save,
         alpha.space = alpha_space.save, 
         beta.space = beta_space.save,
         others = other.save, 
         Y = Y.hat, 
         standardize.param = standardize.param,
         theta.acc = theta.acc,
         tau.acc = tau.acc,
         nngp.info = nngp.info.save,
         discrete.theta.alpha.info = discrete.theta.alpha.info.save,
         discrete.theta.beta.info = discrete.theta.beta.info.save,
         cov_kern = cov_kern)
}

