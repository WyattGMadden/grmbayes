#' Fit the Geostatistical Regression Model (GRM) With Cross Validation
#'
#' This function fits Bayesian Hierarchical Model (BHM) in the form of Y ~ beta X + gamma L + delta M with cross-validation
#'
#' @inheritParams grm
#' @param cv.object A named list containing cv.id, num.folds, and type. Can be created with create_cv function. 
#'
#' @return A data frame containing cross validation predictions
#'
#' @examples
#' # grm_cv()
#' 
#' 
#' @export
grm_cv = function(Y, 
                  X, 
                  cv.object,
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
                  n.iter = 25000,
                  burn = 5000,
                  thin = 4,
                  tau.a = 0.5,
                  tau.b = 0.005,
                  omega.a = 0.5,
                  omega.b = 0.005,
                  theta.tune = 0.2, 
                  theta.a = 5, 
                  theta.b = 0.05,
                  sigma.a = 0.001, 
                  sigma.b = 0.001,
                  verbose = TRUE) {

    CV.id <- cv.object$cv.id
    num.folds <- cv.object$num.folds

    Y.CV = data.frame(time_id = time.id, 
                      space_id = space.id, 
                      obs = Y, 
                      estimate = NA, 
                      sd = NA)
  
  
    for (CV.i in 1:num.folds) {
    
        print(paste0("Performing CV Experiment ---- Fold ", CV.i))
        Y.train = Y[CV.id != CV.i]
        Y.test = Y[CV.id == CV.i]
        X.train = X[CV.id != CV.i]
        X.test = X[CV.id == CV.i]
    
        #Subset of L matrix based on variable s
        L = as.matrix(L)
        L.train = L[CV.id != CV.i, , drop = FALSE]
        L.test = L[CV.id == CV.i, , drop = FALSE]
        M = as.matrix(M)
        M.train = M[CV.id != CV.i, , drop = FALSE]
        M.test = M[CV.id == CV.i, , drop = FALSE]
    
        Time_ID.train = time.id[CV.id != CV.i]
        Time_ID.test = time.id[CV.id == CV.i]
        Space_ID.train = space.id[CV.id != CV.i]
        Space_ID.test = space.id[CV.id == CV.i]
        SpaceTime_ID.train = spacetime.id[CV.id != CV.i]
        SpaceTime_ID.test = spacetime.id[CV.id == CV.i]
        coords.train = coords[CV.id != CV.i, ]
    
        fit.cv = grm(Y = Y.train, 
                     X = X.train, 
                     L = L.train, 
                     M = M.train, 
                     coords = coords.train,
                     space.id = Space_ID.train, 
                     time.id = Time_ID.train, 
                     spacetime.id = SpaceTime_ID.train, 
                     include.additive.weekly.resid = include.additive.weekly.resid,
                     include.additive.annual.resid = include.additive.annual.resid,
                     include.multiplicative.weekly.resid = include.multiplicative.weekly.resid,
                     include.multiplicative.annual.resid = include.multiplicative.annual.resid,
                     n.iter = n.iter,
                     burn = burn,
                     thin = thin,
                     tau.a = tau.a,
                     tau.b = tau.b,
                     omega.a = omega.a,
                     omega.b = omega.b,
                     theta.tune = theta.tune, 
                     theta.a = theta.a, 
                     theta.b = theta.b,
                     sigma.a = sigma.a, 
                     sigma.b = sigma.b,
                     verbose = verbose)
    
        #Standardize
        standardize.param = fit.cv$standardize.param
    
        X.test = (X.test - standardize.param[standardize.param$Type == "X", ]$Mean) / 
            standardize.param[standardize.param$Type == "X", ]$SD
    
        L.var = as.character(standardize.param[standardize.param$Type == "L", ]$Name)

        for (l in L.var) {
            L.test[, colnames(L.test) == l] = (L.test[, colnames(L.test) == l] - 
                                               standardize.param[standardize.param$Name == l, ]$Mean) / 
            standardize.param[standardize.param$Name == l, ]$SD
    }
    
    M.var = as.character(standardize.param[standardize.param$Type == "M", ]$Name)

    for (m in M.var) {
        M.test[, colnames(M.test) == m] = (M.test[, colnames(M.test) == m] - 
                                           standardize.param[standardize.param$Name == m, ]$Mean) / 
        standardize.param[standardize.param$Name == m, ]$SD
    }
    ##Make predictions
    CV.Results = data.frame(Time_ID = time.id[CV.id == CV.i],
                            Space_ID = space.id[CV.id == CV.i],
                            SpaceTime_ID = spacetime.id[CV.id == CV.i])
    CV.Results$Estimate = 0
    CV.Results$SD = 0
    
    id.temp.train = paste0(fit.cv$alpha.space$space.id, 
                           "_", 
                           fit.cv$alpha.space$spacetime.id)

    id.temp = paste0(Space_ID.test, 
                     "_", 
                     SpaceTime_ID.test)
    
    others = fit.cv$others
    delta = as.matrix(fit.cv$delta)
    gamma = as.matrix(fit.cv$gamma)
    
    for (m in 1:nrow(delta)) {
        intercept = others$alpha0[m] + 
            fit.cv$alpha.time[Time_ID.test, m + 1] + 
            fit.cv$alpha.space[match(id.temp, id.temp.train), m + 2]
        slope = others$beta0[m] + 
            fit.cv$beta.time[Time_ID.test, m+1] + 
            fit.cv$beta.space[match(id.temp, id.temp.train), m + 2]
        fix.L = L.test %*% gamma[m, ]
        fix.M = M.test %*% delta[m, ]
      
        pred.mu = intercept + slope * X.test + fix.L + fix.M
        pred.mu = pred.mu + stats::rnorm(length(pred.mu), 0, sqrt(others$sigma2[m])) 
      
        CV.Results$Estimate = CV.Results$Estimate + pred.mu / nrow(delta)
        CV.Results$SD = CV.Results$SD + pred.mu ^ 2 / nrow(delta)
    }

    CV.Results$SD = sqrt((CV.Results$SD - CV.Results$Estimate ^ 2))

    
    Y.CV$estimate[CV.id == CV.i] = CV.Results$Estimate
    Y.CV$sd[CV.id == CV.i] = CV.Results$SD

    }
 
    Y.CV$upper_95 = Y.CV$estimate + 1.96 * Y.CV$sd
    Y.CV$lower_95 = Y.CV$estimate - 1.96 * Y.CV$sd
  
    return(Y.CV)
}
