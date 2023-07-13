#' Fit the Geostatistical Regression Model (GRM) With Cross Validation
#'
#' This function fits Bayesian Hierarchical Model (BHM) in the form of Y ~ beta X + gamma L + delta M with cross-validation
#'
#' @inheritParams grm_pred
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
                  nngp = F,
                  num_neighbors = 10,
                  discrete.theta.alpha.values = NULL,
                  discrete.theta.beta.values = NULL,
                  discrete.theta.gibbs = T,
                  n.iter = 25000,
                  burn = 5000,
                  thin = 4,
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

    cv.id <- cv.object$cv.id

    Y.cv = data.frame(time_id = time.id, 
                      space_id = space.id, 
                      obs = Y, 
                      estimate = NA, 
                      sd = NA)
  
  
    for (cv.i in 1:cv.object$num.folds) {
    
        print(paste0("Performing CV Experiment ---- Fold ", cv.i))
        train.id <- cv.id != cv.i & cv.id != 0
        test.id.temp <- cv.id == cv.i
        #remove any test observations that are not within the training observations time range
        #these will be NA's in final cv predictions
        test.id.remove <- (min(time.id[train.id]) > time.id | max(time.id[train.id]) < time.id) & test.id.temp
        test.id <- test.id.temp & !test.id.remove

        time.id.train = time.id[train.id]
        time.id.test = time.id[test.id]

        Y.train = Y[train.id]
        Y.test = Y[test.id]

        X.train = X[train.id]
        X.test = X[test.id]
    
        #Subset of L matrix based on variable s
        L = as.matrix(L)
        L.train = L[train.id, , drop = FALSE]
        L.test = L[test.id, , drop = FALSE]
        M = as.matrix(M)
        M.train = M[train.id, , drop = FALSE]
        M.test = M[test.id, , drop = FALSE]
    
        space.id.train = space.id[train.id]
        space.id.test = space.id[test.id]

        #grm requires space.id to be from 1:max(space_id)
        #spatial cross validation breaks this assumption (missing space_id values)
        #here we create temporary space.id values that are from 1:length(unique(space.id))
        space.id.train.key = sort(unique(space.id.train))
        space.id.test.key = sort(unique(space.id.test))
        space.id.train.temp = sapply(space.id.train, 
                                     function(x) which(space.id.train.key == x))
        space.id.test.temp = sapply(space.id.test,
                                    function(x) which(space.id.test.key == x))
        spacetime.id.train = spacetime.id[train.id]
        spacetime.id.test = spacetime.id[test.id]
        coords.train = coords[train.id, ]
        coords.test = coords[test.id, ]
   
        fit.cv = grm(Y = Y.train, 
                     X = X.train, 
                     L = L.train, 
                     M = M.train, 
                     coords = coords.train,
                     space.id = space.id.train.temp, 
                     time.id = time.id.train, 
                     spacetime.id = spacetime.id.train, 
                     include.additive.weekly.resid = include.additive.weekly.resid,
                     include.additive.annual.resid = include.additive.annual.resid,
                     include.multiplicative.weekly.resid = include.multiplicative.weekly.resid,
                     include.multiplicative.annual.resid = include.multiplicative.annual.resid,
                     nngp = nngp,
                     num_neighbors = num_neighbors,
                     discrete.theta.alpha.values = discrete.theta.alpha.values,
                     discrete.theta.beta.values = discrete.theta.beta.values,
                     discrete.theta.gibbs = discrete.theta.gibbs,
                     n.iter = n.iter,
                     burn = burn,
                     thin = thin,
                     tau.alpha.tune = tau.alpha.tune,
                     tau.alpha.a = tau.alpha.a,
                     tau.alpha.b = tau.alpha.b,
                     tau.beta.tune = tau.beta.tune,
                     tau.beta.a = tau.beta.a,
                     tau.beta.b = tau.beta.b,
                     omega.a = omega.a,
                     omega.b = omega.b,
                     theta.alpha.tune = theta.alpha.tune, 
                     theta.alpha.a = theta.alpha.a, 
                     theta.alpha.b = theta.alpha.b,
                     theta.alpha.init = theta.alpha.init,
                     theta.beta.tune = theta.beta.tune, 
                     theta.beta.a = theta.beta.a, 
                     theta.beta.b = theta.beta.b,
                     theta.beta.init = theta.beta.init,
                     rho.alpha.init = rho.alpha.init,
                     rho.beta.init = rho.beta.init,
                     sigma.a = sigma.a, 
                     sigma.b = sigma.b,
                     verbose = verbose,
                     verbose.iter = verbose.iter)
    
        in_sample <- if (cv.object$type == "ordinary") {

            TRUE

        } else if (cv.object$type %in% c("spatial", "spatial_clustered")) {

            FALSE

        } else {

            stop("cv.object$type must be either 'ordinary', 'spatial', or 'spatial_clustered'")

        }


        cv.results = grm_pred(grm.fit = fit.cv, 
                              n.iter = (n.iter - burn) / thin,
                              X.pred = X.test, 
                              L.pred = L.test, 
                              M.pred = M.test, 
                              coords.Y = coords.train,
                              coords.pred = coords.test,
                              space.id.Y = space.id.train.temp,
                              space.id = space.id.test.temp,
                              time.id = time.id.test,
                              spacetime.id = spacetime.id.test,
                              in.sample = in_sample)



        Y.cv$estimate[test.id] = cv.results$estimate
        Y.cv$sd[test.id] = cv.results$sd

    }
 
    Y.cv$upper_95 = Y.cv$estimate + 1.96 * Y.cv$sd
    Y.cv$lower_95 = Y.cv$estimate - 1.96 * Y.cv$sd
  
    return(Y.cv)
}
