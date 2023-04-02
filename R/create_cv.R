
#' Create Cross Validate ID's
#'
#' This function creates a vector of cross validation ID's for a given number of folds and type of cross validation. 
#'
#' @inheritParams grm
#' @param num.folds Number of folds used in the cross validation process (default = 10)
#' @param type Type of cross validation to be performed. Options are "spatial", "ordinary", or "spatial_clustered". (default = "spatial")
#'
#' @return A named list containing a vector of cross validation ID's, number of folds, and cross validation type. 
#'
#' @examples
#' # grm_cv()
#' 
#' 
#' @export
create_cv <- function(space.id,
                      time.id,
                      num.folds = 10,
                      type = "spatial",
                      coords = NULL) {

    #remove first/last time observations
    time_obs_1 <- sum(time.id == 1)
    time_obs_max <- sum(time.id == max(time.id))

    space_id_cv <- space.id[(time_obs_1 + 1):(length(space.id) - time_obs_max)]
    cv_id <- rep(NA, length(space.id) - time_obs_1 - time_obs_max)


    if (type == "spatial") {

        #unshuffled cv spatial id's with (almost) equal number of observations in each fold
        cv_spat_id <- (1:num.folds)[(1:max(space_id_cv)) %% num.folds + 1]

        #shuffle cv spatial id's
        cv_spat_id <- sample(cv_spat_id, replace = F)
        
        #assign spatial id's to cv spatial id's
        cv_id <- sapply(space_id_cv, function(i) cv_spat_id[i])

    } else if (type == "ordinary") {

        for (i in 1:max(space_id_cv)) {

            # make sure no spatial location is not in some but not all folds
            if (sum(space_id_cv == i) < num.folds) {

                cv_id[space_id_cv == i] <- rep(0, sum(space_id_cv == 1))

            } else {


                #number of observations in site i for cv
                obs_for_site_i <- sum(space_id_cv == i) 

                #unshuffled cv id's with (almost) equal number of observations in each fold
                cv_id_i <- (1:num.folds)[(1:obs_for_site_i) %% num.folds + 1]

                #shuffle cv id's
                cv_id_i = sample(cv_id_i, replace = F)
            
                cv_id[space_id_cv == i] <- cv_id_i
            }
        
        }
            

    } else if (type == "spatial_clustered") {

        if (is.null(coords)) {
            stop("coords must be provided for spatial_clustered cross validation")
        }

        
        #remove first/last time observations
        coords <- coords[(time_obs_1 + 1):(nrow(coords) - time_obs_max), ]
        cv_id <- stats::kmeans(x = coords, 
                               centers = num.folds)$cluster |>
            unname()



    } else {
        stop("type must be either 'spatial', 'ordinary', or 'spatial_clustered'")
    }

    cv_id_full <- c(rep(0, time_obs_1), cv_id, rep(0, time_obs_max))

    return(list(cv_id = cv_id_full, 
                num.folds = num.folds, 
                type = type))
}
