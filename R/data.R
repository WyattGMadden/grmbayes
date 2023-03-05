#' Chemical Transport Model and PM2.5 Data For Atlanta Metro Area, 2003-2005
#' PM2.5 Monitor and Corresponding Chemical Transport Model Data For Atlanta Metro Area, 2003-2005
#'
#' Daily PM2.5 readings, CMT output, and covariates of interest in Atlanta Metro Area from 2003-2005. 
#' Each row of the data is a unique monitor location/day. 
#' Data included only if CMT output is available for that day and monitoring location.
#' (ie. if the monitor is within a 12km x 12km grid square for which CMT data is available).
#'
#' @format ## `ctm_pm25`
#' A data frame with 11,404 rows and 16 columns:
#' \describe{
#'   \item{time_id}{Week identifier}
#'   \item{space_id}{Location identifier}
#'   \item{spacetime_id}{Year identifier}
#'   \item{pm25}{PM2.5 reading}
#'   \item{ctm}{Chemical transport model value}
#'   \item{elevation}{Elevation spatial covariate (US Geological Survey)}
#'   \item{forestcover}{Percentage of forest cover spatial covariate (2001 National Land Cover database)}
#'   \item{hwy_length}{Sum of major roadway length spatial covariate (2001 National Land Cover database)}
#'   \item{lim_hwy_length}{Sum of major roadway length spatial covariate (2001 National Land Cover database)}
#'   \item{lim_hwy_length}{Sum of local roadway lengths spatial covariate (2001 National Land Cover database)}
#'   \item{point_emi_any}{Indicator of PM2.5 primary emission point source spatial covariate (2002 USEPA National Emissions Inventory)}
#'   \item{tmp}{Temperature spatio-temporal covariate (North American Land Data Assimilation Systems)}
#'   \item{wind}{Wind spatio-temporal covariate (North American Land Data Assimilation Systems)}
#'   \item{date}{Date of measurement}
#'   \item{x}{Location x-coordinate}
#'   \item{y}{Location y-coordinate}
#'   ...
#' }
#' @source <https://www.nature.com/articles/jes201390>
"ctm_pm25"


#' PM2.5 Monitor and Corresponding Satellite Aerial Optical Depth Data For Atlanta Metro Area, 2003-2005
#'
#' Daily PM2.5 readings, AOD, and covariates of interest in Atlanta Metro Area from 2003-2005. 
#' Each row of the data is a unique location/day.
#' Data included only if AOD output is available for that day and monitoring location.
#' (ie. if the monitor is within a 12km x 12km grid square for which AOD data is available).
#'
#' @format ## `maia_pm25`
#' A data frame with 6531 rows and 20 columns:
#' \describe{
#'   \item{time_id}{Week identifier}
#'   \item{space_id}{Location identifier}
#'   \item{spacetime_id}{Year identifier}
#'   \item{pm25}{PM2.5 reading}
#'   \item{aod}{Aerial optical depth value}
#'   \item{elevation}{Elevation spatial covariate (US Geological Survey)}
#'   \item{forestcover}{Percentage of forest cover spatial covariate (2001 National Land Cover database)}
#'   \item{hwy_length}{Sum of major roadway length spatial covariate (2001 National Land Cover database)}
#'   \item{lim_hwy_length}{Sum of major roadway length spatial covariate (2001 National Land Cover database)}
#'   \item{lim_hwy_length}{Sum of local roadway lengths spatial covariate (2001 National Land Cover database)}
#'   \item{point_emi_any}{Indicator of PM2.5 primary emission point source spatial covariate (2002 USEPA National Emissions Inventory)}
#'   \item{tmp}{Temperature spatio-temporal covariate (North American Land Data Assimilation Systems)}
#'   \item{wind}{Wind spatio-temporal covariate (North American Land Data Assimilation Systems)}
#'   \item{cmaq}{CMT value spatio-temporal covariate (USEPA Models-3/Community Multiscale Air Quality)}
#'   \item{date}{Date of measurement}
#'   \item{x}{Location x-coordinate}
#'   \item{y}{Location y-coordinate}
#'   ...
#' }
#' @source <https://pubmed.ncbi.nlm.nih.gov/31465992/>
"maia_pm25"

#' Chemical Transport Model Data For Atlanta Metro Area, June 2004
#'
#' Daily CMT output, and covariates of interest in Atlanta Metro Area for June 2004. 
#' Data included for every 12km x 12km grid square in study area. 
#' Each row of the data is a unique location/day. 
#'
#' @format ## `ctm_preds`
#' A data frame with 72,000 rows and 15 columns:
#' \describe{
#'   \item{time_id}{Week identifier}
#'   \item{space_id}{Location identifier}
#'   \item{spacetime_id}{Year identifier}
#'   \item{ctm}{Chemical transport model value}
#'   \item{elevation}{Elevation spatial covariate (US Geological Survey)}
#'   \item{forestcover}{Percentage of forest cover spatial covariate (2001 National Land Cover database)}
#'   \item{hwy_length}{Sum of major roadway length spatial covariate (2001 National Land Cover database)}
#'   \item{lim_hwy_length}{Sum of major roadway length spatial covariate (2001 National Land Cover database)}
#'   \item{lim_hwy_length}{Sum of local roadway lengths spatial covariate (2001 National Land Cover database)}
#'   \item{point_emi_any}{Indicator of PM2.5 primary emission point source spatial covariate (2002 USEPA National Emissions Inventory)}
#'   \item{tmp}{Temperature spatio-temporal covariate (North American Land Data Assimilation Systems)}
#'   \item{wind}{Wind spatio-temporal covariate (North American Land Data Assimilation Systems)}
#'   \item{date}{Date of measurement}
#'   \item{x}{Location x-coordinate}
#'   \item{y}{Location y-coordinate}
#'   ...
#' }
#' @source <https://www.nature.com/articles/jes201390>
"ctm_preds"


#' Satellite Aerial Optical Depth and PM2.5 Data For Atlanta Metro Area, 2003-2005
#'
#' Daily AOD, and covariates of interest in Atlanta Metro Area for June 2004.
#' Data included for every 12km x 12km grid square in study area for which AOD data is available. 
#' Each row of the data is a unique location/day. 
#'
#' @format ## `maia_pm25`
#' A data frame with 6531 rows and 20 columns:
#' \describe{
#'   \item{time_id}{Week identifier}
#'   \item{space_id}{Location identifier}
#'   \item{spacetime_id}{Year identifier}
#'   \item{aod}{Aerial optical depth value}
#'   \item{elevation}{Elevation spatial covariate (US Geological Survey)}
#'   \item{forestcover}{Percentage of forest cover spatial covariate (2001 National Land Cover database)}
#'   \item{hwy_length}{Sum of major roadway length spatial covariate (2001 National Land Cover database)}
#'   \item{lim_hwy_length}{Sum of major roadway length spatial covariate (2001 National Land Cover database)}
#'   \item{lim_hwy_length}{Sum of local roadway lengths spatial covariate (2001 National Land Cover database)}
#'   \item{point_emi_any}{Indicator of PM2.5 primary emission point source spatial covariate (2002 USEPA National Emissions Inventory)}
#'   \item{tmp}{Temperature spatio-temporal covariate (North American Land Data Assimilation Systems)}
#'   \item{wind}{Wind spatio-temporal covariate (North American Land Data Assimilation Systems)}
#'   \item{cmaq}{CMT value spatio-temporal covariate (USEPA Models-3/Community Multiscale Air Quality)}
#'   \item{date}{Date of measurement}
#'   \item{x}{Location x-coordinate}
#'   \item{y}{Location y-coordinate}
#'   ...
#' }
#' @source <https://pubmed.ncbi.nlm.nih.gov/31465992/>
"maia_preds"

