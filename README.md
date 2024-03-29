
<!-- README.md is generated from README.Rmd. Please edit that file -->

# grmbayes

<!-- badges: start -->

[![R-CMD-check](https://github.com/WyattGMadden/ensembleDownscaleR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/WyattGMadden/ensembleDownscaleR/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The goal of `grmbayes` is to provide functions for fitting efficient
Bayesian geostatistical regression models.

## Installation

You can install the development version of ensembleDownscaleR like so:

``` r
devtools::install_github("WyattGMadden/grmbayes", build_vignettes = TRUE)
```

## Example

### Stage 1

``` r
library(grmbayes)

?cmaq_aqs_matched

cmaq_fit <- grm(Y = cmaq_aqs_matched$pm25,
                X = cmaq_aqs_matched$ctm,
                L = cmaq_aqs_matched[, c("elevation", "forestcover",
                                 "hwy_length", "lim_hwy_length", 
                                 "local_rd_length", "point_emi_any")],
                M = cmaq_aqs_matched[, c("tmp", "wind")],
                coords = cmaq_aqs_matched[, c("x", "y")],
                n.iter = 100,
                burn = 20,
                thin = 4,
                nngp = T,
                num_neighbors = 10,
                covariance = "matern",
                matern.nu = 0.5,
                space.id = cmaq_aqs_matched$space_id,
                time.id = cmaq_aqs_matched$time_id,
                spacetime.id = cmaq_aqs_matched$spacetime_id,
                verbose.iter = 10)


cv_id_ctm_ord <- create_cv(space.id = cmaq_aqs_matched$space_id,
                           time.id = cmaq_aqs_matched$time_id, 
                           type = "ordinary")



cmaq_fit_cv <- grm_cv(Y = cmaq_aqs_matched$pm25,
                      X = cmaq_aqs_matched$ctm,
                      cv.object = cv_id_ctm_ord,
                      L = cmaq_aqs_matched[, c("elevation", "forestcover",
                                       "hwy_length", "lim_hwy_length", 
                                       "local_rd_length", "point_emi_any")],
                      M = cmaq_aqs_matched[, c("tmp", "wind")],
                      n.iter = 1000,
                      burn = 200,
                      thin = 4,
                      coords = cmaq_aqs_matched[, c("x", "y")],
                      nngp = T,
                      num_neighbors = 10,
                      covariance = "matern",
                      matern.nu = 0.5,
                      space.id = cmaq_aqs_matched$space_id,
                      time.id = cmaq_aqs_matched$time_id,
                      spacetime.id = cmaq_aqs_matched$spacetime_id,
                      verbose.iter = 10)


?modis_aqs_matched

cv_id_modis_ord <- create_cv(space.id = modis_aqs_matched$space_id,
                             time.id = modis_aqs_matched$time_id,
                             type = "ordinary")

modis_fit <- grm(Y = modis_aqs_matched$pm25,
                 X = modis_aqs_matched$aod,
                 L = modis_aqs_matched[, c("elevation", "forestcover",
                                   "hwy_length", "lim_hwy_length", 
                                   "local_rd_length", "point_emi_any")],
                 M = modis_aqs_matched[, c("tmp", "wind", "cmaq", "tempaod", 
                                   "windaod", "elevationaod")],
                 n.iter = 500,
                 num_neighbors = 6,
                 burn = 100,
                 thin = 4,
                 coords = modis_aqs_matched[, c("x", "y")],
                 space.id = modis_aqs_matched$space_id,
                 time.id = modis_aqs_matched$time_id,
                 spacetime.id = modis_aqs_matched$spacetime_id)

modis_fit_cv <- grm_cv(Y = modis_aqs_matched$pm25,
                       X = modis_aqs_matched$aod,
                       cv.object = cv_id_maia_ord,
                       L = modis_aqs_matched[, c("elevation", "forestcover",
                                         "hwy_length", "lim_hwy_length", 
                                         "local_rd_length", "point_emi_any")],
                       M = modis_aqs_matched[, c("tmp", "wind", "cmaq", "tempaod", 
                                         "windaod", "elevationaod")],
                       n.iter = 500,
                       burn = 100,
                       thin = 4,
                       coords = modis_aqs_matched[, c("x", "y")],
                       space.id = modis_aqs_matched$space_id,
                       time.id = modis_aqs_matched$time_id,
                       spacetime.id = modis_aqs_matched$spacetime_id)
```

### Stage 2

``` r
?cmaq_full
cmaq_pred <- grm_pred(grm.fit = cmaq_fit,
                      X.pred = cmaq_full$ctm,
                      L.pred = cmaq_full[, c("elevation", "forestcover",
                                             "hwy_length", "lim_hwy_length", 
                                             "local_rd_length", "point_emi_any")],
                      M.pred = cmaq_full[, c("tmp", "wind")],
                      coords.Y = cmaq_aqs_matched[, c("x", "y")],
                      space.id.Y = cmaq_aqs_matched$space_id,
                      coords.pred = cmaq_full[, c("x", "y")],
                      space.id = cmaq_full$space_id,
                      time.id = cmaq_full$time_id,
                      spacetime.id = cmaq_full$spacetime_id,
                      n.iter = 20,
                      verbose = T,
                      include.random.effects = T)



?modis_full

modis_pred <- grm_pred(grm.fit = modis_fit,
                       X.pred = modis_full$aod, 
                       L.pred = modis_full[, c("elevation", "forestcover",
                                               "hwy_length", "lim_hwy_length", 
                                               "local_rd_length", "point_emi_any")],
                       M.pred = modis_full[, c("tmp", "wind", "cmaq", "tempaod", 
                                               "windaod", "elevationaod")],
                       coords.Y = modis_aqs_matched[, c("x", "y")],
                       space.id.Y = modis_aqs_matched$space_id,
                       coords.pred = modis_full[, c("x", "y")],
                       space.id = modis_full$space_id, 
                       time.id = modis_full$time_id, 
                       spacetime.id = modis_full$spacetime_id,
                       n.iter = 100,
                       verbose = T)
```
