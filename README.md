
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ensembleDownscaleR

<!-- badges: start -->

[![R-CMD-check](https://github.com/WyattGMadden/ensembleDownscaleR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/WyattGMadden/ensembleDownscaleR/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The goal of ensembleDownscaleR is to …

## Installation

You can install the development version of ensembleDownscaleR like so:

``` r
devtools::install_github("WyattGMadden/ensembleDownscaleR")
```

## Example

### Stage 1

``` r
library(ensembleDownscaleR)

?ctm_pm25

ctm_fit <- grm(Y = ctm_pm25$pm25,
               X = ctm_pm25$ctm,
               L = ctm_pm25[, c("elevation", "forestcover",
                                "hwy_length", "lim_hwy_length", 
                                "local_rd_length", "point_emi_any")],
               M = ctm_pm25[, c("tmp", "wind")],
               n.iter = 500,
               burn = 100,
               thin = 4,
               coords = ctm_pm25[, c("x", "y")],
               space.id = ctm_pm25$space_id,
               time.id = ctm_pm25$time_id,
               spacetime.id = ctm_pm25$spacetime_id)

ctm_fit_cv <- grm_cv(Y = ctm_pm25$pm25,
                     X = ctm_pm25$pm25,
                     L = ctm_pm25[, c("elevation", "forestcover",
                                      "hwy_length", "lim_hwy_length", 
                                      "local_rd_length", "point_emi_any")],
                     M = ctm_pm25[, c("tmp", "wind")],
                     n.iter = 500,
                     burn = 100,
                     thin = 4,
                     coords = ctm_pm25[, c("x", "y")],
                     space.id = ctm_pm25$space_id,
                     time.id = ctm_pm25$time_id,
                     spacetime.id = ctm_pm25$spacetime_id)


?maia_pm25


maia_fit <- grm(Y = maia_pm25$pm25,
                X = maia_pm25$aod,
                L = maia_pm25[, c("elevation", "forestcover",
                                  "hwy_length", "lim_hwy_length", 
                                  "local_rd_length", "point_emi_any")],
                M = maia_pm25[, c("tmp", "wind", "cmaq", "tempaod", 
                                  "windaod", "elevationaod")],
                n.iter = 500,
                burn = 100,
                thin = 4,
                coords = maia_pm25[, c("x", "y")],
                space.id = maia_pm25$space_id,
                time.id = maia_pm25$time_id,
                spacetime.id = maia_pm25$spacetime_id)

maia_fit_cv <- grm_cv(Y = maia_pm25$pm25,
                      X = maia_pm25$aod,
                      L = maia_pm25[, c("elevation", "forestcover",
                                        "hwy_length", "lim_hwy_length", 
                                        "local_rd_length", "point_emi_any")],
                      M = maia_pm25[, c("tmp", "wind", "cmaq", "tempaod", 
                                        "windaod", "elevationaod")],
                      n.iter = 500,
                      burn = 100,
                      thin = 4,
                      coords = maia_pm25[, c("x", "y")],
                      space.id = maia_pm25$space_id,
                      time.id = maia_pm25$time_id,
                      spacetime.id = maia_pm25$spacetime_id)
```

### Stage 2

``` r
?ctm_preds

ctm_pred <- grm_pred(grm.fit = ctm_fit,
                     X.pred = ctm_preds$ctm,
                     L.pred = ctm_preds[, c("elevation", "forestcover",
                                            "hwy_length", "lim_hwy_length", 
                                            "local_rd_length", "point_emi_any")],
                     M.pred = ctm_preds[, c("tmp", "wind")],
                     coords.Y = ctm_pm25[, c("x", "y")],
                     space.id.Y = ctm_pm25$space_id,
                     coords.pred = ctm_preds[, c("x", "y")],
                     space.id = ctm_preds$space_id,
                     time.id = ctm_preds$time_id,
                     spacetime.id = ctm_preds$spacetime_id,
                     include.additive.annual.resid = T,
                     include.multiplicative.annual.resid = T,
                     n.iter = 100,
                     verbose = T)



?maia_preds

maia_pred <- grm_pred(grm.fit = maia_fit,
                      X.pred = maia_preds$aod, 
                      L.pred = maia_preds[, c("elevation", "forestcover",
                                              "hwy_length", "lim_hwy_length", 
                                              "local_rd_length", "point_emi_any")],
                      M.pred = maia_preds[, c("tmp", "wind", "cmaq", "tempaod", 
                                              "windaod", "elevationaod")],
                      coords.Y = maia_pm25[, c("x", "y")],
                      space.id.Y = maia_pm25$space_id,
                      coords.pred = maia_preds[, c("x", "y")],
                      space.id = maia_preds$space_id, 
                      time.id = maia_preds$time_id, 
                      spacetime.id = maia_preds$spacetime_id,
                      include.additive.annual.resid = T,
                      include.multiplicative.annual.resid = T,
                      n.iter = 100,
                      verbose = T)
```

### Stage 3

``` r

ensemble_fit <- ensemble_spatial(grm.fit.cv.1 = ctm_fit_cv,
                                 grm.fit.cv.2 = maia_fit_cv,
                                 date.Y.1 = ctm_pm25$date,
                                 date.Y.2 = maia_pm25$date,
                                 coords.Y.1 = ctm_pm25[, c("x", "y")],
                                 space.id.Y.1 = ctm_pm25$space_id,
                                 n.iter = 5000, 
                                 burn = 1000, 
                                 thin = 4,
                                 tau.a = 0.001, 
                                 tau.b = 0.001, 
                                 theta.tune = 0.2, 
                                 theta.a = 5, 
                                 theta.b = 0.05)
```

### Stage 4

``` r

weight_preds <- weight_pred(ensemble.fit = ensemble_fit,
                            coords.Y.1 = ctm_pm25[, c("x", "y")], 
                            space.id.Y.1 = ctm_pm25$space_id, 
                            coords.pred.1 = ctm_preds[, c("x", "y")],
                            space.id.pred.1 = ctm_preds$space_id)


results <- gap_fill(grm.pred.1 = ctm_pred,
                    grm.pred.2 = maia_pred,
                    date.pred.1 = ctm_preds$date,
                    date.pred.2 = maia_preds$date, 
                    space.id = ctm_pred$space.id, 
                    weights = weight_preds)
```
