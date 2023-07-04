
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
               n.iter = 1000,
               burn = 100,
               thin = 4,
               nngp = F,
               #discrete.theta.beta.values = seq(1, 150, 1),
#               num_neighbors = 6,
               coords = ctm_pm25[, c("x", "y")],
               space.id = ctm_pm25$space_id,
               time.id = ctm_pm25$time_id,
               spacetime.id = ctm_pm25$spacetime_id,
               verbose.iter = 10)

summary(ctm_fit$others$theta.beta[100:225])

hist(ctm_fit$others$theta.beta, breaks = 100)
summary(ctm_fit$others$theta.beta)

test <- list(matrix(1, ncol = 3, nrow = 3), 
     matrix(2, ncol = 3, nrow = 3),
     matrix(3, ncol = 3, nrow = 3)) 

test 
lapply(test, function(x) x + x)

ctm_fit_dis <- grm(Y = ctm_pm25$pm25,
               X = ctm_pm25$ctm,
               L = ctm_pm25[, c("elevation", "forestcover",
                                "hwy_length", "lim_hwy_length", 
                                "local_rd_length", "point_emi_any")],
               M = ctm_pm25[, c("tmp", "wind")],
               n.iter = 10000,
               burn = 2000,
               thin = 4,
               nngp = F,
               discrete.theta.beta.values = seq(1, 20, 1),
               discrete.theta.gibbs = T,
#               num_neighbors = 6,
               coords = ctm_pm25[, c("x", "y")],
               space.id = ctm_pm25$space_id,
               time.id = ctm_pm25$time_id,
               spacetime.id = ctm_pm25$spacetime_id,
               verbose.iter = 10)


str(ctm_fit_dis$others$theta.beta)
hist(ctm_fit_dis$others$theta.beta, breaks = 100)
summary(ctm_fit_dis$others$theta.beta)
tibble(x = 1:length(ctm_fit_dis$others$theta.beta), 
       y = ctm_fit_dis$others$theta.beta) %>%
    ggplot(aes(x = x, y = y)) +
    geom_line()
library(tidyverse)
tibble(x = 1:length(ctm_fit$others$theta.beta), 
       y = ctm_fit$others$theta.beta) %>%
    ggplot(aes(x = x, y = y)) +
    geom_line()
ctm_fit$theta.acc

cv_id_ctm_spat <- create_cv(space.id = ctm_pm25$space_id,
                            time.id = ctm_pm25$time_id)

cv_id_ctm_ord <- create_cv(space.id = ctm_pm25$space_id,
                           time.id = ctm_pm25$time_id, type = "ordinary")

cv_id_ctm_spcl <- create_cv(space.id = ctm_pm25$space_id,
                            time.id = ctm_pm25$time_id,
                            type = "spatial_clustered",
                            coords = ctm_pm25[, c("x", "y")])



ctm_fit_cv <- grm_cv(Y = ctm_pm25$pm25,
                     X = ctm_pm25$ctm,
                     cv.object = cv_id_ctm_ord,
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


ctm_fit_cv_nngp <- grm_cv(Y = ctm_pm25$pm25,
                          X = ctm_pm25$ctm,
                          cv.object = cv_id_ctm_ord,
                          L = ctm_pm25[, c("elevation", "forestcover",
                                           "hwy_length", "lim_hwy_length", 
                                           "local_rd_length", "point_emi_any")],
                          M = ctm_pm25[, c("tmp", "wind")],
                          n.iter = 500,
                          burn = 100,
                          thin = 4,
                          nngp = T,
                          coords = ctm_pm25[, c("x", "y")],
                          space.id = ctm_pm25$space_id,
                          time.id = ctm_pm25$time_id,
                          spacetime.id = ctm_pm25$spacetime_id)

ctm_fit_cv
ctm_fit_cv$estimate[1000:1010]
ctm_pm25$pm25[1000:1010]
ctm_fit_cv$sd[1000:1010]
ctm_fit_cv$upper_95[1000:1010]

mean((ctm_fit_cv$estimate - ctm_pm25$pm25) ^ 2, na.rm = T)
mean((ctm_fit_cv_nngp$estimate - ctm_pm25$pm25) ^ 2, na.rm = T)


ctm_fit_cv_spat <- grm_cv(Y = ctm_pm25$pm25,
                          X = ctm_pm25$ctm,
                          cv.object = cv_id_ctm_spcl,
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

cv_id_maia_ord <- create_cv(space.id = maia_pm25$space_id,
                            time.id = maia_pm25$time_id,
                            type = "ordinary")

maia_fit <- grm(Y = maia_pm25$pm25,
                X = maia_pm25$aod,
                L = maia_pm25[, c("elevation", "forestcover",
                                  "hwy_length", "lim_hwy_length", 
                                  "local_rd_length", "point_emi_any")],
                M = maia_pm25[, c("tmp", "wind", "cmaq", "tempaod", 
                                  "windaod", "elevationaod")],
                n.iter = 500,
                nngp = T,
                num_neighbors = 6,
                burn = 100,
                thin = 4,
                coords = maia_pm25[, c("x", "y")],
                space.id = maia_pm25$space_id,
                time.id = maia_pm25$time_id,
                spacetime.id = maia_pm25$spacetime_id)

maia_fit_cv <- grm_cv(Y = maia_pm25$pm25,
                      X = maia_pm25$aod,
                      cv.object = cv_id_maia_ord,
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
