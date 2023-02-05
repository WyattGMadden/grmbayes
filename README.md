
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

str(ctm_pm25)

var_names_ctm_L <- c("elevation", "forestcover", "hwy_length",
                 "lim_hwy_length", "local_rd_length", "point_emi_any")
var_names_ctm_M <- c("tmp", "wind")

dist_dat_ctm_pm25 <- unique(ctm_pm25[, c("space_id", "x", "y")])
dist_dat_ctm_pm25 <- dist_dat_ctm_pm25[order(dist_dat_ctm_pm25$space_id), ]
dist_mat_ctm_pm25 <- as.matrix(dist(dist_dat_ctm_pm25[, c("x", "y")], 
                                    diag = TRUE, 
                                    upper = TRUE))

ctm_fit <- grm(Y = ctm_pm25$pm25,
               X = ctm_pm25$ctm,
               L = as.matrix(ctm_pm25[, var_names_ctm_L]),
               M = as.matrix(ctm_pm25[, var_names_ctm_M]),
               n.iter = 500,
               burn = 100,
               thin = 4,
               dist.space.mat = dist_mat_ctm_pm25,
               space.id = ctm_pm25$space_id,
               time.id = ctm_pm25$time_id,
               spacetime.id = ctm_pm25$spacetime_id)

ctm_fit_cv <- grm_cv(Y = ctm_pm25$pm25,
                     X = ctm_pm25$pm25,
                     L = as.matrix(ctm_pm25[, var_names_ctm_L]),
                     M = as.matrix(ctm_pm25[, var_names_ctm_M]),
                     n.iter = 500,
                     burn = 100,
                     thin = 4,
                     dist.space.mat = dist_mat_ctm_pm25, 
                     space.id = ctm_pm25$space_id,
                     time.id = ctm_pm25$time_id,
                     spacetime.id = ctm_pm25$spacetime_id)


str(maia_pm25)

var_names_maia_L <- c("elevation", "forestcover", "hwy_length",
                 "lim_hwy_length", "local_rd_length", "point_emi_any")
var_names_maia_M <- c("tmp", "wind", "cmaq", 
                      "tempaod", "windaod", "elevationaod")

dist_dat_maia_pm25 <- unique(maia_pm25[, c("space_id", "x", "y")])
dist_dat_maia_pm25 <- dist_dat_maia_pm25[order(dist_dat_maia_pm25$space_id), ]
dist_mat_maia_pm25 <- as.matrix(dist(dist_dat_maia_pm25[, c("x", "y")], 
                                    diag = TRUE, 
                                    upper = TRUE))

maia_fit <- grm(Y = maia_pm25$pm25,
                X = maia_pm25$aod,
                L = as.matrix(maia_pm25[, var_names_maia_L]),
                M = as.matrix(maia_pm25[, var_names_maia_M]),
                n.iter = 500,
                burn = 100,
                thin = 4,
                dist.space.mat = dist_mat_maia_pm25, 
                space.id = maia_pm25$space_id,
                time.id = maia_pm25$time_id,
                spacetime.id = maia_pm25$spacetime_id)

maia_fit_cv <- grm_cv(Y = maia_pm25$pm25,
                      X = maia_pm25$aod,
                      L = as.matrix(maia_pm25[, var_names_maia_L]),
                      M = as.matrix(maia_pm25[, var_names_maia_M]),
                      n.iter = 500,
                      burn = 100,
                      thin = 4,
                      dist.space.mat = dist_mat_maia_pm25, 
                      space.id = maia_pm25$space_id,
                      time.id = maia_pm25$time_id,
                      spacetime.id = maia_pm25$spacetime_id)
```

### Stage 2

``` r
str(ctm_preds)

dist_dat_ctm_preds <- unique(ctm_preds[, c("space_id", "x", "y")])
dist_dat_ctm_preds <- dist_dat_ctm_preds[order(dist_dat_ctm_preds$space_id), ]
dist_mat_ctm_preds <- as.matrix(dist(dist_dat_ctm_preds[, c("x", "y")], 
                                    diag = TRUE, 
                                    upper = TRUE))

#standardize
for (i in var_names_ctm_L) {
         ctm_preds[, i] <- (ctm_preds[, i] - mean(ctm_preds[, i])) / 
             sd(ctm_preds[, i])
}

for (i in var_names_ctm_M) {
         ctm_preds[, i] <- (ctm_preds[, i] - mean(ctm_preds[, i])) / 
             sd(ctm_preds[, i])
}

ctm_preds$ctm <- (ctm_preds$ctm - mean(ctm_preds$ctm)) / 
    sd(ctm_preds$ctm)





ctm_pred <- grm_pred(grm.fit = ctm_fit,
                     X.pred = ctm_preds$ctm, 
                     L.pred = as.matrix(ctm_preds[, var_names_ctm_L]),
                     M.pred = as.matrix(ctm_preds[, var_names_ctm_M]),
                     locations.Y = dist_dat_ctm_pm25[, c("x", "y")],
                     locations.pred = dist_dat_ctm_preds[, c("x", "y")],
                     space.id = ctm_preds$space_id,
                     time.id = ctm_preds$time_id,
                     spacetime.id = ctm_preds$spacetime_id,
                     include.additive.annual.resid = T,
                     include.multiplicative.annual.resid = T,
                     n.iter = 100,
                     verbose = T)



str(maia_preds)

dist_dat_maia_preds <- unique(maia_preds[, c("space_id", "x", "y")])
dist_dat_maia_preds <- dist_dat_maia_preds[order(dist_dat_maia_preds$space_id), ]
dist_mat_maia_preds <- as.matrix(dist(dist_dat_maia_preds[, c("x", "y")], 
                                      diag = TRUE, 
                                      upper = TRUE))

#standardize
for (i in var_names_maia_L) {
         maia_preds[, i] <- (maia_preds[, i] - mean(maia_preds[, i])) / 
             sd(maia_preds[, i])
}

for (i in var_names_maia_M) {
         maia_preds[, i] <- (maia_preds[, i] - mean(maia_preds[, i])) / 
             sd(maia_preds[, i])
}

maia_preds$aod <- (maia_preds$aod - mean(maia_preds$aod)) / 
    sd(maia_preds$aod)

maia_pred <- grm_pred(grm.fit = maia_fit,
                      X.pred = maia_preds$aod, 
                      L.pred = as.matrix(maia_preds[, var_names_maia_L]),
                      M.pred = as.matrix(maia_preds[, var_names_maia_M]),
                      locations.Y = dist_dat_maia_pm25[, c("x", "y")],
                      locations.pred = dist_dat_maia_preds[, c("x", "y")],
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

ctm_input <- ctm_fit_cv
maia_input <- maia_fit_cv

ctm_input$date <- ctm_pm25$date
maia_input$date <- maia_pm25$date

# Remove NA's from the first and last time interval
ctm_input <- subset(ctm_input, !is.na(estimate))
maia_input <- subset(maia_input, !is.na(estimate))

# Use only CTM results with AOD is observed
ctm_input$link_id <- paste(ctm_input$date, 
                          ctm_input$space_id, 
                          sep = "_")
maia_input$link_id <- paste(maia_input$date, 
                           maia_input$space_id, 
                           sep = "_")
ctm_input <- subset(ctm_input, 
                   link_id %in% maia_input$link_id)

maia_input$estimate_ctm <- ctm_input$estimate[match(ctm_input$link_id, 
                                                    maia_input$link_id)]
maia_input$sd_ctm <- ctm_input$estimate[match(ctm_input$link_id, 
                                              maia_input$link_id)]

maia_input <- maia_input[order(maia_input$space_id, 
                              maia_input$date), ]


# Calculate densities
d1 <- dnorm(maia_input$obs, maia_input$estimate_ctm, maia_input$sd_ctm)
d2 <- dnorm(maia_input$obs, maia_input$estimate, maia_input$sd)


ensemble_fit <- ensemble_spatial(d1 = d1, 
                                 d2 = d2, 
                                 dist.space.mat = dist_mat_ctm_pm25, 
                                 space.id = maia_input$space_id, 
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

weight_preds <- weight_pred(q = ensemble_fit$q, 
                            theta = ensemble_fit$other$theta,
                            tau2 = ensemble_fit$other$tau2,
                            locations.Y = dist_dat_ctm_pm25, 
                            locations.pred = dist_dat_ctm_preds)

ctm_pred$date <- ctm_preds$date
maia_pred$date <- maia_preds$date

#Merge MAIA predictions onto CTM prediction dataset
ctm_pred$link_id <- paste(ctm_pred$date, 
                          ctm_pred$space_id, 
                          sep = "_")
maia_pred$link_id <- paste(maia_pred$date, 
                           maia_pred$space_id, 
                           sep = "_")

ctm_pred$estimate_maia <- maia_pred$estimate[match(ctm_pred$estimate, 
                                                   maia_pred$estimate)]
ctm_pred$sd_maia <- maia_pred$sd[match(ctm_pred$link_id, 
                                       maia_pred$link_id)]

# Perform Gap Filling
results <- gap_fill(Y.pred.1 = ctm_pred$estimate, 
                    Y.pred.2 = ctm_pred$estimate_maia, 
                    Y.sd.1 = ctm_pred$sd, 
                    Y.sd.2 = ctm_pred$sd_maia, 
                    space.id = ctm_pred$space_id, 
                    weights = weight_preds)
```
