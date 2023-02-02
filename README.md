
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

Y.input.ctm <- read.csv("../onedrive_code/Stage 1 Input Data/CTM/Y.csv")
X.input.ctm <- read.csv("../onedrive_code/Stage 1 Input Data/CTM/X.csv")
L.input.ctm <- as.matrix(read.csv("../onedrive_code/Stage 1 Input Data/CTM/L.csv"))
M.input.ctm <- as.matrix(read.csv ("../onedrive_code/Stage 1 Input Data/CTM/M.csv"))
monitor.locs.ctm <- as.matrix(read.csv("../onedrive_code/Stage 1 Input Data/CTM/Monitor_XY.csv"))
dist.mat.ctm <- as.matrix(dist(monitor.locs.ctm[, -1], 
                          diag = TRUE, 
                          upper = TRUE))

ctm.fit <- grm(Y = Y.input.ctm$pm25,
               X = X.input.ctm$CTM,
               L = L.input.ctm[, 4:ncol(L.input.ctm)],
               M = M.input.ctm[, 4:ncol(M.input.ctm)],
               n.iter = 500,
               burn = 100,
               thin = 4,
               dist.space.mat = dist.mat.ctm, 
               space.id = Y.input.ctm$Space_ID,
               time.id = Y.input.ctm$Time_ID,
               spacetime.id = Y.input.ctm$SpaceTime_ID)

ctm.fit.cv <- grm_cv(Y = Y.input.ctm$pm25,
                     X = X.input.ctm$CTM,
                     L = L.input.ctm[, 4:ncol(L.input.ctm)],
                     M = M.input.ctm[, 4:ncol(M.input.ctm)],
                     n.iter = 500,
                     burn = 100,
                     thin = 4,
                     dist.space.mat = dist.mat.ctm, 
                     space.id = Y.input.ctm$Space_ID,
                     time.id = Y.input.ctm$Time_ID,
                     spacetime.id = Y.input.ctm$SpaceTime_ID)

Y.input.maia <- read.csv("../onedrive_code/Stage 1 Input Data/MAIA/Y.csv")
X.input.maia <- read.csv("../onedrive_code/Stage 1 Input Data/MAIA/X.csv")
L.input.maia <- as.matrix(read.csv("../onedrive_code/Stage 1 Input Data/MAIA/L.csv"))
M.input.maia <- as.matrix(read.csv ("../onedrive_code/Stage 1 Input Data/MAIA/M.csv"))
monitor.locs.maia <- as.matrix(read.csv("../onedrive_code/Stage 1 Input Data/MAIA/Monitor_XY.csv"))
dist.mat.maia <- as.matrix(dist(monitor.locs.maia[, -1], 
                          diag = TRUE, 
                          upper = TRUE))

maia.fit <- grm(Y = Y.input.maia$pm25,
                X = X.input.maia$aod,
                L = L.input.maia[, 4:ncol(L.input.maia)],
                M = M.input.maia[, 4:ncol(M.input.maia)],
                n.iter = 500,
                burn = 100,
                thin = 4,
                dist.space.mat = dist.mat.maia, 
                space.id = Y.input.maia$Space_ID,
                time.id = Y.input.maia$Time_ID,
                spacetime.id = Y.input.maia$SpaceTime_ID)

maia.fit.cv <- grm_cv(Y = Y.input.maia$pm25,
                      X = X.input.maia$aod,
                      L = L.input.maia[, 4:ncol(L.input.maia)],
                      M = M.input.maia[, 4:ncol(M.input.maia)],
                      n.iter = 500,
                      burn = 100,
                      thin = 4,
                      dist.space.mat = dist.mat.maia, 
                      space.id = Y.input.maia$Space_ID,
                      time.id = Y.input.maia$Time_ID,
                      spacetime.id = Y.input.maia$SpaceTime_ID)
```

### Stage 2

``` r
pred.locs.ctm <- as.matrix(read.csv("../onedrive_code/Stage 2 Input data/CTM/Cell_XY.csv"))
L.pred.ctm <- as.matrix(read.csv("../onedrive_code/Stage 2 Input data/CTM/L.csv"))
M.pred.ctm <- as.matrix(read.csv("../onedrive_code/Stage 2 Input data/CTM/M.csv"))
X.pred.ctm <- read.csv("../onedrive_code/Stage 2 Input data/CTM/X.csv")

#standardize
for (i in 4:ncol(L.pred.ctm)) {
         L.pred.ctm[, i] <- (L.pred.ctm[, i] - mean(L.pred.ctm[, i])) / 
             sd(L.pred.ctm[, i])
}

for (i in 4:ncol(M.pred.ctm)) {
         M.pred.ctm[, i] <- (M.pred.ctm[, i] - mean(M.pred.ctm[, i])) / 
             sd(M.pred.ctm[, i])
}

X.pred.ctm$CTM <- (X.pred.ctm$CTM - mean(X.pred.ctm$CTM)) / 
    sd(X.pred.ctm$CTM)





ctm.pred = grm_pred(grm.fit = ctm.fit,
                    X.pred = X.pred.ctm$CTM, 
                    L.pred = L.pred.ctm[, 4:ncol(L.pred.ctm)], 
                    M.pred = M.pred.ctm[, 4:ncol(M.pred.ctm)], 
                    locations.Y = monitor.locs.ctm[, -1], 
                    locations.pred = pred.locs.ctm[, -1], 
                    space.id = X.pred.ctm$Space_ID, 
                    time.id = X.pred.ctm$Time_ID, 
                    spacetime.id = X.pred.ctm$SpaceTime_ID,
                    include.additive.annual.resid = T,
                    include.multiplicative.annual.resid = T,
                    n.iter = 100,
                    verbose = T)

pred.locs.maia = as.matrix(read.csv("../onedrive_code/Stage 2 Input data/MAIA/Cell_XY.csv"))
L.pred.maia = as.matrix(read.csv("../onedrive_code/Stage 2 Input data/MAIA/L.csv"))
M.pred.maia = as.matrix(read.csv("../onedrive_code/Stage 2 Input data/MAIA/M.csv"))
X.pred.maia = read.csv("../onedrive_code/Stage 2 Input data/MAIA/X.csv")

maia.pred = grm_pred(grm.fit = maia.fit,
                     X.pred = X.pred.maia$aod, 
                     L.pred = L.pred.maia[, 4:ncol(L.pred.maia)], 
                     M.pred = M.pred.maia[, 4:ncol(M.pred.maia)], 
                     locations.Y = monitor.locs.maia[, -1], 
                     locations.pred = pred.locs.maia[, -1], 
                     space.id = X.pred.maia$Space_ID, 
                     time.id = X.pred.maia$Time_ID, 
                     spacetime.id = X.pred.maia$SpaceTime_ID,
                     include.additive.annual.resid = T,
                     include.multiplicative.annual.resid = T,
                     n.iter = 100,
                     verbose = T)
#standardize
for (i in 4:ncol(L.pred.maia)) {
         L.pred.maia[, i] <- (L.pred.maia[, i] - mean(L.pred.maia[, i])) / 
             sd(L.pred.maia[, i])
}

for (i in 4:ncol(M.pred.maia)) {
         M.pred.maia[, i] <- (M.pred.maia[, i] - mean(M.pred.maia[, i])) / 
             sd(M.pred.maia[, i])
}

X.pred.maia$aod <- (X.pred.maia$aod - mean(X.pred.maia$aod)) / 
    sd(X.pred.maia$aod)
```

### Stage 3

``` r


# Read in out-of-sample predictions from CV analyses
ctm.input <- ctm.fit.cv
ctm.dateinfo <- read.csv("../onedrive_code/Stage 3 Input Data/CTM_Date_Mon_ID.csv")

# Check and amend date info
if (all (ctm.dateinfo$Time_ID == ctm.input$Time_ID & 
         ctm.dateinfo$Space_ID == ctm.input$Space_ID)){
  ctm.input$date <- ctm.dateinfo$Date
}

maia.input <- maia.fit.cv
maia.dateinfo <- read.csv("../onedrive_code/Stage 3 Input Data/MAIA_Date_Mon_ID.csv")
#Check and amend date info
if (all (maia.dateinfo$Time_ID == maia.input$time_id & 
         maia.dateinfo$Space_ID == maia.input$space_id)){
  maia.input$date = maia.dateinfo$Date
}


# Remove NA's from the first and last time interval
ctm.input <- subset(ctm.input, !is.na(estimate))
maia.input <- subset(maia.input, !is.na(estimate))

# Use only CTM results with AOD is observed
ctm.input$link_id <- paste(ctm.input$date, 
                          ctm.input$space_id, 
                          sep = "_")
maia.input$link_id <- paste(maia.input$date, 
                           maia.input$space_id, 
                           sep = "_")
ctm.input <- subset(ctm.input, 
                   link_id %in% maia.input$link_id)

maia.input$estimate_ctm <- ctm.input$estimate[match(ctm.input$link_id, 
                                              maia.input$link_id)]
maia.input$sd_ctm <- ctm.input$estimate[match(ctm.input$link_id, 
                                             maia.input$link_id)]

maia.input <- maia.input[order(maia.input$space_id, 
                              maia.input$date), ]


# Calculate densities
d1 <- dnorm(maia.input$obs, maia.input$estimate_ctm, maia.input$sd_ctm)
d2 <- dnorm(maia.input$obs, maia.input$estimate, maia.input$sd)


ensemble_fit <- ensemble_spatial(d1 = d1, 
                                 d2 = d2, 
                                 dist.space.mat = dist.mat.ctm, 
                                 space.id = maia.input$space_id, 
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
                            locations.Y = monitor.locs.ctm, 
                            locations.pred = pred.locs.ctm)

#Stage 2 predictions and Date info
ctm.dateinfo.pred <- read.csv ("../onedrive_code/Stage 4 Input Data/CTM_Pred_Date_Cell_ID.csv")
if (all (ctm.dateinfo.pred$Time_ID == ctm.pred$time.id & 
         ctm.dateinfo.pred$space.id == ctm.pred$space.id)) {
    ctm.pred$date <- ctm.dateinfo.pred$Date
}

maia.dateinfo.pred <- read.csv("../onedrive_code/Stage 4 Input Data/MAIA_Pred_Date_Cell_ID.csv")

if (all (maia.dateinfo.pred$time.id == maia.pred$time.id & 
         maia.dateinfo.pred$space.id == maia.pred$space.id)) {
    maia.pred$date <- maia.dateinfo.pred$Date
}

#Merge MAIA predictions onto CTM prediction dataset
ctm.pred$link_id <- paste(ctm.pred$date, 
                          ctm.pred$space.id, 
                          sep = "_")
maia.pred$link_id <- paste(maia.pred$date, 
                           maia.pred$space.id, 
                           sep = "_")

ctm.pred$estimate_maia <- maia.pred$estimate[match(ctm.pred$estimate, 
                                                   maia.pred$estimate)]
ctm.pred$sd_maia <- maia.pred$sd[match(ctm.pred$link_id, 
                                       maia.pred$link_id)]

# Perform Gap Filling
results <- gap_fill(Y.pred.1 = ctm.pred$estimate, 
                    Y.pred.2 = ctm.pred$estimate_maia, 
                    Y.sd.1 = ctm.pred$sd, 
                    Y.sd.2 = ctm.pred$sd_maia, 
                    space.id = ctm.pred$space.id, 
                    weights = weight_preds)
```
