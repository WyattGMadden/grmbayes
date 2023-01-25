
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ensembleDownscaleR

<!-- badges: start -->

[![R-CMD-check](https://github.com/WyattGMadden/ensembleDownscaleR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/WyattGMadden/ensembleDownscaleR/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The goal of ensembleDownscaleR is to …

## Installation

You can install the development version of ensembleDownscaleR like so:

``` r
install_github("WyattGMadden/ensembleDownscaleR")
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
pred.locs.ctm = as.matrix(read.csv("../onedrive_code/Stage 2 Input data/CTM/Cell_XY.csv"))
L.pred.ctm = as.matrix(read.csv("../onedrive_code/Stage 2 Input data/CTM/L.csv"))
M.pred.ctm = as.matrix(read.csv("../onedrive_code/Stage 2 Input data/CTM/M.csv"))
X.pred.ctm = read.csv("../onedrive_code/Stage 2 Input data/CTM/X.csv")

test_pred = grm_pred(grm.fit = ctm.fit,
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

test_pred = grm_pred(grm.fit = maia.fit,
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
```

### Stage 3

``` r


##Read in out-of-sample predictions from CV analyses
CTM.input = ctm.fit.cv
CTM.DateInfo = read.csv("../onedrive_code/Stage 3 Input Data/CTM_Date_Mon_ID.csv")

#Check and amend date info
if (all (CTM.DateInfo$Time_ID == CTM.input$Time_ID & CTM.DateInfo$Space_ID == CTM.input$Space_ID)){
  CTM.input$date = CTM.DateInfo$Date
}

MAIA.input = maia.fit.cv
MAIA.DateInfo = read.csv("../onedrive_code/Stage 3 Input Data/MAIA_Date_Mon_ID.csv")
#Check and amend date info
if (all (MAIA.DateInfo$Time_ID == MAIA.input$time_id & MAIA.DateInfo$Space_ID == MAIA.input$space_id)){
  MAIA.input$date = MAIA.DateInfo$Date
}


#Remove NA's from the first and last time interval
CTM.input = subset(CTM.input, !is.na(estimate))
MAIA.input = subset(MAIA.input, !is.na(estimate))

#Use only CTM results with AOD is observed
CTM.input$link_id = paste(CTM.input$date, 
                          CTM.input$space_id, 
                          sep = "_")
MAIA.input$link_id = paste(MAIA.input$date, 
                           MAIA.input$space_id, 
                           sep = "_")
CTM.input = subset(CTM.input, 
                   link_id %in% MAIA.input$link_id)

MAIA.input$estimate_ctm = CTM.input$estimate[match(CTM.input$link_id, 
                                              MAIA.input$link_id)]
MAIA.input$sd_ctm = CTM.input$estimate[match(CTM.input$link_id, 
                                             MAIA.input$link_id)]

MAIA.input = MAIA.input[order(MAIA.input$space_id, 
                              MAIA.input$date), ]


#Calculate densities
d1 = dnorm(MAIA.input$obs, MAIA.input$estimate_ctm, MAIA.input$sd_ctm)
d2 = dnorm(MAIA.input$obs, MAIA.input$estimate, MAIA.input$sd)


ensemble_fit = ensemble_spatial(d1 = d1, 
                                d2 = d2, 
                                dist.space.mat = dist.mat.ctm, 
                                space.id = MAIA.input$space_id, 
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

#Monitoring locations and all grid cell locations
Monitor_XY <- monitor.locs.ctm
Cell_XY <- pred.locs.ctm


q.samples <- ensemble_fit$q


theta.samples <- ensemble_fit$other$theta
tau2.samples <- ensemble_fit$other$tau2

weight_preds <- weight_pred(q = ensemble_fit$q, 
                            theta = ensemble_fit$other$theta,
                            tau2 = ensemble_fit$other$tau2,
                            locations.Y = monitor.locs.ctm, 
                            locations.pred = pred.locs.ctm)
```
