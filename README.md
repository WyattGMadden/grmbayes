
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ensembleDownscaleR

<!-- badges: start -->

[![R-CMD-check](https://github.com/WyattGMadden/ensembleDownscaleR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/WyattGMadden/ensembleDownscaleR/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The goal of ensembleDownscaleR is to …

## Installation

You can install the development version of ensembleDownscaleR like so:

``` r
# repo is currently private
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
