
###STAGE 1###

#ctm 
Y.input.ctm <- read.csv("../onedrive_code/Stage 1 Input Data/CTM/Y.csv")
X.input.ctm <- read.csv("../onedrive_code/Stage 1 Input Data/CTM/X.csv")
L.input.ctm <- as.matrix(read.csv("../onedrive_code/Stage 1 Input Data/CTM/L.csv"))
M.input.ctm <- as.matrix(read.csv ("../onedrive_code/Stage 1 Input Data/CTM/M.csv"))
monitor.locs.ctm <- as.matrix(read.csv("../onedrive_code/Stage 1 Input Data/CTM/Monitor_XY.csv"))
ctm.dateinfo <- read.csv("../onedrive_code/Stage 3 Input Data/CTM_Date_Mon_ID.csv")

ctm_pm25 <- Y.input.ctm
ctm_pm25$ctm <- X.input.ctm$CTM

L_ctm_names <- colnames(L.input.ctm)[4:ncol(L.input.ctm)]
ctm_pm25[, L_ctm_names] <- L.input.ctm[, L_ctm_names]

M_ctm_names <- colnames(M.input.ctm)[4:ncol(M.input.ctm)]
ctm_pm25[, M_ctm_names] <- M.input.ctm[, M_ctm_names]

ctm_pm25$Date <- ctm.dateinfo$Date

monitor.locs.ctm <- as.data.frame(monitor.locs.ctm)

ctm_pm25 <- merge(ctm_pm25, 
                 monitor.locs.ctm, 
                 by.x = "Space_ID", 
                 by.y = "ID",
                 all.x = T,
                 sort = F)

names(ctm_pm25) <- tolower(names(ctm_pm25))
names(ctm_pm25) <- gsub("\\.", "_", names(ctm_pm25))

ctm_pm25 <- ctm_pm25[, c(c("time_id", "space_id", "spacetime_id"), 
                         names(ctm_pm25)[4:ncol(ctm_pm25)])]
ctm_pm25 <- ctm_pm25[with(ctm_pm25, order(time_id, space_id, spacetime_id)), ]



usethis::use_data(ctm_pm25, overwrite = TRUE)


#maia
Y.input.maia <- read.csv("../onedrive_code/Stage 1 Input Data/MAIA/Y.csv")
X.input.maia <- read.csv("../onedrive_code/Stage 1 Input Data/MAIA/X.csv")
L.input.maia <- as.matrix(read.csv("../onedrive_code/Stage 1 Input Data/MAIA/L.csv"))
M.input.maia <- as.matrix(read.csv ("../onedrive_code/Stage 1 Input Data/MAIA/M.csv"))
monitor.locs.maia <- as.matrix(read.csv("../onedrive_code/Stage 1 Input Data/MAIA/Monitor_XY.csv"))
maia.dateinfo <- read.csv("../onedrive_code/Stage 3 Input Data/MAIA_Date_Mon_ID.csv")

maia_pm25 <- Y.input.maia
maia_pm25$aod <- X.input.maia$aod

L_maia_names <- colnames(L.input.maia)[4:ncol(L.input.maia)]
maia_pm25[, L_maia_names] <- L.input.maia[, L_maia_names]

M_maia_names <- colnames(M.input.maia)[4:ncol(M.input.maia)]
maia_pm25[, M_maia_names] <- M.input.maia[, M_maia_names]

maia_pm25$Date <- maia.dateinfo$Date

monitor.locs.maia <- as.data.frame(monitor.locs.maia)

maia_pm25 <- merge(maia_pm25, 
                 monitor.locs.ctm, 
                 by.x = "Space_ID", 
                 by.y = "ID",
                 all.x = T,
                 sort = F)


names(maia_pm25) <- tolower(names(maia_pm25))
names(maia_pm25) <- gsub("\\.", "_", names(maia_pm25))

maia_pm25 <- maia_pm25[, c(c("time_id", "space_id", "spacetime_id"), 
                         names(maia_pm25)[4:ncol(maia_pm25)])]
maia_pm25 <- maia_pm25[with(maia_pm25, order(time_id, space_id, spacetime_id)), ]

usethis::use_data(maia_pm25, overwrite = TRUE)


###STAGE 2###

#ctm
pred.locs.ctm <- as.matrix(read.csv("../onedrive_code/Stage 2 Input data/CTM/Cell_XY.csv"))
L.pred.ctm <- as.matrix(read.csv("../onedrive_code/Stage 2 Input data/CTM/L.csv"))
M.pred.ctm <- as.matrix(read.csv("../onedrive_code/Stage 2 Input data/CTM/M.csv"))
X.pred.ctm <- read.csv("../onedrive_code/Stage 2 Input data/CTM/X.csv")
ctm.dateinfo.pred <- read.csv ("../onedrive_code/Stage 4 Input Data/CTM_Pred_Date_Cell_ID.csv")

ctm_preds <- X.pred.ctm

ctm_preds[, M_ctm_names] <- M.pred.ctm[, M_ctm_names]
ctm_preds[, L_ctm_names] <- L.pred.ctm[, L_ctm_names]
ctm_preds$Date <- ctm.dateinfo.pred$Date

pred.locs.ctm <- as.data.frame(pred.locs.ctm)

ctm_preds <- merge(ctm_preds, 
                 pred.locs.ctm, 
                 by.x = "Space_ID", 
                 by.y = "Cell",
                 all.x = T,
                 sort = F)


names(ctm_preds) <- tolower(names(ctm_preds))
names(ctm_preds) <- gsub("\\.", "_", names(ctm_preds))

ctm_preds <- ctm_preds[, c(c("time_id", "space_id", "spacetime_id"), 
                         names(ctm_preds)[4:ncol(ctm_preds)])]
ctm_preds <- ctm_preds[with(ctm_preds, order(time_id, space_id, spacetime_id)), ]

usethis::use_data(ctm_preds, overwrite = TRUE)



#maia
pred.locs.maia = as.matrix(read.csv("../onedrive_code/Stage 2 Input data/MAIA/Cell_XY.csv"))
L.pred.maia = as.matrix(read.csv("../onedrive_code/Stage 2 Input data/MAIA/L.csv"))
M.pred.maia = as.matrix(read.csv("../onedrive_code/Stage 2 Input data/MAIA/M.csv"))
X.pred.maia = read.csv("../onedrive_code/Stage 2 Input data/MAIA/X.csv")
maia.dateinfo.pred <- read.csv("../onedrive_code/Stage 4 Input Data/MAIA_Pred_Date_Cell_ID.csv")

maia_preds <- X.pred.maia

maia_preds[, M_maia_names] <- M.pred.maia[, M_maia_names]
maia_preds[, L_maia_names] <- L.pred.maia[, L_maia_names]
maia_preds$Date <- maia.dateinfo.pred$Date

pred.locs.maia <- as.data.frame(pred.locs.maia)

maia_preds <- merge(maia_preds, 
                 pred.locs.maia, 
                 by.x = "Space_ID", 
                 by.y = "Cell",
                 all.x = T,
                 sort = F)


names(maia_preds) <- tolower(names(maia_preds))
names(maia_preds) <- gsub("\\.", "_", names(maia_preds))

maia_preds <- maia_preds[, c(c("time_id", "space_id", "spacetime_id"), 
                         names(maia_preds)[4:ncol(maia_preds)])]
maia_preds <- maia_preds[with(maia_preds, order(time_id, space_id, spacetime_id)), ]


usethis::use_data(maia_preds, overwrite = TRUE)



