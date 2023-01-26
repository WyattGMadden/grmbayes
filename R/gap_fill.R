#' Make predictions for locations with two observations
#'
#' This functions uses the weights and predictions from the previous separate grm's to make ensemble predictions
#'
#' @inheritParams grm
#' @param Y.pred.1 Predictions made with grm_pred and first primary independent variable
#' @param Y.pred.2 Predictions made with grm_pred and second primary independent variable
#' @param Y.sd.1 Prediction standard deviation output from grm_pred and first primary independent variable
#' @param Y.sd.2 Prediction standard deviation output from grm_pred and second primary independent variable
#' @param weights Ensemble weights output from weight_pred 
#'
#' @return A list including ensemble predictions and standard deviations
#'
#' @examples
#' # gap_fill()
#' 
#' 
#' @export

gap_fill = function(Y.pred.1, 
                    Y.pred.2, 
                    Y.sd.1, 
                    Y.sd.2, 
                    space.id, 
                    weights) {
  
  W = 1 / (exp(-weights) + 1)
  W.mean =  apply(W, 1, mean)
  W.sd =  apply(1 / (exp(-weights) + 1), 1, stats::sd)
  
  which.use = which(!is.na(Y.pred.2))
  
  W.space = W.mean[space.id[which.use]] 
  
  Est = (Y.pred.1[which.use]) * W.space + Y.pred.2[which.use] * (1 - W.space)
  Est.SD = sqrt((Y.sd.1[which.use])^2 * W.space + 
                Y.sd.2[which.use]^2 * (1 - W.space) + 
                (Y.pred.1[which.use]^2 * W.space + 
                 Y.pred.2[which.use]^2 * (1 - W.space) - Est^2))
  
  
  Est.ensemb = Y.pred.1
  Est.ensemb[which.use] = Est
  SD.ensemb = Y.sd.1
  SD.ensemb[which.use] = Est.SD
  
  list(ensemble.estimate = Est.ensemb,
       ensemble.sd = SD.ensemb)
}
