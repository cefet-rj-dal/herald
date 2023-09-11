# Installing and loading Harbinger -------------------------------------------------------

# install.packages("devtools")
#library(devtools)
#devtools::install_github("cefet-rj-dal/harbinger", force=TRUE, dependencies=FALSE, upgrade="never", build_vignettes = TRUE)
#Loading harbinger
require(harbinger)


#'@description Ancestor class for time series event detection
#'@details The Harbinger class establishes the basic interface for time series event detection.
#'  Each method should be implemented in a descendant class of Harbinger
#'@return Harbinger object
#'@examples detector <- harbinger()
#'@export
har_herald <- function(lag_pred = 1, online_step=1, 
                       decomp_fun, decomp_par,
                       pred_fun, pred_par,
                       detect_fun, detect_par) {
  obj <- harbinger()
  obj$n.ahead <- lag_pred + online_step
  obj$decomp_fun <- decomp_fun
  obj$decomp_par <- decomp_par
  obj$pred_fun <- pred_fun
  obj$pred_par <- pred_par
  obj$detect_fun <- detect_fun
  obj$detect_par <- detect_par
  class(obj) <- append("har_herald", class(obj))
  return(obj)
}


#'@export
detect.har_herald <- function(obj, serie) {
  if(is.null(serie)) stop("No data was provided for computation",call. = FALSE)
  
  n <- length(serie)
  non_na <- which(!is.na(serie))
  
  serie <- na.omit(serie)
  
  #inertial component
  ic <- do.call(decomp_fun, args=c(list(serie), decomp_par))
  ic <- as.data.frame(na.omit(ic))
  
  train_ic <- head(ic,nrow(ic)-obj$n.ahead)
  test_ic <- tail(ic,obj$n.ahead)
  
  #plot(serie,type='l')
  #lines(ic,type='l',col='blue')
  
  #prediction of inertial component
  pred_ic <- do.call(pred_fun, args=c(list(ic), n.ahead=obj$n.ahead, pred_par))
  
  ic_w_pred <- as.data.frame(c(unlist(train_ic),unlist(pred_ic)))
  
  serie <- as.data.frame(serie)
  #detection over predicted inertial component
  detection <- do.call(detect_fun, args=c(list(serie),list(ic_w_pred), detect_par))
  
  return(detection)
}
