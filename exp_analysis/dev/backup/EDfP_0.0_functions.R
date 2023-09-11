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
har_edfp <- function(online_step=1,
                     pred_fun, pred_par,
                     detect_fun, detect_par) {
  obj <- harbinger()
  obj$n.ahead <- online_step
  obj$pred_fun <- pred_fun
  obj$pred_par <- pred_par
  obj$detect_fun <- detect_fun
  obj$detect_par <- detect_par
  class(obj) <- append("har_edfp", class(obj))
  return(obj)
}


#'@export
detect.har_edfp <- function(obj, serie) {
  if(is.null(serie)) stop("No data was provided for computation",call. = FALSE)
  
  n <- length(serie)
  non_na <- which(!is.na(serie))
  
  serie <- na.omit(serie)
  serie <- as.data.frame(serie)
  
  train <- head(serie,nrow(serie)-1)#-obj$n.ahead
  
  #prediction of inertial component
  pred <- do.call(pred_fun, args=c(list(serie), n.ahead=obj$n.ahead, pred_par))
  #browser()
  serie_w_pred <- as.data.frame(c(unlist(train),tail(unlist(pred),1)))
  
  
  #detection over predicted inertial component
  detection <- do.call(detect_fun, args=c(list(serie),list(serie_w_pred), detect_par))
  
  return(detection)
}
