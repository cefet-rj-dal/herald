#Case A:
# - model: MA
# - prediction: MA of last #order known points (input: serie)
# - detection: model+prediction

#Case C:
# - model: MA
# - prediction: MA of last #order known points (input: serie + last prediction)
# - detection: model+prediction+prediction

#Case B:
# - model: lag
# - prediction: last known point (input: serie)
# - detection: model+prediction

#Case D:
# - model: lag
# - prediction: last known point (input: serie + last prediction)
# - detection: model+prediction+prediction


#General Case:

# - prediction: pred over known points
pred <- function(x) tail(x,1) # Direct prediction
pred_ma <- function(x,order) mean(tail(x,order)) # MA prediction

# - model: pred called in loop over the serie (finish to end)
model <- function(serie, pred_func, ...){
  mdl <- list()
  for(i in length(serie):1){
    mdl[i] <- do.call(pred_func, args=c(list(serie[1:(i-1)]), ...))
  }
  return(unlist(mdl))
}
# - lagged forecast: pred called in loop (for each lag) over last known serie+last prediction points
lag_pred <- function(serie, pred_func, lag, ...){
  
  if(length(serie)-lag-1<=0) lag <- 0
  serie <- serie[1:(length(serie)-lag-1)]
  serie_mdl <- model(serie,pred_func,...)
  
  serie_w_pred <- serie
  for(i in 1:lag) serie_w_pred <- c(serie_w_pred, do.call(pred_func, args=c(list(serie_w_pred), ...)))
  
  new_mdl <- c(serie_mdl, tail(serie_w_pred,lag), do.call(pred_func, args=c(list(serie_w_pred), ...)) )
  
  return(new_mdl)
}
# - lagged forecast recursive: lag_pred called in loop over the serie (finish to end)
model_lag_pred <- function(serie, pred_func, lag, ...){
  mdl <- list()
  for(i in length(serie):1){
    mdl[i] <- tail( do.call(lag_pred, args=c(list(serie[1:i]), pred_func=pred_func, lag=lag, ...)) ,1)
  }
  return(unlist(mdl))
}

serie <- c(0.7,2.8,2.7,6.2,8.5)

serie
model(serie,pred)
lag_pred(serie, pred_ma, lag=0,order=3)
model_lag_pred(serie, pred_ma, lag=2,order=3)

plot(1:length(serie), serie, pch=19, col="blue")
points(1:length(serie), model(serie,pred), pch=19, col="red")
points(1:length(serie), model(serie,pred_ma, order=3), pch=19, col="green")

plot(1:length(serie),serie-model(serie,pred_ma, order=3), type = "l", col="green")
lines(1:length(serie),serie-model(serie,pred), type = "l", col="red")



herald_mdl <- function(serie, lag=1, ic=FALSE, recursive=FALSE, ...){
  
  pred_ma <- function(x,order) mean(tail(x,order)) # MA prediction
  pred <- function(x) tail(x,1) # Direct prediction
  
  if(ic) pred_func <- pred_ma
  else pred_func <- pred
  
  model <- function(serie, pred_func, ...){
    mdl <- list()
    for(i in length(serie):1){
      mdl[i] <- do.call(pred_func, args=c(list(serie[1:(i-1)]), ...))
    }
    return(unlist(mdl))
  }
  
  lag_pred <- function(serie, pred_func, lag, ...){
    
    if(length(serie)-lag-1<=0) lag <- 0
    serie <- serie[1:(length(serie)-lag-1)]
    serie_mdl <- model(serie,pred_func,...)
    
    serie_w_pred <- serie
    for(i in 1:lag) serie_w_pred <- c(serie_w_pred, do.call(pred_func, args=c(list(serie_w_pred), ...)))
    
    new_mdl <- c(serie_mdl, tail(serie_w_pred,lag), do.call(pred_func, args=c(list(serie_w_pred), ...)) )
    
    return(new_mdl)
  }
  
  model_lag_pred <- function(serie, pred_func, lag, ...){
    mdl <- list()
    for(i in length(serie):1){
      mdl[i] <- tail( do.call(lag_pred, args=c(list(serie[1:i]), pred_func=pred_func, lag=lag, ...)) ,1)
    }
    return(unlist(mdl))
  }
  
  if(recursive) return(model_lag_pred(serie, pred_func, lag=lag,...))
  else return(lag_pred(serie, pred_func, lag=lag,...))
  
}

# Detection over predicted inertial component -----------------------------------------
detector <- function(serie,pred,alpha = 1.5){
  n <- nrow(serie)
  non_na <- which(!is.na(serie))
  
  serie <- na.omit(serie)
  
  #===== Boxplot analysis of results ======
  outliers.boxplot <- function(data, alpha = 1.5){
    org = length(data)
    cond <- rep(FALSE, org)
    q = quantile(data, na.rm=TRUE)
    IQR = q[4] - q[2]
    lq1 = as.double(q[2] - alpha*IQR)
    hq3 = as.double(q[4] + alpha*IQR)
    #cond = data < lq1 | data > hq3
    cond = data > hq3
    return (cond)
  }
  
  outliers.boxplot.index <- function(data, alpha = 1.5){
    cond <- outliers.boxplot(data, alpha)
    index.cp = which(cond)
    return (index.cp)
  }
  
  #Adjustment error on the entire series
  s <- (tail(serie,nrow(pred))-pred)^2
  outliers <- outliers.boxplot.index(s, alpha)
  group_outliers <- split(outliers, cumsum(c(1, diff(outliers) != 1)))
  outliers <- rep(FALSE, nrow(s))
  for (g in group_outliers) {
    if (length(g) > 0) {
      i <- min(g)
      outliers[i] <- TRUE
    }
  }
  outliers <- c(rep(FALSE, n-nrow(s)),outliers)
  i_outliers <- rep(NA, n)
  i_outliers[non_na] <- outliers
  
  detection <- data.frame(idx=1:n, event = i_outliers, type="")
  detection$type[i_outliers] <- "change point"
  
  return(detection)
}



serie <- c(0.7,2.8,2.7,6.2,8.5)

pred <- herald_tests(serie, lag=2, ic=TRUE, recursive=FALSE,order=3) #Working

detector(serie,pred,alpha = 1.5)
