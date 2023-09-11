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

# Installing and loading Harbinger -------------------------------------------------------

#install.packages("devtools")
#library(devtools)
#devtools::install_github("cefet-rj-dal/harbinger", force=TRUE, dependencies=TRUE, upgrade="never", build_vignettes = TRUE)
#Loading harbinger
#require(harbinger)
source("https://raw.githubusercontent.com/cefet-rj-dal/harbinger-examples/main/jupyter.R")
load_library("daltoolbox")
load_github("cefet-rj-dal/harbinger")


#'@description Ancestor class for time series event detection
#'@details The Harbinger class establishes the basic interface for time series event detection.
#'  Each method should be implemented in a descendant class of Harbinger
#'@return Harbinger object
#'@examples detector <- harbinger()
#'@export
har_herald <- function(lag_pred = 1, online_step=1, 
                       pred_ic=FALSE, pred_par=NULL,
                       pred_recursive=FALSE,
                       alpha=1.5) {
  obj <- harbinger()
  obj$lag_pred <- lag_pred
  obj$online_step <- online_step
  obj$n.ahead <- lag_pred + online_step
  obj$pred_ic <- pred_ic
  obj$pred_par <- pred_par
  obj$pred_recursive <- pred_recursive
  obj$alpha <- alpha
  class(obj) <- append("har_herald", class(obj))
  return(obj)
}


#'@export
detect.har_herald <- function(obj, serie) {
  if(is.null(serie)) stop("No data was provided for computation",call. = FALSE)
  
  n <- length(serie)
  non_na <- which(!is.na(serie))
  
  serie <- na.omit(serie)
  
  # Model and prediction
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
  
  # Detection
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
  
  mdl_serie <- herald_mdl(serie, lag=obj$lag_pred, ic=obj$pred_ic, recursive=obj$pred_recursive, obj$pred_par)
  
  serie <- as.data.frame(serie)
  mdl_serie <- as.data.frame(mdl_serie)
  
  #detection over predicted inertial component
  detection <- detector(serie,mdl_serie,alpha = obj$alpha)
  
  return(detection)
}


# Nexus function ---------------------------------------------------------------

# Loading Nexus
source("https://raw.githubusercontent.com/cefet-rj-dal-dev/nexus/main/R/nexus.R?token=GHSAT0AAAAAACEGUDHHIWADD4XLKTDV6L3GZFFZI7A")

run_nexus <- function(model,data, warm_size = 30, batch_size = 15, mem_batches = 0, png_folder="dev/plots/", verbose=TRUE){
  datasource <- nex_simulated_datasource("data", data$serie)
  
  online_detector <- nexus(datasource, model, warm_size = warm_size, batch_size = batch_size, mem_batches = mem_batches)
  
  online_detector <- warmup(online_detector)
  browser()
  bt <- 1
  event_happened <- FALSE
  bt_event_happened <- 0
  event_idx <- 0
  while (!is.null(online_detector$datasource)) {
    online_detector <- detect(online_detector)
    if(verbose) print(paste0("Batch ",bt))

    #for plotting
    if(any(data$event[online_detector$detection$idx]) & !event_happened){
      if(verbose) print("Event happened!")
      event_happened <- TRUE
      event_idx <- which(data$event[online_detector$detection$idx])
      bt_event_happened <- bt
      time_event_happened <- Sys.time()
    }
    #browser()
    herald_serie <- herald_mdl(online_detector$serie, lag=model$lag_pred, ic=model$pred_ic, recursive=model$pred_recursive, model$pred_par)
    png(file=paste0(png_folder,"batch_",bt,".png"),
        width=600, height=350)
    grf <- 
      plot.harbinger(online_detector$detector, data$serie[online_detector$detection$idx], online_detector$detection, data$event[online_detector$detection$idx])+
      geom_line(aes(x=online_detector$detection$idx,y=c(rep(NA,length(online_detector$detection$idx)-length(herald_serie)),herald_serie)),col="orange")+
      geom_point(aes(x=online_detector$detection$idx,y=c(rep(NA,length(online_detector$detection$idx)-model$n.ahead),tail(herald_serie,model$n.ahead))),col="orange")
    
    if(any(which(as.logical(online_detector$detection$event))>event_idx) & event_happened) {
      if(verbose) print("Detected event!")
      online_detector$detection_lag <- list()
      online_detector$detection_lag$batch <- bt-bt_event_happened
      online_detector$detection_lag$obs <- min(which(as.logical(online_detector$detection$event)))-event_idx
      online_detector$detection_lag$time <- round(as.numeric(difftime(time1 = Sys.time(), time2 = time_event_happened, units = "secs")), 3)
    }
    plot(grf)
    dev.off()
    bt <- bt + 1
    if(any(which(as.logical(online_detector$detection$event))>event_idx) & event_happened) {
      break
    }#comment this line for full result 
    #end
  }
  
  return(online_detector)
}




#Experiments ===================================================================

# Simulated time series --------------------------------------------------------
nonstationarity_sym <- function(ts.len,ts.mean=0,ts.var=1,ts.trend=0.04,seed=12345){#1234
  require(tseries)
  require(forecast)
  
  #x variable (time)
  t <- c(1:ts.len)
  
  #stationary time series
  set.seed(seed)
  sta <- arima.sim(model=list(ar=0.5), n=ts.len, mean=ts.mean, sd=sqrt(ts.var)) #AR(1)
  
  #test of stationarity OK!
  #adf.test(sta)
  
  #trend-stationary time series
  trend <- ts.trend*t
  tsta <- sta + trend
  
  c(sta,tsta)
}

nonstat_ts <- nonstationarity_sym(ts.len=100,ts.mean=0,ts.var=0.0001,ts.trend=0.0005,seed=12)
data <- data.frame(serie=nonstat_ts,event=FALSE)
data[100, "event"] <- TRUE

require(ggplot2)
#Plotting time series
ggplot(data, aes(y=serie, x=1:nrow(data))) +
  geom_line()+
  xlab("Time")+
  ylab("")+#ylab(serie_name)+
  theme_bw()


run_test <- function(data,lag_pred=0, online_step=1, pred_ic=FALSE, pred_par=NULL, pred_recursive=FALSE, alpha=3,
                     warm_size=30, batch_size=1, mem_batches=31, png_folder="plots_test/", verbose=TRUE){
  
  model <- har_herald(lag_pred=lag_pred, online_step=online_step,
                      pred_ic=pred_ic, pred_par=pred_par, pred_recursive=pred_recursive,
                      alpha=alpha)
  herald_detection <- run_nexus(model,data, warm_size=warm_size, batch_size=batch_size, mem_batches=mem_batches, png_folder=png_folder, verbose=verbose)
  
  if(verbose) 
    cat(paste0("Detection lag:\n",herald_detection$detection_lag$batch," batches\n",
               herald_detection$detection_lag$obs," observations\n",
               herald_detection$detection_lag$time," seconds"))
  
  return(herald_detection)
}

online_step <- 1
alpha <- 3
ma_order <- 3
warm_size <- 30
batch_size <- 1
mem_batches <- 31
png_folder <- "plots_test/"
verbose <- FALSE



# Running tests ----------------------------------------------------------------

#Exp A - No lag, no IC
herald_detection_A <- run_test(data,lag_pred=0, pred_ic=FALSE, pred_par=NULL, pred_recursive=FALSE, png_folder=paste0(png_folder,"ExpA/"), 
                             online_step=online_step, alpha=alpha, warm_size=warm_size, batch_size=batch_size, mem_batches=mem_batches, verbose=verbose)

#Exp B - No lag, MA IC
herald_detection_B <- run_test(data,lag_pred=0, pred_ic=TRUE, pred_par=list(order=ma_order), pred_recursive=FALSE, png_folder=paste0(png_folder,"ExpB/"), 
                             online_step=online_step, alpha=alpha, warm_size=warm_size, batch_size=batch_size, mem_batches=mem_batches, verbose=verbose)

#Exp C - Lag, no IC
herald_detection_C <- run_test(data,lag_pred=1, pred_ic=FALSE, pred_par=NULL, pred_recursive=FALSE, png_folder=paste0(png_folder,"ExpC/"), 
                             online_step=online_step, alpha=alpha, warm_size=warm_size, batch_size=batch_size, mem_batches=mem_batches, verbose=verbose)

#Exp D - Lag, MA IC
herald_detection_D <- run_test(data,lag_pred=1, pred_ic=TRUE, pred_par=list(order=ma_order), pred_recursive=FALSE, png_folder=paste0(png_folder,"ExpD/"), 
                             online_step=online_step, alpha=alpha, warm_size=warm_size, batch_size=batch_size, mem_batches=mem_batches, verbose=verbose)

#Exp E - Lag, recursive, no IC
herald_detection_E <- run_test(data,lag_pred=1, pred_ic=FALSE, pred_par=NULL, pred_recursive=TRUE, png_folder=paste0(png_folder,"ExpE/"), 
                             online_step=online_step, alpha=alpha, warm_size=warm_size, batch_size=batch_size, mem_batches=mem_batches, verbose=verbose)

#Exp F - Lag, recursive, MA IC
herald_detection_F <- run_test(data,lag_pred=1, pred_ic=TRUE, pred_par=list(order=ma_order), pred_recursive=TRUE, png_folder=paste0(png_folder,"ExpF/"), 
                             online_step=online_step, alpha=alpha, warm_size=warm_size, batch_size=batch_size, mem_batches=mem_batches, verbose=verbose)


# Printing lags ----------------------------------------------------------------
print_lags <- function(herald_detection){
  cat(paste0("Detection lag:\n",herald_detection$detection_lag$batch," batches\n",
             herald_detection$detection_lag$obs," observations\n",
             herald_detection$detection_lag$time," seconds"))
}

cat("\nExperiment A: No lag, no IC (Baseline)\n")
print_lags(herald_detection_A)

cat("\nExperiment B: No lag, MA IC\n")
print_lags(herald_detection_B)

cat("\nExperiment C: Lag, no IC\n")
print_lags(herald_detection_C)

cat("\nExperiment D: Lag, MA IC\n")
print_lags(herald_detection_D)

cat("\nExperiment E: Lag, recursive, no IC\n")
print_lags(herald_detection_E)

cat("\nExperiment F: Lag, recursive, MA IC\n")
print_lags(herald_detection_F)