# Model and prediction
herald_mdl <- function(serie, lag=30, ic=FALSE, recursive=FALSE, ...){
  
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
    
    if(lag==0){
      new_mdl <- c(serie_mdl, do.call(pred_func, args=c(list(serie), ...)) )
    }
    else{
      serie_w_pred <- serie
      for(i in 1:lag) serie_w_pred <- c(serie_w_pred, do.call(pred_func, args=c(list(serie_w_pred), ...)))
      
      new_mdl <- c(serie_mdl, tail(serie_w_pred,lag), do.call(pred_func, args=c(list(serie_w_pred), ...)) )
    }
    
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
source("https://raw.githubusercontent.com/cefet-rj-dal/daltoolbox-examples/main/jupyter.R")
load_library("daltoolbox")
source("https://raw.githubusercontent.com/cefet-rj-dal/harbinger-examples/main/jupyter.R")
load_library("daltoolbox")
load_github("cefet-rj-dal/harbinger")


#'@description Ancestor class for time series event detection
#'@details The Harbinger class establishes the basic interface for time series event detection.
#'  Each method should be implemented in a descendant class of Harbinger
#'@return Harbinger object
#'@examples detector <- harbinger()
#'@export
har_herald <- function(pred=TRUE, lag_pred = 30, mdl_train=30, online_step=1,
                       pred_ic=FALSE, pred_par=NULL, pred_recursive=FALSE,
                       conf=0.05, dist_func=function(v1, v2) sqrt(sum((v1 - v2)^2)),
                       mem_dist=30, seq_dist_path=paste0(getwd(),"/seq_dist.rds")) {
  
  obj <- harbinger()
  obj$pred <- pred
  obj$lag_pred <- lag_pred
  obj$mdl_train <- mdl_train
  obj$online_step <- online_step
  obj$n.ahead <- lag_pred + online_step
  obj$pred_ic <- pred_ic
  obj$pred_par <- pred_par
  obj$pred_recursive <- pred_recursive
  obj$conf <- conf
  obj$dist_func <- dist_func
  obj$seq_dist_path <- seq_dist_path
  obj$mem_dist <- mem_dist
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
  herald_mdl <- function(serie, lag=30, ic=FALSE, recursive=FALSE, ...){
    
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
      
      if(lag==0){
        new_mdl <- c(serie_mdl, do.call(pred_func, args=c(list(serie), ...)) )
      }
      else{
        serie_w_pred <- serie
        for(i in 1:lag) serie_w_pred <- c(serie_w_pred, do.call(pred_func, args=c(list(serie_w_pred), ...)))
        
        new_mdl <- c(serie_mdl, tail(serie_w_pred,lag), do.call(pred_func, args=c(list(serie_w_pred), ...)) )
      }
      
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
  detector <- function(data,pred,seq_dist_path,mem_dist,lag_pred,dist_func,conf = 0.05){
    n <- nrow(data)
    non_na <- which(!is.na(data))
    
    serie <- na.omit(data)
    
    #===== Structural change analysis ======
    cp.index <- function(data, conf=0.05, cp_from=0.15){
      require(strucchange)
      ## compute F statistics
      fs <- strucchange::Fstats(data ~ 1, from = cp_from, to = NULL)
      
      ## perform the supF test
      supF_test <- strucchange::sctest(fs, type="supF")
      
      if(supF_test$p.value < conf) {
        ## F statistics indicate one breakpoint
        bp <- breakpoints(fs)
        bp_point <- bp$breakpoints
        return(bp_point)
      }
      else return(NA)
    }
    

    if(!is.null(pred)){
      #Sequence distance between prediction and serie
      d <- dist_func( tail(serie,lag_pred+1), tail(pred,lag_pred+1) )
      
      if(file.exists(seq_dist_path)) seq_dist <- readRDS(seq_dist_path)
      else seq_dist <- NULL
      
      distr_dist <- seq_dist
      
      if(length(distr_dist) < mem_dist){
        changepoints <- rep(FALSE, n)
      }
      else{
        distr_dist <- c(seq_dist,d)
        
        cp <- cp.index(distr_dist, conf)
        changepoints <- rep(FALSE, length(distr_dist))
        if(!is.na(cp))  changepoints[cp] <- TRUE
        
        lookback_distr <- min(c( n, length(distr_dist) ))
        changepoints <- c(rep(FALSE, n-lookback_distr),tail(changepoints,lookback_distr))
      }
      
      distr_dist <- c(seq_dist,d)
      seq_dist <- tail(distr_dist, mem_dist)
      saveRDS(seq_dist,file=seq_dist_path)
    }
    else{
      if(n < mem_dist){
        changepoints <- rep(FALSE, n)
      }
      else{
        cp <- cp.index(serie[[1]], conf)
        changepoints <- rep(FALSE, n)
        if(!is.na(cp))  changepoints[cp] <- TRUE
      }
    }

    i_cps <- rep(NA, n)
    i_cps[non_na] <- changepoints 
    
    detection <- data.frame(idx=1:n, event = i_cps, type="")
    detection$type[i_cps] <- "changepoint"
    
    return(detection)
  }
  
  if(obj$pred){
    mdl_serie <- herald_mdl(serie, lag=obj$lag_pred, ic=obj$pred_ic, recursive=obj$pred_recursive, obj$pred_par)
    mdl_serie <- as.data.frame(mdl_serie)
  }
  else mdl_serie <- NULL
  
  serie <- as.data.frame(serie)
  
  #detection over predicted inertial component
  detection <- detector(serie, mdl_serie, seq_dist_path=obj$seq_dist_path, mem_dist=obj$mem_dist, lag_pred=obj$lag_pred, dist_func=obj$dist_func, conf=obj$conf)
  
  return(detection)
}


# Herald Nexus function ---------------------------------------------------------------

# Loading Nexus
#source("https://raw.githubusercontent.com/cefet-rj-dal-dev/nexus/main/R/nexus.R?token=GHSAT0AAAAAACEGUDHGSHHPJXOZH2FHDJ42ZFFZRXA")
source("dev/nexus/nexus.R")

herald_nexus <- function(model,data, full=FALSE, png_folder="dev/plots/",title=TRUE,exp_descr="", verbose=TRUE){
  
  
  batch_size <- model$online_step
  
  if(model$pred){
    mem_batches <- model$mdl_train + model$lag_pred + 1
    warm_size <- mem_batches
    if(full) mem_batches <- 0
  }
  else{
    mem_batches <- model$mem_dist + 1
    warm_size <- mem_batches
    if(full) mem_batches <- 0
  }
  
  datasource <- nex_simulated_datasource("data", data$serie)
  
  online_detector <- nexus(datasource, model, warm_size = warm_size, batch_size = batch_size, mem_batches = mem_batches)
  
  online_detector <- warmup(online_detector)
  
  bt <- 1
  event_happened <- FALSE
  already_detected <- FALSE
  bt_event_happened <- 0
  event_idx <- 0
  
  while (!is.null(online_detector$datasource)) {
    
    online_detector <- detect(online_detector)
    
    #if(verbose) print(paste0("Batch ",bt))

    #first event happened
    if(any(data$event[online_detector$detection$idx]) & !event_happened){
      if(verbose) print("Event happened!")
      event_happened <- TRUE
      event_idx <- which(data$event[online_detector$detection$idx])
      bt_event_happened <- bt
      time_event_happened <- Sys.time()
    }
    
    #first detection after event happened
    if(any(which(as.logical(online_detector$detection$event))>event_idx) & event_happened) {
      
      if(!already_detected){
        if(verbose) print("Detected event!")
        online_detector$detection_lag <- list()
        online_detector$detection_lag$batch <- bt-bt_event_happened
        det.index <- which(as.logical(online_detector$detection$event))
        online_detector$detection_lag$obs <- min(det.index[det.index>event_idx])-event_idx
        online_detector$detection_lag$time <- round(as.numeric(difftime(time1 = Sys.time(), time2 = time_event_happened, units = "secs")), 3)
        already_detected <- TRUE
      }
      
      #for plotting (move out of this if for gif generation)
      #exp_descr <- paste0(exp_descr,", batch:",bt)
      
      modifyString <- function(inputString) {
        modifiedString <- gsub(":", "", inputString)  # Remove ":"
        modifiedString <- gsub(", ", "_", modifiedString)  # Substitute "," with "_"
        return(modifiedString)
      }
      png_file <- paste0(png_folder,modifyString(exp_descr),".png")
      #png(file=png_file,width=600, height=350)
      
      if(title) {
        det_lag <- paste0("Detection lag: ",online_detector$detection_lag$batch," batches, ",
                          online_detector$detection_lag$obs," observations, ",
                          online_detector$detection_lag$time," seconds")
        plot_title <- paste0("Settings: ",exp_descr,"\n",det_lag,"\n")
      }
      
      
      grf <- 
        har_plot(online_detector$detector, data$serie[online_detector$detection$idx], online_detector$detection, data$event[online_detector$detection$idx])
        
      if(model$pred){
        herald_serie <- herald_mdl(online_detector$serie, lag=model$lag_pred, ic=model$pred_ic, recursive=model$pred_recursive, model$pred_par)  
        grf <- grf+
          geom_line(aes(x=online_detector$detection$idx,y=c(rep(NA,length(online_detector$detection$idx)-length(herald_serie)),herald_serie)),col="orange")+
          geom_point(aes(x=online_detector$detection$idx,y=c(rep(NA,length(online_detector$detection$idx)-model$n.ahead),tail(herald_serie,model$n.ahead))),col="orange")
      }
        
      require(ggtext)
      if(title) grf <- grf + ggtitle(plot_title) + theme(plot.title = element_textbox_simple())
      
      #plot(grf)
      #dev.off()
      #plot(grf)
      online_detector$plot<- grf
      
      break #comment this line for full result 
    }
    
    bt <- bt + 1
  }
  
  return(online_detector)
}


# Experiment function --------------------------------------------------------
run_exp <- function(data,id,pred=TRUE,lag_pred=0, mdl_train=30, online_step=1, pred_ic=FALSE, pred_par=NULL, pred_recursive=FALSE, 
                    dist_func=function(v1, v2) sqrt(sum((v1 - v2)^2)), conf=0.05, full=FALSE,
                    mem_dist=30, seq_dist_path=paste0(getwd(),"/seq_dist.rds"), res_path=paste0(getwd(),"/res_exp.rds"),
                    png_folder="plots/", title=TRUE, verbose=TRUE){
  
  #harbinger model definition
  model <- har_herald(pred=pred,lag_pred=lag_pred, mdl_train=mdl_train, online_step=online_step,
                      pred_ic=pred_ic, pred_par=pred_par, pred_recursive=pred_recursive,
                      conf=conf, dist_func=dist_func, mem_dist=mem_dist, seq_dist=seq_dist_path)
  
  #Experimental settings description 
  exp_settings <- paste0("id:",id,", pred:",pred,", lag:",lag_pred,", ic:",pred_ic,", recursive:",pred_recursive,", mem:",mem_dist,
                         ", full:",full,", mdl:",mdl_train,", conf:",conf)
  if(pred_ic) exp_settings <- paste0(exp_settings,", par:",pred_par[[1]])
  
  #running online detection
  herald_detection <- herald_nexus(model,data,title=title,exp_descr=exp_settings,full=full, png_folder=png_folder, verbose=verbose)
  
  if(!is.null(herald_detection$detection_lag)){
    detlag_batch <- herald_detection$detection_lag$batch
    detlag_obs <- herald_detection$detection_lag$obs
    detlag_time <- herald_detection$detection_lag$time
  }
  else detlag_batch <- detlag_obs <- detlag_time <- NA
  
  if(verbose) 
    cat(paste0("Detection lag:\n",detlag_batch," batches\n",
               detlag_obs," observations\n",
               detlag_time," seconds"))
  
  #Writing results into data.frame
  if(file.exists(res_path)) results <- readRDS(res_path)
  else results <- data.frame()
  
  exp_res <- cbind(id=id,pred=pred,lag=lag_pred,mdl=mdl_train,ic=pred_ic,par=ifelse(pred_ic,pred_par[[1]],NA),
                   recursive=pred_recursive,mem=mem_dist,conf=conf,
                   detlag_batch=detlag_batch,detlag_obs=detlag_obs,detlag_time=detlag_time)
  results <- rbind(results,exp_res)
  saveRDS(results,file=res_path)
  
  return(list(detection=herald_detection,results=results))
}




#Experiments ===================================================================

# Running experiments function----------------------------------------------------------------
experiments <- function(data,parameters,verbose=TRUE){
  
  settings <- expand.grid(parameters[lengths(parameters) > 1])
  par_fixed <- parameters[lengths(parameters) == 1]
  
  getpar <- function(v1,v2,name){
    if(name %in% names(v1)) v1[[name]]
    else v2[[name]]
  }
  
  det_results <- list()
  p <- par_fixed
  res_path <- getpar(settings,p,"res_path")
  if(file.exists(res_path)) file.remove(res_path) #Reset memory of test results

  nexp <- nrow(settings)
  if(nrow(settings)==0) nexp <- 1
  
  for(setting in 1:nexp){
    
    if(verbose) cat("Exp ",setting,"/",nrow(settings),"\n")
    s <- settings[setting,, drop = FALSE]
    
    if(getpar(s,p,"pred_ic")) pred_par <- list(order=getpar(s,p,"ma_order"))
    else pred_par <- NULL
    
    seq_dist_path <- getpar(s,p,"seq_dist_path")
    if(file.exists(seq_dist_path)) file.remove(seq_dist_path) #Reset memory of sequence distances
    
    det_results[[setting]] <- run_exp(data,id=setting, pred=getpar(s,p,"pred"), mdl_train=getpar(s,p,"mdl_train"), lag_pred=getpar(s,p,"lag_pred"), pred_ic=getpar(s,p,"pred_ic"), 
                                      pred_par=pred_par, pred_recursive=getpar(s,p,"pred_recursive"), 
                                      conf=getpar(s,p,"conf"), mem_dist=getpar(s,p,"mem_dist"), full=getpar(s,p,"full_mem"), 
                                      
                                      dist_func=getpar(s,p,"dist_func"), png_folder=getpar(s,p,"png_folder"), title=getpar(s,p,"plot_title"), 
                                      online_step=getpar(s,p,"online_step"), seq_dist_path=seq_dist_path, 
                                      verbose=getpar(s,p,"verbose"), res_path=res_path) 
  }

  return(list(detection=det_results,results=det_results[[setting]]$results))
}


#Baseline method based on the Chow test-----------------------------------------
#Create my own version of online -> BASELINE!
baseline_cp <- function(data, mem=27, conf=0.05, cp_from=0.15, full=FALSE){
  
  bt <- 1
  event_happened <- FALSE
  already_detected <- FALSE
  bt_event_happened <- 0
  event_idx <- 0
  break_detected <- FALSE
  
  bps <- NULL
  for(i in mem:length(serie)){
    
    begin_window <- ifelse(full,1,(i-mem+1))
    mem_data <- data[begin_window:i,]
    mem_serie <- mem_data$serie
    
    #first event happened
    if(any(mem_data$event) & !event_happened){ #Event happened!
      event_happened <- TRUE
      event_idx <- min(which(data$event))
      bt_event_happened <- bt
      time_event_happened <- Sys.time()
    }
    
    ## compute F statistics
    fs <- Fstats(mem_serie ~ 1, from = cp_from, to = NULL)
    ## perform the supF test
    supF_test <- sctest(fs, type="supF")
    if(supF_test$p.value < conf) break_detected <- TRUE
    
    #first detection after event happened
    if(break_detected & event_happened) {
      
      ## F statistics indicate one breakpoint
      bp <- breakpoints(fs)
      bp_point <- bp$breakpoints+begin_window-1
      bps <- c(bps, bp_point)
      #cat(paste0("Caught breakpoint at ",bp_point," in time ",i,"\n"))
      
      if(!already_detected){
        detection_lag <- list()
        detection_lag$batch <- bt-bt_event_happened
        detection_lag$obs <- bp_point-event_idx
        detection_lag$time <- round(as.numeric(difftime(time1 = Sys.time(), time2 = time_event_happened, units = "secs")), 3)
        already_detected <- TRUE
      }
      
      break #comment this line for full result
    }
    break_detected <- FALSE
    bt <- bt + 1
  }
  
  detection <- rep(FALSE,nrow(data))
  detection[bps] <- TRUE
  
  return(list(detection=detection,lag=detection_lag))
}


# Simulated time series function --------------------------------------------------------

# Function to generate simulated time series with different trend functions
trend_sym <- function(n=100,ts.var=0.1,seed=12345,plot=TRUE,
                      trend_function=c("linear","log","sqrt","poly2","poly3","exp","sigmoid")) {
  
  trend_function <- match.arg(trend_function)
  
  if (trend_function == "linear") {
    x <- tail(seq(from = 1, to = 2, length.out=n+1),n)
    trend <- x
  } else if (trend_function == "log") {
    x <- tail(seq(from = exp(1), to = exp(2), length.out=n+1),n)
    trend <- log(x)
  } else if (trend_function == "sqrt") {
    x <- tail(seq(from = 1, to = 2^2, length.out=n+1),n)
    trend <- sqrt(x)
  } else if (trend_function == "poly2") {
    x <- tail(seq(from = 1, to = sqrt(2), length.out=n+1),n)
    trend <- (x)^2
  } else if (trend_function == "poly3") {
    x <- tail(seq(from = 1, to = 2^(1/3), length.out=n+1),n)
    trend <- (x)^3
  } else if (trend_function == "exp") {
    x <- tail(seq(from = log(1), to = log(2), length.out=n+1),n)
    trend <- exp(x)
  } else if (trend_function == "sigmoid") {
    x <- 1:n
    sigmoid <- function(x, midpoint, slope) {
      y <- 1 / (1 + exp(-slope * (x - midpoint)))
      y <- pmin(pmax(y, 0), 1)  # Clip values between 0 and 1
      return(y)
    }
    midpoint <- n / 2
    slope <- 0.1*(2*(100/n))
    trend <- sigmoid(x, midpoint, slope)
  } else {
    stop("Invalid trend function specified.")
  }
  
  # Combine the horizontal line and trend
  if(trend_function == "sigmoid") time_series <- c(rep(0, n), trend)
  else time_series <- c(rep(1, n), trend)
  
  # Add volatility by introducing random noise
  set.seed(seed)
  time_series <- time_series + rnorm(n, sd = ts.var)
  
  data <- data.frame(serie=time_series,event=FALSE)
  data[n, "event"] <- TRUE
  
  if(plot) print(ggplot(data, aes(y=serie, x=1:nrow(data))) +
    geom_line()+xlab("Time")+ ylab("")+theme_bw())
  
  return(data)
}



level_sym <- function(ts.len,ts.mean=0,ts.var=1,increase.level=5,seed=12345){
  require(tseries)
  require(forecast)
  
  #x variable (time)
  t <- c(1:ts.len)
  
  #stationary time series
  set.seed(1234)
  sta <- arima.sim(model=list(ar=0.5), n=ts.len, mean=ts.mean, sd=sqrt(ts.var)) #AR(1)
  
  #test of stationarity OK!
  #adf.test(sta)
  
  #level-stationary time series
  level <- c(rep(ts.mean,ts.len/2),rep(ts.mean+increase.level,ts.len/2))
  lsta <- sta + level
  
  nonstat_ts <- c(sta,tsta)
  data <- data.frame(serie=nonstat_ts,event=FALSE)
  data[ts.len+ts.len/2, "event"] <- TRUE
  
  return(data)
}


volatility_sym <- function(ts.len,ts.mean=0,ts.var=1,increase.var=2,seed=12345){
  require(tseries)
  require(forecast)
  
  #x variable (time)
  t <- c(1:ts.len)
  
  #stationary time series
  set.seed(1234)
  sta <- arima.sim(model=list(ar=0.5), n=ts.len, mean=ts.mean, sd=sqrt(ts.var)) #AR(1)
  
  #test of stationarity OK!
  #adf.test(sta)
  
  #heteroscedastic time series (nonstationary in variance)
  var.level <- c(rep(1,ts.len/2),rep(increase.var,ts.len/2))
  hsta <- sta * var.level
  
  nonstat_ts <- c(sta,hsta)
  data <- data.frame(serie=nonstat_ts,event=FALSE)
  data[ts.len+ts.len/2, "event"] <- TRUE
  
  return(data)
}



#Experiment results
# Generate the simulated time series
# Specify the trend function ("linear", "log", "poly2", "poly3", "exp", "sigmoid")
#data <- trend_sym(n=100,ts.var=0.001,seed=12345, trend_function="sigmoid")
#herald_results <- experiments(data,parameters)
