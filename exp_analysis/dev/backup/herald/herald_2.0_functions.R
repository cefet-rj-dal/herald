herald_mdl <- function(serie, mdl_func, mdl_par=NULL, pred_func, lag=30, recursive=FALSE, ...){
  
  lag_pred <- function(serie, mdl_func, mdl_parr=NULL, pred_func, lag, ...){
    
    if(length(serie)-lag-1<=0) lag <- 0
    serie <- serie[1:(length(serie)-lag-1)]
    serie_mdl <- do.call(mdl_func, args=c(list(serie), mdl_par))
    
    serie_w_pred <- serie
    serie_w_pred_mdl <- serie_mdl
    for(i in 1:lag) {
      serie_w_pred <- c(serie_w_pred, do.call(pred_func, args=c(list(serie_w_pred),list(serie_w_pred_mdl), ...)))
      serie_w_pred_mdl <- do.call(mdl_func, args=c(list(serie_w_pred), mdl_par))
    }
    
    new_mdl <- c(fitted(serie_mdl), tail(serie_w_pred,lag), do.call(pred_func, args=c(list(serie_w_pred),list(serie_w_pred_mdl), ...)) )
    
    return(new_mdl)
  }
  
  model_lag_pred <- function(serie, mdl_func, mdl_parr=NULL, pred_func, lag, ...){
    mdl <- list()
    for(i in length(serie):1){
      mdl[i] <- tail( do.call(lag_pred, args=c(list(serie[1:i]), mdl_func=mdl_func, mdl_par=mdl_par, 
                                               pred_func=pred_func, lag=lag, ...)) ,1)
    }
    return(unlist(mdl))
  }
  
  if(recursive) return(model_lag_pred(serie, mdl_func, mdl_par, pred_func, lag=lag,...))
  else return(lag_pred(serie, mdl_func, mdl_par, pred_func, lag=lag,...))
  
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
har_herald <- function(lag_pred = 30, mdl_train=30, online_step=1,
                       pred_ic=FALSE, ic_func=NULL, ic_par=NULL, mdl_func=NULL, mdl_par=NULL,
                       pred_func=NULL, pred_par=NULL, pred_recursive=FALSE,
                       lookback=FALSE, alpha=1.5, dist_func=function(v1, v2) sqrt(sum((v1 - v2)^2)),
                       mem_dist=30, seq_dist_path=paste0(getwd(),"/seq_dist.rds")) {
  obj <- harbinger()
  obj$lag_pred <- lag_pred
  obj$mdl_train <- mdl_train
  obj$online_step <- online_step
  obj$n.ahead <- lag_pred + online_step
  obj$pred_ic <- pred_ic
  obj$ic_func <- ifelse(is.null(ic_func),function(serie, order) forecast::ma(serie, order, centre=FALSE),ic_func)
  obj$ic_par <- ic_par
  obj$mdl_func <- ifelse(is.null(mdl_func),function(serie) forecast::auto.arima(serie),mdl_func)
  obj$mdl_par <- mdl_par
  obj$pred_func <- ifelse(is.null(pred_func),function(serie, mdl) forecast::forecast(serie, model=mdl, h=1)$mean,pred_func)
  obj$pred_par <- pred_par
  obj$pred_recursive <- pred_recursive
  obj$lookback <- lookback
  obj$alpha <- alpha
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
  input_data <- serie
  
  # Model and prediction
  herald_mdl <- function(serie, mdl_func, mdl_par=NULL, pred_func, lag=30, recursive=FALSE, ...){
    
    lag_pred <- function(serie, mdl_func, mdl_par=NULL, pred_func, lag, ...){
      
      if(length(serie)-lag-1<=0) lag <- 0
      serie <- serie[1:(length(serie)-lag-1)]
      serie_mdl <- do.call(mdl_func, args=c(list(serie), mdl_par))
      
      if(lag==0){
        new_mdl <- c(fitted(serie_mdl), do.call(pred_func, args=c(list(serie),list(serie_mdl), ...)) )
      }
      else{
        serie_w_pred <- serie
        serie_w_pred_mdl <- serie_mdl
        
        for(i in 1:lag) {
          serie_w_pred <- c(serie_w_pred, do.call(pred_func, args=c(list(serie_w_pred),list(serie_w_pred_mdl), ...)))
          serie_w_pred_mdl <- do.call(mdl_func, args=c(list(serie_w_pred), mdl_par))
        }
        
        new_mdl <- c(fitted(serie_mdl), tail(serie_w_pred,lag), do.call(pred_func, args=c(list(serie_w_pred),list(serie_w_pred_mdl), ...)) )
      }
      
      return(new_mdl)
    }
    
    model_lag_pred <- function(serie, mdl_func, mdl_par=NULL, pred_func, lag, ...){
      mdl <- list()
      for(i in length(serie):1){
        mdl[i] <- tail( do.call(lag_pred, args=c(list(serie[1:i]), mdl_func=mdl_func, mdl_par=mdl_par, 
                                                 pred_func=pred_func, lag=lag, ...)) ,1)
      }
      return(unlist(mdl))
    }
    
    if(recursive) return(model_lag_pred(serie, mdl_func, mdl_par, pred_func, lag=lag,...))
    else return(lag_pred(serie, mdl_func, mdl_par, pred_func, lag=lag,...))
    
  }
  
  # Detection
  detector <- function(serie,pred,seq_dist_path,mem_dist,lag_pred,dist_func,lookback=FALSE,alpha = 1.5){
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
      cond = data < lq1 | data > hq3
      #cond = data > hq3
      return(list(cond=cond,lq1=lq1,hq3=hq3))
    }
    
    outliers.boxplot.index <- function(data, alpha = 1.5){
      cond <- outliers.boxplot(data, alpha)$cond
      index.cp = which(cond)
      return (index.cp)
    }
    
    #Sequence distance between prediction and serie
    d <- dist_func( tail(serie,lag_pred+1), tail(pred,lag_pred+1) )
    
    if(file.exists(seq_dist_path)) seq_dist <- readRDS(seq_dist_path)
    else seq_dist <- NULL
    
    distr_dist <- seq_dist
    
    if(length(distr_dist) < mem_dist){
      outliers <- rep(FALSE, n)
    }
    else{
      if(!lookback) { #test only the last d against the boxplot of distr_dist
        boxplot <- outliers.boxplot(distr_dist, alpha)
        d_is_outlier <- d < boxplot$lq1 | d > boxplot$hq3
        outliers <- c(rep(FALSE, n-1),d_is_outlier)
      }
      else{
        distr_dist <- c(seq_dist,d)
        
        outliers <- outliers.boxplot.index(distr_dist, alpha)
        group_outliers <- split(outliers, cumsum(c(1, diff(outliers) != 1)))
        outliers <- rep(FALSE, length(distr_dist))
        for (g in group_outliers) {
          if (length(g) > 0) {
            i <- min(g)
            outliers[i] <- TRUE
          }
        }
        
        lookback_distr <- min(c( n, length(distr_dist) ))
        outliers <- c(rep(FALSE, n-lookback_distr),tail(outliers,lookback_distr))
      }
    }

    i_outliers <- rep(NA, n)
    i_outliers[non_na] <- outliers 
    
    detection <- data.frame(idx=1:n, event = i_outliers, type="")
    detection$type[i_outliers] <- "change point"
    
    distr_dist <- c(seq_dist,d)
    seq_dist <- tail(distr_dist, mem_dist)
    saveRDS(seq_dist,file=seq_dist_path)
    
    return(detection)
  }
  
  #inertial component
  if(obj$pred_ic) {
    ic <- do.call(obj$ic_func, args=c(list(serie), obj$ic_par))
    ic <- na.omit(ic)
    
    input_data <- ic
  }
  
  mdl_serie <- herald_mdl(input_data, mdl_func=obj$mdl_func, mdl_par=obj$mdl_par, pred_func=obj$pred_func, 
                          lag=obj$lag_pred, recursive=obj$pred_recursive, obj$pred_par)
  
  serie <- as.data.frame(serie)
  mdl_serie <- as.data.frame(mdl_serie)
  
  #detection over predicted inertial component
  detection <- detector(serie, mdl_serie, seq_dist_path=obj$seq_dist_path, mem_dist=obj$mem_dist, lag_pred=obj$lag_pred, dist_func=obj$dist_func, lookback=obj$lookback, alpha=obj$alpha)
  
  return(detection)
}


# Herald Nexus function ---------------------------------------------------------------

# Loading Nexus
#source("https://raw.githubusercontent.com/cefet-rj-dal-dev/nexus/main/R/nexus.R?token=GHSAT0AAAAAACEGUDHGSHHPJXOZH2FHDJ42ZFFZRXA")
source("dev/nexus.R")

herald_nexus <- function(model,data, full=FALSE, png_folder="dev/plots/",title=TRUE,exp_descr="", verbose=TRUE){
  
  mem_batches <- model$mdl_train + model$lag_pred +1
  warm_size <- mem_batches
  batch_size <- model$online_step
  
  if(full) mem_batches <- 0
  
  datasource <- nex_simulated_datasource("data", data$serie)
  
  online_detector <- nexus(datasource, model, warm_size = warm_size, batch_size = batch_size, mem_batches = mem_batches)
  
  online_detector <- warmup(online_detector)
  
  bt <- 1
  event_happened <- FALSE
  bt_event_happened <- 0
  event_idx <- 0
  
  while (!is.null(online_detector$datasource)) {
    
    online_detector <- detect(online_detector)
    
    herald_serie <- herald_mdl(online_detector$serie, mdl_func=model$mdl_func, mdl_par=model$mdl_par, pred_func=model$pred_func, 
                               lag=model$lag_pred, recursive=model$pred_recursive, model$pred_par)
    
      
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
      
      if(verbose) print("Detected event!")
      online_detector$detection_lag <- list()
      online_detector$detection_lag$batch <- bt-bt_event_happened
      online_detector$detection_lag$obs <- min(which(as.logical(online_detector$detection$event)))-event_idx
      online_detector$detection_lag$time <- round(as.numeric(difftime(time1 = Sys.time(), time2 = time_event_happened, units = "secs")), 3)
      
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
        har_plot(online_detector$detector, data$serie[online_detector$detection$idx], online_detector$detection, data$event[online_detector$detection$idx])+
        geom_line(aes(x=online_detector$detection$idx,y=c(rep(NA,length(online_detector$detection$idx)-length(herald_serie)),herald_serie)),col="orange")+
        geom_point(aes(x=online_detector$detection$idx,y=c(rep(NA,length(online_detector$detection$idx)-model$n.ahead),tail(herald_serie,model$n.ahead))),col="orange")
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
run_exp <- function(data,id,lag_pred=0, mdl_train=30, online_step=1, pred_ic=FALSE, ic_func=NULL, ic_par=NULL, 
                    mdl_func=NULL, mdl_par=NULL, pred_func=NULL, pred_par=NULL, pred_recursive=FALSE, 
                    dist_func=function(v1, v2) sqrt(sum((v1 - v2)^2)), lookback=FALSE, alpha=3, full=FALSE,
                    mem_dist=30, seq_dist_path=paste0(getwd(),"/seq_dist.rds"), res_path=paste0(getwd(),"/res_exp.rds"),
                    png_folder="plots_test/", title=TRUE, verbose=TRUE){
  
  #harbinger model definition
  model <- har_herald(lag_pred=lag_pred, mdl_train=mdl_train, online_step=online_step,
                      pred_ic=pred_ic, ic_func=ic_func, ic_par=ic_par, mdl_func=mdl_func, mdl_par=mdl_par,
                      pred_func=pred_func, pred_par=pred_par, pred_recursive=pred_recursive,
                      lookback=lookback,alpha=alpha, dist_func=dist_func, mem_dist=mem_dist, seq_dist=seq_dist_path)
  
  
  #Experimental settings description 
  exp_settings <- paste0("id:",id,", lag:",lag_pred,", ic:",pred_ic,", recursive:",pred_recursive,", mem:",mem_dist,", lookback:",lookback,
                         ", full:",full,", mdl:",mdl_train,", alpha:",alpha)
  if(pred_ic) exp_settings <- paste0(exp_settings,", par:",pred_par[[1]])
  
  
  #Emptying memory of distances
  if(file.exists(seq_dist_path)) file.remove(seq_dist_path)
  
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
  
  exp_res <- cbind(id=id,lag=lag_pred,mdl=mdl_train,ic=pred_ic,par=ifelse(pred_ic,ic_par[[1]],NA),
                   recursive=pred_recursive,mem=mem_dist,alpha=alpha,lookback=lookback,
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
  if(file.exists(getpar(settings,p,"res_path"))) file.remove(getpar(settings,p,"res_path"))
  
  nexp <- nrow(settings)
  if(nrow(settings)==0) nexp <- 1
    
  for(setting in 1:nexp){
    
    if(verbose) cat("Exp ",setting,"/",nexp,"\n")
    s <- settings[setting,, drop = FALSE]
    
    if(getpar(s,p,"pred_ic")) ic_par <- list(order=getpar(s,p,"ma_order"))
    else ic_par <- NULL
    
    det_results[[setting]] <- run_exp(data,id=setting, mdl_train=getpar(s,p,"mdl_train"), lag_pred=getpar(s,p,"lag_pred"), pred_ic=getpar(s,p,"pred_ic"), 
                                      ic_func=getpar(s,p,"ic_func"), ic_par=ic_par, mdl_func=getpar(s,p,"mdl_func"), mdl_par=getpar(s,p,"mdl_par"), 
                                      pred_func=getpar(s,p,"pred_func"), pred_par=getpar(s,p,"pred_par"), pred_recursive=getpar(s,p,"pred_recursive"), 
                                      alpha=getpar(s,p,"alpha"), mem_dist=getpar(s,p,"mem_dist"), lookback=getpar(s,p,"lookback"), full=getpar(s,p,"full_mem"), 
                                      
                                      dist_func=getpar(s,p,"dist_func"), png_folder=getpar(s,p,"png_folder"), title=getpar(s,p,"plot_title"), 
                                      online_step=getpar(s,p,"online_step"), seq_dist_path=getpar(s,p,"seq_dist_path"), 
                                      verbose=getpar(s,p,"verbose"), res_path=getpar(s,p,"res_path"))
  }

  return(list(detection=det_results,results=det_results[[setting]]$results))
}

# General experimental settings --------------------------------------------------------

# parameters <- list(
#   #Herald hiperparameters
#   mdl_train = c(27,81,243), #m
#   lag_pred = c(0,1,2,3,9,27), #q = 0,1,2 (anom); 3,9,27 (cp)
#   mem_dist = c(9,27,81), #k
#   ma_order = c(3,9,27), #o
#   pred_ic = c(FALSE,TRUE),
#   pred_recursive = c(FALSE,TRUE),
#   alpha = c(1.5,3), #a
#   lookback = c(FALSE,TRUE), #look at past outliers in memory?
# 
#   #Herald hiperparameters (fixed)
#   online_step = 1,
#   dist_func = function(v1, v2) sqrt(sum((v1 - v2)^2)), #euclidean
#   seq_dist_path = paste0(getwd(),"/dev/seq_dist.rds"),
# 
#   #plot hiperparameters (fixed)
#   png_folder = paste0(getwd(),"/dev/plots_test/"),
#   verbose = FALSE,
#   plot_title = TRUE,
# 
#   #experiment results file
#   res_path = paste0(getwd(),"/dev/res_exp/res_exp.rds")
# )



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
