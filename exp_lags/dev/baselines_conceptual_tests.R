#############################################################################
## Integration between Nexus + Harbinger + DAL Events
#############################################################################
# install.packages("devtools")
source("https://raw.githubusercontent.com/cefet-rj-dal/daltoolbox-examples/main/jupyter.R")
load_library("daltoolbox")
source("https://raw.githubusercontent.com/cefet-rj-dal/harbinger-examples/main/jupyter.R")
load_library("daltoolbox")
load_github("cefet-rj-dal/harbinger")


# Nexus -------------------------------------------------------------------
# Run Nexus ---------------------------------------------------------------
run_nexus <- function(model, data, warm_size = 30, batch_size = 30, mem_batches = 0, png_folder="dev/plots/") {
  #Create auxiliary batch and slide counters
  bt_num <- 1
  sld_bt <- 1
  ef_start <- FALSE
  
  event_happened <- FALSE
  bt_event_happened <- 0
  event_idx <- 0
  
  #Prepare data to experiment
  datasource <- nex_simulated_datasource("data", data$serie)
  online_detector <- nexus(datasource, model, warm_size = warm_size, batch_size = batch_size, mem_batches = mem_batches)
  online_detector <- warmup(online_detector)
  
  #Sliding batches through series
  while (!is.null(online_detector$datasource)) {
    online_detector <- detect(online_detector)
    
    
    #first event happened
    if(any(data$event[online_detector$detection$idx]) & !event_happened){
      event_happened <- TRUE
      event_idx <- which(data$event[online_detector$detection$idx])
      bt_event_happened <- bt_num
      time_event_happened <- Sys.time()
    }
    
    
    #Update batch and slide counters
    print(paste("Current position:", sld_bt+warm_size))
    sld_bt <- sld_bt + 1
    

    if(event_happened &
       any(as.logical(online_detector$detection$event) &
           seq_along(online_detector$detection$event) > event_idx &
           !online_detector$detection$type %in% c("","anomaly") ) ) {
      
      online_detector$detection_lag <- list()
      online_detector$detection_lag$batch <- bt_num-bt_event_happened
      online_detector$detection_lag$obs <- min(which(as.logical(online_detector$detection$event) &
                                                       seq_along(online_detector$detection$event) > event_idx &
                                                       !online_detector$detection$type %in% c("","anomaly") ))-event_idx
      online_detector$detection_lag$time <- round(as.numeric(difftime(time1 = Sys.time(), time2 = time_event_happened, units = "secs")), 3)
      break #comment this line for full result 
    }
    
    #Batch number update
    bt_num <- bt_num + 1
    
    #Print partial results
    print("Results:")
    print(table(online_detector$detection$event))
    print("--------------------------")
    print(paste("Batch:", bt_num))
    print("==========================")
    
  }
  
  return(online_detector)
}


lag_pred <- 27 #c(0,1,2,3,9,27) #q = 0,1,2 (anom); 3,9,27 (cp)
mdl_train <- 27 #c(27,81,243), #m
online_step <- 1
ma_order <- 18 #max sem dar erro #c(3,9,27), #o
alpha <- 1.5 #fixo no Harbinger #c(1.5,3), #a
mem_batches <- mdl_train + lag_pred + 1
warm_size <- mem_batches
batch_size <- online_step



## --------------------------------------------------------
# Run Experiments

# Generate the simulated time series
# Specify the trend function ("linear", "log", "sqrt", "poly2", "poly3", "exp", "sigmoid")
source("dev/herald_conceptual_tests_functions.R")
library(ggplot2)
library(dplyr)
n <- 100
seed <- 123
volatility <- c(low=0.025,medium=0.05,high=0.1) #volatility levels
exp_series <- sapply(c("linear","log","sqrt","poly2","poly3","exp","sigmoid"), function(curve) {
  sapply(names(volatility), function(v)
    trend_sym(n, ts.var=volatility[v], seed=seed, trend_function=curve),simplify = FALSE, USE.NAMES = TRUE)
},simplify = FALSE, USE.NAMES = TRUE
)
test_series <- exp_series[1]

exp_res_CF <- list()
exp_res_CF_df <- data.frame()
for(curve in names(exp_series)){
  exp_res_CF[[curve]] <- list()
  
  for(var in names(exp_series[[curve]])){
    data <- exp_series[[curve]][[var]]
    name_data <- paste0(curve,"_",var)
    
    # establishing method
    model <- hcp_cf_arima(sw_size = ma_order) #CF using Arima
    result <- run_nexus(model=model, data=data, warm_size=warm_size, batch_size=batch_size, mem_batches=mem_batches, png_folder="dev/plots_baseline/")
    
    # plotting the results
    grf <- har_plot(result$detector, head(data$serie,nrow(result$detection)), result$detection, head(data$event,nrow(result$detection)))
    exp_res_CF[[curve]][[var]] <- list(res=result, plot=grf)
    
    res_df <- cbind(curve=curve,var=var,lag=NA,mdl=NA,ic=TRUE,par=ma_order,
                     recursive=NA,mem=NA,alpha=alpha,lookback=TRUE,
                     detlag_batch=result$detection_lag$batch,
                     detlag_obs=result$detection_lag$obs,
                     detlag_time=result$detection_lag$time,exp="CF")
    exp_res_CF_df <- rbind(exp_res_CF_df,res_df)
    
    folder <- paste0(getwd(),"/dev/plots_test/",name_data,"/")
    if (!dir.exists(folder)) dir.create(folder)
    saveRDS(list(res=exp_res_CF,df_res=exp_res_CF_df),paste0(folder,"res_exp_CF.rds"))
  }
}


exp_res_SCP <- list()
exp_res_SCP_df <- data.frame()
for(curve in names(exp_series)){
  exp_res_SCP[[curve]] <- list()
  
  for(var in names(exp_series[[curve]])){
    data <- exp_series[[curve]][[var]]
    name_data <- paste0(curve,"_",var)
    
    # establishing method
    model <- hcp_scp(sw = 50) #SCP
    result <- run_nexus(model=model, data=data, warm_size=warm_size, batch_size=batch_size, mem_batches=mem_batches, png_folder="dev/plots_baseline/")
    
    # plotting the results
    grf <- har_plot(result$detector, head(data$serie,nrow(result$detection)), result$detection, head(data$event,nrow(result$detection)))
    exp_res_SCP[[curve]][[var]] <- list(res=result, plot=grf)
    
    res_df <- cbind(curve=curve,var=var,lag=NA,mdl=NA,ic=TRUE,par=mdl_train,
                    recursive=NA,mem=NA,alpha=alpha,lookback=TRUE,
                    detlag_batch=result$detection_lag$batch,
                    detlag_obs=result$detection_lag$obs,
                    detlag_time=result$detection_lag$time,exp="SCP")
    exp_res_SCP_df <- rbind(exp_res_SCP_df,res_df)
    
    folder <- paste0(getwd(),"/dev/plots_test/",name_data,"/")
    if (!dir.exists(folder)) dir.create(folder)
    saveRDS(list(res=exp_res_SCP,df_res=exp_res_SCP_df),paste0(folder,"res_exp_SCP.rds"))
  }
}

