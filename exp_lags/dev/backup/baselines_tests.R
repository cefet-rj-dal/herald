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
ma_order <- 9 #c(3,9,27), #o
alpha <- 3 #c(1.5,3), #a
mem_batches <- mdl_train + lag_pred + 1
warm_size <- mem_batches
batch_size <- online_step



## --------------------------------------------------------
# Run Experiments

# Generate the simulated time series
# Specify the trend function ("linear","log","poly2","poly3","exp","sigmoid")

n <- 100
trend_function <- "sigmoid"
seed <- 123
ts.var <- c(0.001,0.01,0.1)

#Gradual trend with little volatility
data <- trend_sym(n,ts.var=0.001,seed=seed, trend_function=trend_function)

#Gradual trend with volatility
data <- trend_sym(n,ts.var=0.01,seed=seed, trend_function=trend_function)

#Gradual trend with a lot of volatility
data <- trend_sym(n,ts.var=0.1,seed=seed, trend_function=trend_function)


# establishing method
model <- hcp_cf_arima(w = ma_order, alpha = alpha) #CF using Arima
model <- hcp_scp(sw = mdl_train, alpha = alpha) #SCP


result <- run_nexus(model=model, data=data, warm_size=warm_size, batch_size=batch_size, mem_batches=mem_batches, png_folder="dev/plots_baseline/")
# plotting the results
grf <- har_plot(result$detector, head(data$serie,nrow(result$detection)), result$detection, head(data$event,nrow(result$detection)))
plot(grf)
cat(paste0("Detection lag:\n",result$detection_lag$batch," batches\n",
           result$detection_lag$obs," observations\n",
           result$detection_lag$time," seconds"))


