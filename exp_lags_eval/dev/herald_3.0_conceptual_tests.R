source("dev/herald_3.0_conceptual_tests_functions.R")

library(ggplot2)
library(dplyr)


# General experimental settings --------------------------------------------------------
parameters <- list(
  #Herald hiperparameters
  pred = TRUE, #perform predictions
  mdl_train = c(27,81),#c(27,81,243), #m
  lag_pred = c(1,2,3,9,27), #q = 0,1,2 (anom); 3,9,27 (cp)?
  mem_dist = c(9,27,81),#c(9,27,81), #k
  ma_order = c(3,9,27), #o
  pred_ic = c(FALSE,TRUE),
  pred_recursive = c(FALSE,TRUE),
  conf = 0.05,#c(0.01,0.05), #a
  full_mem = TRUE, #c(FALSE,TRUE), #online execution in full memory
  
  #Herald hiperparameters (fixed)
  online_step = 1,
  dist_func = function(v1, v2) sqrt(sum((v1 - v2)^2)), #euclidean
  seq_dist_path = paste0(getwd(),"/dev/seq_dist.rds"),
  
  #plot hiperparameters (fixed)
  png_folder = paste0(getwd(),"/dev/plots/"),
  verbose = FALSE,
  plot_title = TRUE,
  
  #experiment results file
  res_path = paste0(getwd(),"/dev/res_exp/res_exp.rds")
)


#Exp A - No lag, no IC
#Exp B - No lag, MA IC
#Exp C - Lag, no IC
#Exp D - Lag, MA IC
#Exp E - Lag, recursive, no IC
#Exp F - Lag, recursive, MA IC
tests <- data.frame(pred=c(FALSE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE),
                    lag_pred=c(0,0,0,NA,NA,NA,NA),
                    pred_ic=c(FALSE,FALSE,TRUE,FALSE,TRUE,FALSE,TRUE),
                    pred_recursive=c(FALSE,FALSE,FALSE,FALSE,FALSE,TRUE,TRUE),
                    exp=c("baseline","A","B","C","D","E","F"))

run_test <- function(data,data_name,parameters,test,verbose=TRUE){
  par_test <- parameters
  if(test$pred){
    par_test$pred <- test$pred
    if(!is.na(test$lag_pred)) par_test$lag_pred <- test$lag_pred
    par_test$pred_ic <- test$pred_ic
    par_test$pred_recursive <- test$pred_recursive
  }
  else{
    par_test$pred <- test$pred
    par_test$lag_pred <- test$lag_pred
    par_test$mdl_train <- 0
    par_test$ma_order <- 0
    par_test$pred_ic <- test$pred_ic
    par_test$pred_recursive <- test$pred_recursive
    par_test$full_mem <- FALSE
  }
  
  par_test$png_folder <- paste0(getwd(),"/dev/plots/",data_name,"/")
  par_test$res_path <- paste0(par_test$png_folder,"res_exp_",test$exp,".rds")
  
  if (!dir.exists(par_test$png_folder)) dir.create(par_test$png_folder)
  
  res <- experiments(data,par_test,verbose=verbose)
  return(list(res))
}

data_tests <- function(data,data_name,parameters,tests,verbose=TRUE){
  results <- sapply(tests$exp, function(t) run_test(data,data_name,parameters,tests[tests$exp==t,],verbose=verbose),USE.NAMES = TRUE)
  saveRDS(results,file=paste0(getwd(),"/dev/plots/",data_name,"/","det_res_exp_",data_name,".rds"))
  
  coerced_df_lags <- do.call(rbind, lapply(names(results), function(x) {
    data <- results[[x]]$results$lags
    data$exp <- x
    data
  }))
  coerced_df_eval <- do.call(rbind, lapply(names(results), function(x) {
    data <- results[[x]]$results$eval
    data$exp <- x
    data
  }))
  
  saveRDS(coerced_df_lags,file=paste0(getwd(),"/dev/plots/",data_name,"/","res_exp_",data_name,"_lags.rds"))
  saveRDS(coerced_df_eval,file=paste0(getwd(),"/dev/plots/",data_name,"/","res_exp_",data_name,"_eval.rds"))
  
  return(list(results=results,df_results_lags=coerced_df_lags,df_results_eval=coerced_df_eval))
}


#First conceptual tests ========================================================

# Generate the simulated time series
# Specify the trend function ("linear", "log", "sqrt", "poly2", "poly3", "exp", "sigmoid")
n <- 100
curves <- c("linear","log","sqrt","poly2","poly3","exp","sigmoid")
volatility <- c(low=0.025,medium=0.05,high=0.1) #volatility levels
seed <- as.character(seq(from=123,by=1,length.out=30))

exp_series <- sapply(curves, function(curve) {
    sapply(names(volatility), function(v){
      sapply(seed, function(s) {
        trend_sym(n, ts.var=volatility[v], seed=as.numeric(s), trend_function=curve, flip=TRUE)
      },simplify = FALSE, USE.NAMES = TRUE)
    },simplify = FALSE, USE.NAMES = TRUE)
  },simplify = FALSE, USE.NAMES = TRUE)

test_series <- exp_series
exp_tests <- tests[c(1,5),]

for(curve in names(test_series)){
  for(vol in names(test_series[[curve]])){
    for(seed in names(test_series[[curve]][[vol]])){
      data <- exp_series[[curve]][[vol]][[seed]]
      name_data <- paste0(curve,"_",vol,"_",seed)
      cat(paste0("Serie: ",name_data,", started: ",Sys.time(),"\n"))
      data_tests(data,name_data,parameters,exp_tests,verbose=FALSE)
    }
  }
}

#Gradual trend with little volatility
#data <- trend_sym(n,ts.var=0.001,seed=seed, trend_function="linear")
#results_data1 <- data_tests(data,"data1",parameters,tests)

#Gradual trend with volatility
#data <- trend_sym(n,ts.var=0.01,seed=seed, trend_function=trend_function)
#results_data2 <- data_tests(data,"data2",parameters,tests)

#Gradual trend with a lot of volatility
#data <- trend_sym(n,ts.var=0.1,seed=seed, trend_function=trend_function)
#results_data3 <- data_tests(data,"data3",parameters,tests)
