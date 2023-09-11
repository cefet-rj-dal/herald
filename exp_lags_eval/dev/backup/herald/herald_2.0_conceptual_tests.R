source("dev/herald_conceptual_tests_functions.R")

library(ggplot2)
library(dplyr)


# General experimental settings --------------------------------------------------------
parameters <- list(
  #Herald hiperparameters
  mdl_train = c(27,81),#c(27,81,243), #m
  lag_pred = c(1,2,3,9,27), #q = 0,1,2 (anom); 3,9,27 (cp)?
  mem_dist = 27,#c(9,27,81), #k
  ma_order = c(3,9,27), #o
  pred_ic = c(FALSE,TRUE),
  pred_recursive = c(FALSE,TRUE),
  alpha = 1.5,#c(1.5,3), #a
  lookback = FALSE, #c(FALSE,TRUE), #look at past outliers in memory?
  full_mem = TRUE, #c(FALSE,TRUE), #online execution in full memory?
  
  #Herald hiperparameters (fixed)
  online_step = 1,
  dist_func = function(v1, v2) sqrt(sum((v1 - v2)^2)), #euclidean
  seq_dist_path = paste0(getwd(),"/dev/seq_dist.rds"),
  
  #plot hiperparameters (fixed)
  png_folder = paste0(getwd(),"/dev/plots_test/"),
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
tests <- data.frame(lag_pred=c(0,0,NA,NA,NA,NA),
               pred_ic=c(FALSE,TRUE,FALSE,TRUE,FALSE,TRUE),
               pred_recursive=c(FALSE,FALSE,FALSE,FALSE,TRUE,TRUE),
               exp=c("A","B","C","D","E","F"))

run_test <- function(data,data_name,parameters,test,verbose=TRUE){
  par_test <- parameters
  if(!is.na(test$lag_pred)) par_test$lag_pred <- test$lag_pred
  par_test$pred_ic <- test$pred_ic
  par_test$pred_recursive <- test$pred_recursive
  par_test$png_folder <- paste0(getwd(),"/dev/plots_test/",data_name,"/")
  par_test$res_path <- paste0(par_test$png_folder,"res_exp_",test$exp,".rds")
  
  if (!dir.exists(par_test$png_folder)) dir.create(par_test$png_folder)
  
  res <- experiments(data,par_test,verbose=verbose)
  return(list(res))
}

data_tests <- function(data,data_name,parameters,tests,verbose=TRUE){
  results <- sapply(tests$exp, function(t) run_test(data,data_name,parameters,tests[tests$exp==t,],verbose=verbose),USE.NAMES = TRUE)
  coerced_df <- do.call(rbind, lapply(names(results), function(x) {
    data <- results[[x]]$results
    data$exp <- x
    data
  }))
  
  saveRDS(coerced_df,file=paste0(getwd(),"/dev/plots_test/",data_name,"/","res_exp_",data_name,".rds"))
  
  return(list(results=results,df_results=coerced_df))
}

best_results <- function(results,data_name,parameters){
  best_res <- results$df_results %>%
    filter(mem >=27) %>%
    group_by(exp) %>%
    filter(detlag_obs == min(detlag_obs, na.rm = TRUE)) %>%
    slice(1) %>%
    ungroup()
  
  saveRDS(best_res,file=paste0(getwd(),"/dev/plots_test/",data_name,"/","best_res_",data_name,".rds"))
  
  plots <- lapply(1:nrow(best_res), function(r) {
    exp <- best_res[r,]$exp
    id <- best_res[r,]$id
    p <- results$results[[exp]]$detection[[id]]$detection$plot
    png_file <- paste0(parameters$png_folder,data_name,"/bestExp",exp,".png")
    png(file=png_file,width=600, height=350)
    plot(p)
    dev.off()
    return(p)
  })
  
  return(list(df=best_res,plots=plots))
}

#First conceptual tests ========================================================

# Generate the simulated time series
# Specify the trend function ("linear", "log", "sqrt", "poly2", "poly3", "exp", "sigmoid")
n <- 100
seed <- 123
volatility <- c(low=0.025,medium=0.05,high=0.1) #volatility levels
exp_series <- sapply(c("linear","log","sqrt","poly2","poly3","exp","sigmoid"), function(curve) {
  sapply(names(volatility), function(v)
    trend_sym(n, ts.var=volatility[v], seed=seed, trend_function=curve),simplify = FALSE, USE.NAMES = TRUE)
  },simplify = FALSE, USE.NAMES = TRUE
)

exp_res <- sapply(names(exp_series), 
                function(curve){
                  sapply(names(exp_series[[curve]]), 
                         function(vol) {
                            data <- exp_series[[curve]][[vol]]
                            name_data <- paste0(curve,"_",vol)
                            results_data <- data_tests(data,name_data,parameters,tests,verbose=FALSE)
                            best_res_data <- best_results(results_data,name_data,parameters)
                            list(res=results_data,best=best_res_data)
                          },simplify = FALSE, USE.NAMES = TRUE)
                },simplify = FALSE, USE.NAMES = TRUE)

saveRDS(exp_res,file=paste0(getwd(),"/dev/plots_test/exp_res.rds"))


#Gradual trend with little volatility
#data <- trend_sym(n,ts.var=0.001,seed=seed, trend_function="linear")
#results_data1 <- data_tests(data,"data1",parameters,tests)
#best_res_data1 <- best_results(results_data1,"data1",parameters)

#Gradual trend with volatility
#data <- trend_sym(n,ts.var=0.01,seed=seed, trend_function=trend_function)
#results_data2 <- data_tests(data,"data2",parameters,tests)
#best_res_data2 <- best_results(results_data2,"data2",parameters)

#Gradual trend with a lot of volatility
#data <- trend_sym(n,ts.var=0.1,seed=seed, trend_function=trend_function)
#results_data3 <- data_tests(data,"data3",parameters,tests)
#best_res_data3 <- best_results(results_data3,"data3",parameters)



#Analysis of results ========================================================

# Concatenate results from all series
exp_res_df <- data.frame()
for(curve in names(exp_res)){
  for(var in names(exp_res[[curve]])){
    df <- exp_res[[curve]][[var]]$res$df_results
    exp_res_df <- rbind(exp_res_df, cbind(curve=curve,var=var,df))
  }
}

# Create the ggplot with faceted boxplots
ggplot(exp_res_df, aes(x = exp, y = detlag_obs)) +
  geom_boxplot(fill="deepskyblue") +
  stat_summary(fun = "mean", geom = "point", shape = 17, size = 3, color = "red") +
  stat_summary(fun = "mean", geom = "line", aes(group = 1), size = 1, color = "red") +
  facet_wrap(~ factor(var, levels=c('low', 'medium', 'high')), scales = "free_x", ncol = 1) +
  labs(title = "Boxplots of detection lags grouped by data volatility",
       x = "Herald settings",
       y = "Detection lag")+
  theme_bw()


# Concatenate best results from all series
best_res_df <- data.frame()
for(curve in names(exp_res)){
  for(var in names(exp_res[[curve]])){
    df <- exp_res[[curve]][[var]]$best$df
    best_res_df <- rbind(best_res_df, cbind(curve=curve,var=var,df))
  }
}

# Create the ggplot with faceted barplots
best_res_df %>%
  group_by(curve,var) %>%
  arrange(detlag_obs) %>%
  slice_head(n = 3) %>%
  group_by(var, exp) %>%
  summarise(frequency = n()) %>%
  ggplot(aes(x = exp, y = frequency)) +
    geom_bar(stat = "identity",fill="deepskyblue") +
    facet_wrap(~ factor(var, levels=c('low', 'medium', 'high')), ncol = 1) +
    labs(title = "Frequency of settings among the top 3",
         x = "Herald settings",
         y = "Frequency")+
  theme_bw()
