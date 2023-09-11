source("dev/herald_conceptual_tests_functions.R")

library(ggplot2)
library(dplyr)

#install.packages("strucchange")
library(strucchange)



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


# Toy example!
df <- exp_series[["linear"]][["low"]]
serie <- df$serie




# Tests and breakpoints ========================================================
#pure Chow test
sctest(serie ~ 1, type = "Chow", point = 100)
## test the model null hypothesis that the trend remains 
## constant over time for potential break points between
## the observations in 0.15 of data (from = 0.15) and 0.85 of data (to = 0.85)
## compute F statistics
fs <- Fstats(serie ~ 1, from = 0.15, to = 0.85)
## plot the F statistics
plot(fs, alpha = 0.05)
## add the boundary in another colour
lines(boundary(fs, alpha = 0.05), col = 2)
## and the corresponding p values
#plot(fs, pval = TRUE, alpha = 0.05)
## add the boundary in another colour
#lines(boundary(fs, pval = TRUE, alpha = 0.05), col = 2)
## F statistics indicate one breakpoint
bp <- breakpoints(fs)
bp
lines(bp)
plot(ts(serie))
abline(v=bp$breakpoints, lty=2)
## perform the supF test
sctest(fs, type="supF")
## confidence interval
#ci.serie <- confint( breakpoints(serie ~ 1))
#ci.serie
#lines(ci.serie)


# Summary:
## compute F statistics
fs <- Fstats(serie ~ 1, from = 0.15, to = 0.85)
## plot the F statistics
plot(fs, alpha = 0.05)
## F statistics indicate one breakpoint
bp <- breakpoints(fs)
bp
lines(bp)
plot(ts(serie))
abline(v=bp$breakpoints, lty=2)
## perform the supF test
sctest(fs, type="supF")



# Online monitoring ============================================================
#Not sure it uses the Chow test!
## use the first 50 observations as history period
warmup <- 50
batch_size <- 1
## create efp (Empirical Fluctuation Process) object
e1 <- efp(serie~1, data=df[1:warmup,,drop=FALSE], type="ME", h=1)
## monitoring object for efp (Empirical Fluctuation Processes)
me <- mefp(e1, alpha=0.05)
for(i in seq(warmup,nrow(df),by=batch_size)){
  ## monitor the next observation
  me <- monitor(me, data=df[1:i,,drop=FALSE])
  if(!is.na(me$breakpoint)) {
    browser()
    cat(paste0("Caught breakpoint at ",me$breakpoint," in time ",i))
    plot(me)
    #break
  }
}





# Sensibility analysis =========================================================
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
      
      plot(ts(serie[1:i]))
      points(x=100,y=serie[100], pch = 19,col="blue")
      points(x=bp_point,y=serie[bp_point], pch = 20,col="red")
      abline(v=bps, lty=2)
      
      break #comment this line for full result
    }
    break_detected <- FALSE
    bt <- bt + 1
  }
  
  return(list(cp=bps,lag=detection_lag))
}

data <- exp_series[["linear"]][["high"]]
cps <- baseline_cp(data, mem=27, conf=0.05, cp_from=0.15, full=FALSE)
cat(paste0("Detection lag:\n",cps$lag$batch," batches\n",
           cps$lag$obs," observations\n",
           cps$lag$time," seconds"))


