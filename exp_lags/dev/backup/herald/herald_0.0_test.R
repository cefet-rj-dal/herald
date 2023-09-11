# Installing and loading Harbinger -------------------------------------------------------

# install.packages("devtools")
library(devtools)
devtools::install_github("cefet-rj-dal/harbinger", force=TRUE, dependencies=FALSE, upgrade="never", build_vignettes = TRUE)
devtools::install_github("cefet-rj-dal/daltoolbox", force=TRUE, dependencies=FALSE, upgrade="never", build_vignettes = TRUE)
devtools::install_github("cefet-rj-dal/event_datasets", force=TRUE)

#Loading harbinger
library(harbinger)
# Loading Nexus
source("https://raw.githubusercontent.com/cefet-rj-dal/harbinger/master/develop/nexus.R")


# Simulated time series -------------------------------------------------------------------
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

nonstat_ts <- nonstationarity_sym(ts.len=200,ts.mean=0,ts.var=0.1,ts.trend=0.01)
data <- data.frame(serie=nonstat_ts,event=FALSE)
data[200, "event"] <- TRUE

require(ggplot2)
#Plotting time series
ggplot(data, aes(y=serie, x=1:nrow(data))) +
  geom_line()+
  xlab("Time")+
  ylab("")+#ylab(serie_name)+
  theme_bw()



# Installing and loading Harbinger -------------------------------------------------------

# install.packages("TSPred")
require("TSPred")


# Function definitions  -------------------------------------------------------------------

#inertial component decomposition
decomp_fun <- function(serie,kernel=c("box", "normal"),bandwidth=0.5){
  stats::ksmooth(time(serie),serie,n.points=length(serie),kernel=kernel,bandwidth=bandwidth)$y
}
#plot(data$serie,type='l')
#lines(ksmooth(time(data$serie),data$serie,'normal',bandwidth=10),type='l',col='blue')
decomp_fun <- function(serie) serie


#prediction of inertial component
pred_fun <- function(serie, n.ahead, window_len, size, onestep=FALSE){
  require("TSPred")

  proc_subset <- subsetting()
  proc_mas <- MAS()
  proc_sw <- SW( window_len = window_len )
  proc_mm <- MinMax()
  
  #Obtaining objects of the modeling class
  #modl <- NNET( size = size, sw = proc_sw, proc = list(MM = proc_mm) )
  modl <- ARIMA()

  tspred_obj <- tspred( subsetting = proc_subset, 
                        processing = list(MAS = proc_mas), 
                        modeling = modl,
                        n.ahead=n.ahead)
  
  tspred_res <- tspred_obj %>% 
    subset(data = serie) %>%
    preprocess(prep_test = TRUE) %>% 
    train() %>%
    predict(input_test_data = TRUE, onestep=onestep) %>% 
    postprocess()
  
  return(tspred_res$pred$raw)
}

#detection over predicted inertial component
detect_fun <- function(serie,pred_ic,alpha = 1.5){
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
  s <- (serie-pred_ic)^2
  outliers <- outliers.boxplot.index(s, alpha)
  group_outliers <- split(outliers, cumsum(c(1, diff(outliers) != 1)))
  outliers <- rep(FALSE, nrow(s))
  for (g in group_outliers) {
    if (length(g) > 0) {
      i <- min(g)
      outliers[i] <- TRUE
    }
  }
  i_outliers <- rep(NA, n)
  i_outliers[non_na] <- outliers
  
  detection <- data.frame(idx=1:n, event = i_outliers, type="")
  detection$type[i_outliers] <- "change point"
  
  return(detection)
}


# Parameter definitions  -------------------------------------------------------------------

#decomp_par <- list(kernel="normal",bandwidth=0.1)
pred_par <- list(window_len=15,size = 5,onestep=TRUE)
detect_par <- list(alpha = 1.5)

lag_pred <- 1
online_step <- 1

# establishing method
model <- har_herald(lag_pred=lag_pred, online_step=online_step, 
                    decomp_fun, list(),
                    pred_fun, pred_par,
                    detect_fun, detect_par)


# Herald detection test  -------------------------------------------------------------------

# making detections using method
detection <- model |>
  detect(data$serie)

# filtering detected events
print(detection |> dplyr::filter(event==TRUE))

# evaluating the detections
evaluation <- harbinger::evaluate(model, detection$event, data$event)
print(evaluation$confMatrix)

# ploting the results
grf <- plot.harbinger(model, data$serie, detection, data$event)
plot(grf)


# Herald Nexus detection -------------------------------------------------------------------

run_nexus <- function(model, data, warm_size = 30, batch_size = 15, mem_batches = 0, png_folder="dev/plots/"){
  
  datasource <- nex_simulated_datasource("data", data$serie)
  
  online_detector <- nexus(datasource, model, warm_size = warm_size, batch_size = batch_size, mem_batches = mem_batches)
  
  online_detector <- warmup(online_detector)
  
  bt <- 1
  event_happened <- FALSE
  bt_event_happened <- 0
  event_idx <- 0
  
  while (!is.null(online_detector$datasource)) {
    
    online_detector <- detect(online_detector)
    #print(table(online_detector$detection$event))
    
    #for plotting
    if(any(data$event[online_detector$detection$idx]) & !event_happened){
      event_happened <- TRUE
      event_idx <- which(data$event[online_detector$detection$idx])
      bt_event_happened <- bt
    }
    
    png(file=paste0(png_folder,"batch_",bt,".png"),
        width=600, height=350)
    grf <- plot.harbinger(online_detector$detector, data$serie[online_detector$detection$idx], online_detector$detection, data$event[online_detector$detection$idx])
    plot(grf)
    dev.off()
    
    bt <- bt + 1
    
    if(any(which(as.logical(online_detector$detection$event))>event_idx) & event_happened) {
      break
    }#comment this line for full result 
    #end
  }
  
  return(list(detector=online_detector$detector,detection=online_detector$detection))
}


# establishing method
model <- har_herald(lag_pred=lag_pred, online_step=online_step, 
                    decomp_fun, decomp_par,
                    pred_fun, pred_par,
                    detect_fun, detect_par)

#FULL MEMORY
herald_full <- run_nexus(model,data[100:250,], warm_size = 30, batch_size = 15, mem_batches = 0, png_folder="rebecca_dal/dev/plots_full/")

#library(gifski)
#invisible(save_gif(herald_nexus(model,data, warm_size = 30, batch_size = 15, mem_batches = 0),
#                   "herald_full.gif", delay = 0.05, width = 800, 
#                   height = 300, progress = FALSE))
#include_graphics("herald_full.gif")


#BATCH MEMORY
herald_batch <- run_nexus(model,data[100:250,], warm_size = 30, batch_size = 15, mem_batches = 3, png_folder="rebecca_dal/dev/plots_batch/")

#library(gifski)
#invisible(save_gif(herald_nexus(model,data, warm_size = 30, batch_size = 15, mem_batches = 3),
#                   "herald_batch.gif", delay = 0.05, width = 800, 
#                   height = 300, progress = FALSE))
#include_graphics("herald_batch.gif")


# EDfP Nexus detection -------------------------------------------------------------------

# establishing method
model <- har_edfp(online_step=2,
                  pred_fun, pred_par=pred_par,
                  detect_fun, detect_par)


edfp_full <- run_nexus(model,data[100:250,], warm_size = 30, batch_size = 15, mem_batches = 3, png_folder="rebecca_dal/dev/plots_edpf_full/")
