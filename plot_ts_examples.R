source("dev/herald_3.0_conceptual_tests_functions.R")
library(ggpubr)

n <- 100
curves <- c("linear","log","sqrt","poly2","poly3","exp","sigmoid")
volatility <- c(low=0.025,medium=0.05,high=0.1) #volatility levels
seed <- as.character(seq(from=123,by=1,length.out=30))

exp_series <- sapply(curves, function(curve) {
  sapply(names(volatility), function(v){
    sapply(seed, function(s) {
      trend_sym(n, ts.var=volatility[v], seed=as.numeric(s), trend_function=curve)
    },simplify = FALSE, USE.NAMES = TRUE)
  },simplify = FALSE, USE.NAMES = TRUE)
},simplify = FALSE, USE.NAMES = TRUE)

plots <- list()
for(curve in c("linear","sigmoid")){
  for(vol in c("low","medium","high")){
    data <- exp_series[[curve]][[vol]][[1]]
    plots[[paste0(curve,"_",vol)]] <- ggplot(data, aes(y=serie, x=1:nrow(data))) +
      geom_line()+xlab("Time")+ ylab("")+theme_bw()
  }
}

ggarrange(plotlist=plots, ncol = 3, nrow = 2)
ggsave("CP1_examples.pdf", width = 20, height = 13, units = "cm")


exp_series <- sapply(curves, function(curve) {
  sapply(names(volatility), function(v){
    sapply(seed, function(s) {
      trend_sym(n, ts.var=volatility[v], seed=as.numeric(s), trend_function=curve, flip=TRUE)
    },simplify = FALSE, USE.NAMES = TRUE)
  },simplify = FALSE, USE.NAMES = TRUE)
},simplify = FALSE, USE.NAMES = TRUE)

plots <- list()
for(curve in c("linear","sigmoid")){
  for(vol in c("low","medium","high")){
    data <- exp_series[[curve]][[vol]][[1]]
    plots[[paste0(curve,"_",vol)]] <- ggplot(data, aes(y=serie, x=1:nrow(data))) +
      geom_line()+xlab("Time")+ ylab("")+theme_bw()
  }
}

ggarrange(plotlist=plots, ncol = 3, nrow = 2)
ggsave("CP2_examples.pdf", width = 20, height = 13, units = "cm")