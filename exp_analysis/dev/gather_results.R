# Gathering of results (from folders) ==========================================

## Gathering of results (lags) -------------------------------------------------

curves <- c("linear","log","sqrt","poly2","poly3","exp","sigmoid")
volatility <- c(low=0.025,medium=0.05,high=0.1) #volatility levels
seeds <- as.character(seq(from=123,by=1,length.out=30))

# Concatenate results from all series folders...................................
exp_res_df <- data.frame()
for(curve in curves){
  for(var in names(volatility)){
    for(seed in seeds){
      data_name <- paste0(curve,"_",var,"_",seed)
      df <- readRDS(paste0(getwd(),"/dev/plots/",data_name,"/","res_exp_",data_name,".rds"))
      exp_res_df <- rbind(exp_res_df, cbind(curve=curve,var=var,seed=seed,df))
    }
  }
}
saveRDS(exp_res_df,file=paste0(getwd(),"/dev/plots/exp_res_df.rds"))


## Gathering of plots  (lags) --------------------------------------------------

df_baseline <- exp_res_df_mem  %>%
  filter(exp %in% c("baseline") & var == "high" & seed == 123 & mem==27) %>%
  distinct(curve,var,seed,mem,lag,mdl,par, .keep_all= TRUE)
df_herald <- exp_res_df_mem  %>%
  filter(exp %in% c("D") & var == "high" & seed == 123 & mem==27 & lag==27 & mdl==27 & par==27) %>%
  distinct(curve,var,seed,mem,lag,mdl,par, .keep_all= TRUE)
df_res <- rbind(df_baseline,df_herald)

exp_res_plots <- list()
for(crv in unique(df_res$curve)){
  exp_res_plots[[crv]] <- list()
  for(e in unique(df_res$exp)){
    b <- df_res %>% filter(curve == crv, exp == e)
    data_name <- paste0(b$curve,"_",b$var,"_",b$seed)
    det <- readRDS(paste0(getwd(),"/dev/plots/plots/",data_name,"/","det_res_exp_",data_name,".rds"))
    exp_res_plots[[crv]][[e]] <- det[[as.character(b$exp)]]$detection[[b$id]]$detection$plot
  }
}

df_herald <- exp_res_df  %>%
  filter(exp %in% c("A") & var == "high" & seed == 123 & mem==27 & mdl==27) %>%
  distinct(curve,var,seed,mem,lag,mdl,par, .keep_all= TRUE)
df_res <- rbind(df_herald)

for(crv in unique(df_res$curve)){
  for(e in unique(df_res$exp)){
    b <- df_res %>% filter(curve == crv, exp == e)
    data_name <- paste0(b$curve,"_",b$var,"_",b$seed)
    det <- readRDS(paste0(getwd(),"/dev/plots/plots/plots/",data_name,"/","det_res_exp_",data_name,".rds"))
    exp_res_plots[[crv]][[e]] <- det[[as.character(b$exp)]]$detection[[b$id]]$detection$plot
  }
}
saveRDS(exp_res_plots,file=paste0(getwd(),"/dev/plots/exp_res_plots.rds"))


## Gathering of results (lags + eval) ------------------------------------------

curves <- c("linear","log","sqrt","poly2","poly3","exp","sigmoid")
volatility <- c(low=0.025,medium=0.05,high=0.1) #volatility levels
seeds <- as.character(seq(from=123,by=1,length.out=30))

# Concatenate results from all series folders...................................
exp_res_df_lags <- data.frame()
for(curve in curves){
  for(var in names(volatility)){
    for(seed in seeds){
      data_name <- paste0(curve,"_",var,"_",seed)
      df <- readRDS(paste0(getwd(),"/dev/plots/",data_name,"/","res_exp_",data_name,"_lags.rds"))
      exp_res_df_lags <- rbind(exp_res_df_lags, cbind(curve=curve,var=var,seed=seed,df))
    }
  }
}
exp_res_df_lags[c(4,6,7,9,11:16)] <- sapply(exp_res_df_lags[c(4,6,7,9,11:16)],as.numeric)
exp_res_df_lags[c(5,8,10)] <- sapply(exp_res_df_lags[c(5,8,10)],function(v) as.numeric(as.logical(v)))
saveRDS(exp_res_df_lags,file=paste0(getwd(),"/dev/plots/exp_res_df_lags.rds"))

exp_res_df_eval <- data.frame()
for(curve in curves){
  for(var in names(volatility)){
    for(seed in seeds){
      data_name <- paste0(curve,"_",var,"_",seed)
      df <- readRDS(paste0(getwd(),"/dev/plots/",data_name,"/","res_exp_",data_name,"_eval.rds"))
      exp_res_df_eval <- rbind(exp_res_df_eval, cbind(curve=curve,var=var,seed=seed,df))
    }
  }
}
saveRDS(exp_res_df_eval,file=paste0(getwd(),"/dev/plots/exp_res_df_eval.rds"))


## Gathering of plots (lags + eval) --------------------------------------------

df_baseline <- exp_res_df_mem_eval  %>%
  filter(exp %in% c("baseline") & var == "high" & seed == 123 & mem==27) %>%
  distinct(curve,var,seed,mem,lag,mdl,par, .keep_all= TRUE)
df_herald <- exp_res_df_mem_eval  %>%
  filter(exp %in% c("D") & var == "high" & seed == 123 & mem==27 & lag==27 & mdl==27 & par==27) %>%
  distinct(curve,var,seed,mem,lag,mdl,par, .keep_all= TRUE)
df_res <- rbind(df_baseline,df_herald)

exp_res_plots_eval <- list()
for(crv in unique(df_res$curve)){
  exp_res_plots_eval[[crv]] <- list()
  for(e in unique(df_res$exp)){
    b <- df_res %>% filter(curve == crv, exp == e)
    data_name <- paste0(b$curve,"_",b$var,"_",b$seed)
    det <- readRDS(paste0(getwd(),"/dev/plots/plots/",data_name,"/","det_res_exp_",data_name,".rds"))
    exp_res_plots_eval[[crv]][[e]] <- det[[as.character(b$exp)]]$detection[[b$id]]$detection$plot$final
  }
}

saveRDS(exp_res_plots_eval,file=paste0(getwd(),"/dev/plots/exp_res_plots_eval.rds"))



# Get consolidated saved results ===============================================

## Gathering of results (lags) -------------------------------------------------

#Contains results for baseline and settings A-F
exp_res_df <- readRDS(file=paste0(getwd(),"/dev/plots/exp_res_df.rds"))

exp_res_df$curve <- factor(exp_res_df$curve,levels=curves)
exp_res_df$var <- factor(exp_res_df$var, ordered = TRUE, levels=names(volatility))
exp_res_df$seed <- factor(exp_res_df$seed, levels=seeds)
exp_res_df$exp <- factor(exp_res_df$exp, ordered = TRUE, levels=c("baseline","A","B","C","D","E","F"))


#Contains results for baseline and D for memory sensibility analysis
exp_res_df_mem <- readRDS(paste0(getwd(),"/dev/plots/exp_res_df_0D_mem.rds"))

exp_res_df_mem$curve <- factor(exp_res_df_mem$curve,levels=curves)
exp_res_df_mem$var <- factor(exp_res_df_mem$var, ordered = TRUE, levels=names(volatility))
exp_res_df_mem$seed <- factor(exp_res_df_mem$seed, levels=seeds)
exp_res_df_mem$exp <- factor(exp_res_df_mem$exp, ordered = TRUE, levels=c("baseline","D"))


#Contains results for baseline and D for memory sensibility analysis
exp_res_df_mem_lags <- readRDS(paste0(getwd(),"/dev/plots/exp_res_df_0D_mem_lags.rds"))

exp_res_df_mem_lags$curve <- factor(exp_res_df_mem_lags$curve,levels=curves)
exp_res_df_mem_lags$var <- factor(exp_res_df_mem_lags$var, ordered = TRUE, levels=names(volatility))
exp_res_df_mem_lags$seed <- factor(exp_res_df_mem_lags$seed, levels=seeds)
exp_res_df_mem_lags$exp <- factor(exp_res_df_mem_lags$exp, ordered = TRUE, levels=c("baseline","D"))

#Contains example plots for baseline and settings A and D
exp_res_plots <- readRDS(paste0(getwd(),"/dev/plots/exp_res_plots.rds"))



## Gathering of results (eval) ------------------------------------------

#Contains results for baseline and D for memory sensibility analysis
exp_res_df_mem_eval <- readRDS(paste0(getwd(),"/dev/plots/exp_res_df_0D_mem_eval.rds"))

exp_res_df_mem_eval$curve <- factor(exp_res_df_mem_eval$curve,levels=curves)
exp_res_df_mem_eval$var <- factor(exp_res_df_mem_eval$var, ordered = TRUE, levels=names(volatility))
exp_res_df_mem_eval$seed <- factor(exp_res_df_mem_eval$seed, levels=seeds)
exp_res_df_mem_eval$exp <- factor(exp_res_df_mem_eval$exp, ordered = TRUE, levels=c("baseline","D"))

#Contains example plots for baseline and D
exp_res_plots_eval <- readRDS(paste0(getwd(),"/dev/plots/exp_res_plots_eval.rds"))