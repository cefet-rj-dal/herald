library(effectsize)

df <- exp_res_df %>%
  filter(exp %in% c("baseline"))
par_comb_baseline <- expand_grid(lag=unique(df$lag),mdl=unique(df$mdl),par=unique(df$par))
par_comb_exp <- list()
for(experiment in c("A","B","C","D","E","F")){
  df <- exp_res_df %>%
    filter(exp %in% c(experiment))
  par_comb_exp[[experiment]] <- expand_grid(lag=unique(df$lag),mdl=unique(df$mdl),par=unique(df$par))
}


df_baseline <- exp_res_df  %>%
  filter(exp %in% c("baseline"))

exp_par_stats <- list()
for(experiment in c("A","B","C","D","E","F")){
  par_comb <- par_comb_exp[[experiment]]
  exp_par_stats[[experiment]] <- list()
  for(comb in 1:nrow(par_comb)){
    if(is.na(par_comb[comb,]$par)){
      df_herald <- exp_res_df  %>%
        filter(exp %in% c(experiment) & 
               lag==par_comb[comb,]$lag & mdl==par_comb[comb,]$mdl & is.na(par_comb[comb,]$par)) %>%
        distinct(curve,var,seed,lag,mdl,par, .keep_all= TRUE)
    }
    else{
      df_herald <- exp_res_df  %>%
        filter(exp %in% c(experiment) & 
               lag==par_comb[comb,]$lag & mdl==par_comb[comb,]$mdl & par==par_comb[comb,]$par)%>%
        distinct(curve,var,seed,lag,mdl,par, .keep_all= TRUE)
    }
   
    df_res <- rbind(df_baseline,df_herald)  
    #browser()
    stats_signif <-  df_res %>%
      group_by(curve,var,seed) %>%
      filter(n() > 1 & !any(is.na(detlag_batch))) %>%
      group_by(var) %>%
      wilcox_test(detlag_batch ~ exp, paired = TRUE, alternative = "greater")  %>%
      adjust_pvalue(method = "bonferroni") %>%
      add_significance("p.adj")
    
    stats_effsize <- df_res %>%
      group_by(curve,var,seed) %>%
      filter(n() > 1 & !any(is.na(detlag_batch))) %>%
      group_by(var) %>%
      group_map(~data.frame(var=unique(.x$var), 
                            effsize=interpret(rank_biserial(.x$detlag_batch, factor(.x$exp,levels=c("baseline",experiment)), paired = TRUE, alternative = "greater"), rules = "funder2019")
      ),.keep=TRUE) %>%
      bind_rows()
    
    exp_par_stats[[experiment]][[as.character(comb)]] <- merge(stats_signif,stats_effsize, sort = F)
  }
}

best_par <- list()
for(experiment in c("A","B","C","D","E","F")){
  par_comb <- par_comb_exp[[experiment]]
  par_stats <- exp_par_stats[[experiment]]
  min_pvalue <- min(sapply(par_stats, function(comb) comb[["p.adj"]]))
  if(min_pvalue>0.05) best_par[[experiment]] <- NULL
  else{
    best_par[[experiment]] <- list()
    #best_par[[experiment]][["min_pvalue"]] <- min_pvalue
    best_par[[experiment]][["stats"]] <- par_stats[sapply(par_stats, function(comb) any(comb[["p.adj"]] < 0.05))]
    best_par[[experiment]][["par"]] <- par_comb[c(as.numeric(names(best_par[[experiment]][["stats"]]))),]
  }
}



best_res_D <- exp_res_df %>%
  filter(exp %in% c("D")) %>%
  group_by(curve,var,seed) %>%
  filter(detlag_batch == min(detlag_batch, na.rm = TRUE)) %>%
  slice(1) %>%
  ungroup() #%>%
  #group_by(lag,mdl,par) %>%
  #summarize(n = n())
max_lag <- names(which.max(table(best_res_D$lag)))
max_mdl <- names(which.max(table(best_res_D$mdl)))
max_par <- names(which.max(table(best_res_D$par)))

par_comb <- par_comb_exp[["D"]]
comb <- which(par_comb$lag==max_lag & par_comb$mdl==max_mdl & par_comb$par==max_par)
par_comb_exp[["D"]][comb,]
exp_par_stats[["D"]][[as.character(27)]]




#SoftED check
df_baseline <- exp_res_df_eval  %>%
  filter(exp %in% c("baseline") & mem ==27)
df_herald <- exp_res_df_eval  %>%
  filter(exp %in% c("D") & mem ==27 & lag==27 & mdl==27 & par==27)
df_res <- rbind(df_baseline,df_herald)
df_res$curve <- factor(df_res$curve,levels=curves)
df_res$var <- factor(df_res$var, ordered = TRUE,levels=names(volatility))
df_res$exp <- factor(df_res$exp,levels=c("baseline","D"))

cols_soft <- names(df_res[grep("_soft", names(df_res))])
stats_soft <- list()
for(col in cols_soft){
   stats <- df_res  %>%
    group_by(curve,var,seed) %>%
    filter(n() > 1) %>%
    group_by(var) %>%
    wilcox_test(as.formula(paste0(col," ~ exp")), paired = TRUE, alternative = "less")  %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")
  
  if(any(stats[["p.adj"]]<0.05))
    stats_soft[[col]] <- stats
}
names(stats_soft)
stats_soft