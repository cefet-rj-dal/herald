library(ggplot2)
library(ggpubr)
library(tidyverse)
library(rstatix)
#install.packages("effectsize") # for rank_biserial
library(effectsize)
#https://easystats.github.io/effectsize/reference/interpret_r.html
#https://easystats.github.io/effectsize/articles/interpret.html
#https://easystats.github.io/effectsize/reference/rank_biserial.html



# Gathering of results =========================================================

curves <- c("linear","log","sqrt","poly2","poly3","exp","sigmoid")
volatility <- c(low=0.025,medium=0.05,high=0.1) #volatility levels
seeds <- as.character(seq(from=123,by=1,length.out=30))

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



# Baseline versus Herald =======================================================

## Statistical tests -----------------------------------------------------------

df_baseline <- exp_res_df  %>%
  filter(exp %in% c("baseline") & mem==27)
df_herald <- exp_res_df  %>%
  filter(exp %in% c("D") & mem==27 & lag==27 & mdl==27 & par==27) %>%
  distinct(curve,var,seed,mem,lag,mdl,par, .keep_all= TRUE)

df_res <- rbind(df_baseline,df_herald)
df_res$exp <- factor(df_res$exp, levels=c("baseline","D"))

stats_signif <- df_res %>%
  wilcox_test(detlag_batch ~ exp, paired = TRUE, alternative = "greater")  %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")%>%
  add_xy_position(x = "exp")
stats_effsize <- df_res %>%
  group_map(~data.frame(effsize=interpret(rank_biserial(.x$detlag_batch, factor(.x$exp,levels=c("baseline","D")), paired = TRUE, alternative = "greater"), rules = "funder2019")
  ),.keep=TRUE) %>%
  bind_rows()
herald_test_total <- merge(stats_signif,stats_effsize, sort = F)
#View(herald_test_total)

stats_signif <-  df_res %>%
  group_by(curve,var,seed) %>%
  filter(n() > 1 & !any(is.na(detlag_batch))) %>%
  group_by(var) %>%
  wilcox_test(detlag_batch ~ exp, paired = TRUE, alternative = "greater")  %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")%>% 
  add_xy_position(x = "exp")

stats_effsize <- df_res %>%
  group_by(curve,var,seed) %>%
  filter(n() > 1 & !any(is.na(detlag_batch))) %>%
  group_by(var) %>%
  group_map(~data.frame(var=unique(.x$var), 
                        effsize=interpret(rank_biserial(.x$detlag_batch, factor(.x$exp,levels=c("baseline","D")), paired = TRUE, alternative = "greater"), rules = "funder2019")
  ),.keep=TRUE) %>%
  bind_rows()
herald_test <- merge(stats_signif,stats_effsize, sort = F)
#View(herald_test)


# Plot .........................................................................
# Box plot
bp <- df_res %>%
  ggplot(aes(x = factor(exp,levels=c("baseline","D")),y=detlag_batch)) +
  geom_boxplot(aes(fill=exp), show.legend = FALSE) +
  scale_fill_manual(values=c("white","#007FFF")) +
  labs(y = "Detection lags")+
  scale_x_discrete(name="Method",labels = c('EDfR','Herald'))+
  theme_bw()+ 
  stat_pvalue_manual(herald_test_total, label = "Wilcoxon, p = {p.adj} {p.adj.signif}, Effect = {round(effsize.r_rank_biserial,2)} ({effsize.Interpretation})", tip.length = 0.01)+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))

bp_var <- df_res %>%
  ggplot(aes(x = factor(exp,levels=c("baseline","D")),y=detlag_batch)) +
  geom_boxplot(aes(fill=exp), show.legend = FALSE) +
  facet_wrap(~factor(var, levels=c('high', 'medium', 'low'), labels=c('High volatility', 'Medium volatility', 'Low volatility')), scales = "free_x",ncol=1) +
  scale_fill_manual(values=c("white","#007FFF")) +
  labs(y = "Detection lags")+
  scale_x_discrete(name="Method",labels = c('EDfR','Herald'))+
  theme_bw()+ 
  stat_pvalue_manual(herald_test, label = "Wilcoxon, p = {p.adj} {p.adj.signif}, Effect = {round(effsize.r_rank_biserial,2)} ({effsize.Interpretation})", tip.length = 0.01)+
  scale_y_continuous(expand = expansion(mult = c(0.3, 0.3)))

ggarrange(bp, bp_var,common.legend = TRUE, legend = "bottom")
ggsave("herald_stats.pdf", width = 20, height = 15, units = "cm")



## Lag (Obs) comparison tests --------------------------------------------------

df_baseline <- exp_res_df  %>%
  filter(exp %in% c("baseline") & mem==27)
df_herald <- exp_res_df  %>%
  filter(exp %in% c("D") & mem==27 & lag==27 & mdl==27 & par==27) %>%
  distinct(curve,var,seed,mem,lag,mdl,par, .keep_all= TRUE)

df_res <- rbind(df_baseline,df_herald)
df_res <- df_res %>% filter(var == "high")
df_res$exp <- factor(df_res$exp, levels=c("baseline","D"))

stats_signif <- df_res %>%
  wilcox_test(detlag_obs ~ exp, paired = TRUE, alternative = "greater")  %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")%>% 
  add_xy_position(x = "exp")
stats_effsize <- df_res %>%
  group_map(~data.frame(effsize=interpret(rank_biserial(.x$detlag_obs, factor(.x$exp,levels=c("baseline","D")), paired = TRUE, alternative = "greater"), rules = "funder2019")
  ),.keep=TRUE) %>%
  bind_rows()
herald_test <- merge(stats_signif,stats_effsize, sort = F)
#View(herald_test)


# Plot .........................................................................
# Box plot
bp_obs <- df_res %>%
  ggplot(aes(x = factor(exp,levels=c("baseline","D")),y=detlag_obs)) +
  geom_boxplot(aes(fill=exp), show.legend = FALSE) +
  scale_fill_manual(values=c("white","#007FFF")) +
  labs(y = "Detection lags (observations)")+
  scale_x_discrete(name="Method",labels = c('EDfR','Herald'))+
  theme_bw()+ 
  stat_pvalue_manual(herald_test, label = "Wilcoxon, p = {p.adj} {p.adj.signif}, Effect = {round(effsize.r_rank_biserial,2)} ({effsize.Interpretation})", tip.length = 0.01)+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2)))

#bp_obs



## Time comparison tests -------------------------------------------------------

df_baseline <- exp_res_df  %>%
  filter(exp %in% c("baseline") & mem==27)
df_herald <- exp_res_df  %>%
  filter(exp %in% c("D") & mem==27 & lag==27 & mdl==27 & par==27) %>%
  distinct(curve,var,seed,mem,lag,mdl,par, .keep_all= TRUE)

df_res <- rbind(df_baseline,df_herald)
df_res <- df_res %>% filter(var == "high")
df_res$exp <- factor(df_res$exp, levels=c("baseline","D"))

stats_signif <- df_res %>%
  wilcox_test(detlag_time ~ exp, paired = TRUE, alternative = "less")  %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")%>% 
  add_xy_position(x = "exp")
stats_effsize <- df_res %>%
  group_map(~data.frame(effsize=interpret(rank_biserial(.x$detlag_time, factor(.x$exp,levels=c("baseline","D")), paired = TRUE, alternative = "less"), rules = "funder2019")
  ),.keep=TRUE) %>%
  bind_rows()
herald_test <- merge(stats_signif,stats_effsize, sort = F)
#View(herald_test)


# Plot .........................................................................
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

df_res <- df_res %>% data_summary(varname="detlag_time", groupnames=c("exp"))
df_res$exp <- as.factor(df_res$exp)
herald_test$y.position <- max(df_res$sd)*1.4

# Box plot
bp_time <- df_res %>%
  ggplot(aes(x = factor(exp,levels=c("baseline","D")),y=detlag_time)) +
  geom_bar(aes(fill=exp),stat="identity", color="black", position=position_dodge(), show.legend = FALSE) +
  geom_errorbar(aes(ymin=detlag_time, ymax=detlag_time+sd), width=.2, position=position_dodge(.9)) +
  #geom_boxplot(aes(fill=exp), show.legend = FALSE) +
  scale_fill_manual(values=c("white","#007FFF")) +
  labs(y = "Detection lags (time (sec))")+
  scale_x_discrete(name="Method",labels = c('EDfR','Herald'))+
  theme_bw()+ 
  stat_pvalue_manual(herald_test, label = "Wilcoxon, p = {p.adj} {p.adj.signif}, Effect = {round(effsize.r_rank_biserial,2)} ({effsize.Interpretation})", tip.length = 0.02)+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2)))

#bp_time

ggarrange(bp_obs,bp_time,common.legend = TRUE, legend = "bottom")
ggsave("herald_obs_time.pdf", width = 20, height = 10, units = "cm")


## Example plots ---------------------------------------------------------------

list_plots <- lapply(exp_res_plots, function(curve) curve[c("A","baseline","D")])
list_plots <- unlist(list_plots,recursive=FALSE)
lags <- sapply(list_plots, function(p) str_match(p$labels$title, "\nDetection lag: \\s*(.*?)\\s* batches,") )[2,]

labels <- names(list_plots)
labels <- gsub("^.*\\.", "", labels)
labels <- gsub("D", "Herald", labels)
labels <- gsub("A", "EDfP", labels)
labels <- gsub("baseline", "EDfR", labels)

labels <- paste(labels," (det. lag: ",lags,")",sep="")

list_plots <- lapply(1:length(list_plots),function(p) list_plots[[p]] + ggtitle(labels[p]))

#ggarrange(plotlist=list_plots, nrow=7,ncol=3)

library(gridExtra)
p <- grid.arrange(grobs=list_plots[1:3], layout_matrix = matrix(c(1, 3, 2, 3), nrow = 2))
ggsave("plots_lags.pdf", p, width = 20, height = 10, units = "cm")


# Baseline versus Ablation =====================================================

## Statistical tests -----------------------------------------------------------

df_baseline <- exp_res_df  %>%
  filter(exp %in% c("baseline") & mem==27)

ablation_tests <- list()
for(experiment in c("A","B","C","D")){
  
  df_herald <- exp_res_df  %>%
    filter(exp %in% c(experiment) & mem==27 & (lag==27 | lag==0) & mdl==27 & (par==27 | is.na(par)) )%>%
    distinct(curve,var,seed,mem,lag,mdl,par, .keep_all= TRUE)
  
  df_res <- rbind(df_baseline,df_herald)  
  
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
  
  ablation_tests[[experiment]] <- merge(stats_signif,stats_effsize, sort = F)
}

ablation_tests_df <- do.call(rbind,ablation_tests)
#View(ablation_tests_df)

# Getting data.frame of lag differences ........................................

df_diff <- data.frame()
df_baseline <- exp_res_df  %>%
  filter(exp %in% c("baseline") & mem==27)
for(experiment in c("A","B","C","D")){
  df_herald <- exp_res_df  %>%
    filter(exp %in% c(experiment) & mem==27 & (lag==27 | lag==0) & mdl==27 & (par==27 | is.na(par)) )%>%
    distinct(curve,var,seed,mem,lag,mdl,par, .keep_all= TRUE)
  
  df_res <- inner_join(df_baseline,df_herald, by = c("curve","var","seed","mem")) %>% 
    mutate(diff_detlag_batch = detlag_batch.x - detlag_batch.y,
           diff_detlag_obs = detlag_obs.x - detlag_obs.y,
           diff_detlag_time = detlag_time.x - detlag_time.y)
  
  df_diff <- rbind(df_diff,df_res)
}
df_diff <- inner_join(df_diff,ablation_tests_df, by = join_by(var, exp.y == group2))
df_diff <- df_diff %>%
  filter(var == "high")
#View(df_diff)

#Confirming stats
#df_diff %>% filter(exp.y %in% c("D")) %>%
#  group_by(var) %>% 
#  wilcox_test(diff_detlag_batch~1, alternative = "greater") %>%
#  adjust_pvalue(method = "bonferroni") %>%
#  add_significance("p.adj")


# Plot .........................................................................
# create the discrete variables using cut
df_diff$effsize_cut <- cut(df_diff$effsize.r_rank_biserial, 
                           labels= c("Very large (neg)", "Large (neg)", "Medium (neg)", "Small (neg)", "Very small (neg)", "Tiny (neg)", "Tiny (pos)", "Very small (pos)", "Small (pos)", "Medium (pos)", "Large (pos)", "Very large (pos)"),
                           breaks = c(-1, -0.4, -0.3, -0.2, -0.1, -0.05, 0, 0.05, 0.1, 0.2, 0.3, 0.4, 1),
                           include.lowest = TRUE, right = FALSE)
levels_count <- length(levels(df_diff$effsize_cut))

# Create the colors scale coresponded to the levels in cut
color_scales_fn <- colorRampPalette(c("red", "#007FFF"))
manual_color <- color_scales_fn(levels_count)
names(manual_color) <- levels(df_diff$effsize_cut)
# Create the sizes scale 
manual_size <-  c(3, 2.5, 2, 1.5, 1, 0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3)
names(manual_size) <- levels(df_diff$effsize_cut)

bp <- df_diff %>%
  ggplot(aes(x = factor(exp.y,levels=c("D","C","B","A")))) +
  geom_boxplot(aes(y=diff_detlag_batch,fill=effsize_cut)) +
  geom_hline(linetype = "dashed",yintercept=0, size = 0.5, color = "gray28") +
  #facet_wrap(~factor(var, levels=c('high', 'medium', 'low'), labels=c('High volatility', 'Medium volatility', 'Low volatility')), scales = "free_x",ncol=1) +
  scale_fill_manual(values=manual_color,name="Effect size",limits = force) +
  labs(y = "Detection lag differences")+
  scale_x_discrete(name = "Ablation settings",labels = c('Herald','- IC','- LP','- IC, - LP (EDfP)'))+
  theme_bw()+  
  geom_text(aes(y=max(diff_detlag_batch,na.rm=TRUE),label = p.adj.signif), position = position_dodge(width=0.9),vjust=-0.7)+
  scale_y_continuous(expand = expansion(mult = c(0.2, 0.3)))

sp <- df_diff %>%
  ggplot(aes(x = factor(exp.y,levels=c("D","C","B","A")))) +
  geom_hline(linetype = "dashed",yintercept=0, size = 0.5, color = "gray28") +
  geom_path(aes(y=effsize.r_rank_biserial, group = 1),color="grey",size=1) +
  geom_point(aes(y=effsize.r_rank_biserial, fill=effsize_cut, color=effsize_cut, size=effsize_cut),shape=21) +
  #facet_wrap(~factor(var, levels=c('high', 'medium', 'low'), labels=c('High volatility', 'Medium volatility', 'Low volatility')), scales = "free_x",ncol=1) +
  scale_size_manual(values = manual_size,name="Effect size",limits = force) +
  scale_color_manual(values = manual_color,name="Effect size",limits = force) +
  scale_fill_manual(values=manual_color,name="Effect size",limits = force) +
  scale_y_continuous(position = "right", name = "Effect size",expand = expansion(mult = c(0.2, 0.3)))+ 
  scale_x_discrete(name = "Ablation settings",labels = c('Herald','- IC','- LP','- IC, - LP (EDfP)'))+
  theme_bw()+  
  geom_text(aes(y=max(effsize.r_rank_biserial),label = p.adj.signif), position = position_dodge(width=0.9),vjust=-0.7)

ggarrange(bp, sp,common.legend = FALSE, legend = "bottom")
ggsave("herald_ablation.pdf", width = 20, height = 8, units = "cm")



# Sensibility analysis =========================================================

## Parameters ------------------------------------------------------------------

df_baseline <- exp_res_df  %>%
  filter(exp %in% c("baseline") & var == "high" & mem==27)

df_herald <- exp_res_df  %>%
  filter(exp %in% c("D") & var == "high" & mem==27) %>%
  distinct(curve,var,seed,mem,lag,mdl,par, .keep_all= TRUE)
par_comb <- list(lag=unique(df_herald$lag),mdl=unique(df_herald$mdl),par=unique(df_herald$par))

exp_par_stats <- list()

for(name_par in names(par_comb)){
  exp_par_stats[[name_par]] <- list()
  for(p in par_comb[[name_par]]){
    l <- ifelse(name_par=="lag",p,27)
    m <- ifelse(name_par=="mdl",p,27)
    o <- ifelse(name_par=="par",p,27)
    
    df_herald <- exp_res_df  %>%
      filter(exp %in% c("D") & var == "high" & mem==27 & lag==l & mdl==m & par==o) %>%
      distinct(curve,var,seed,mem,lag,mdl,par, .keep_all= TRUE)
    df_res <- rbind(df_baseline,df_herald)  
    df_res$exp <- factor(df_res$exp, levels=c("baseline","D"))
    stats_signif <-  df_res %>%
      wilcox_test(detlag_batch ~ exp, paired = TRUE, alternative = "greater")  %>%
      adjust_pvalue(method = "bonferroni") %>%
      add_significance("p.adj")
    stats_effsize <- df_res %>%
      group_map(~data.frame(parameter=name_par,lag=l, mdl=m, par=o,
                            effsize=interpret(rank_biserial(.x$detlag_batch, factor(.x$exp,levels=c("baseline",c("D"))), paired = TRUE, alternative = "greater"), rules = "funder2019")
      ),.keep=TRUE) %>%
      bind_rows()
    exp_par_stats[[name_par]][[as.character(p)]] <- merge(stats_signif,stats_effsize, sort = F)
  }
  exp_par_stats[[name_par]] <- do.call(rbind,exp_par_stats[[name_par]])
}

exp_par_stats_df <- do.call(rbind,exp_par_stats)
#View(exp_par_stats_df)

# Getting data.frame of lag differences ........................................

df_diff <- data.frame()
df_baseline <- exp_res_df  %>%
  filter(exp %in% c("baseline") & var == "high" & mem==27)

df_herald <- exp_res_df  %>%
  filter(exp %in% c("D") & var == "high" & mem==27) %>%
  distinct(curve,var,seed,mem,lag,mdl,par, .keep_all= TRUE)
par_comb <- list(lag=unique(df_herald$lag),mdl=unique(df_herald$mdl),par=unique(df_herald$par))

for(name_par in names(par_comb)){
  for(p in par_comb[[name_par]]){
    l <- ifelse(name_par=="lag",p,27)
    m <- ifelse(name_par=="mdl",p,27)
    o <- ifelse(name_par=="par",p,27)
    
    df_herald <- exp_res_df  %>%
      filter(exp %in% c("D") & var == "high" & mem==27 & lag==l & mdl==m & par==o) %>%
      distinct(curve,var,seed,mem,lag,mdl,par, .keep_all= TRUE)
    
    df_res <- inner_join(df_baseline,df_herald, by = c("curve","var","seed","mem")) %>% 
      mutate(parameter = name_par,
             diff_detlag_batch = detlag_batch.x - detlag_batch.y,
             diff_detlag_obs = detlag_obs.x - detlag_obs.y,
             diff_detlag_time = detlag_time.x - detlag_time.y)
    df_diff <- rbind(df_diff,df_res)
  }
}

df_diff <- inner_join(df_diff,exp_par_stats_df, by = join_by(lag.y == lag, mdl.y == mdl, par.y == par, parameter))
#View(df_diff)

# Plot .........................................................................
# create the discrete variables using cut
df_diff$effsize_cut <- cut(df_diff$effsize.r_rank_biserial, 
                           labels= c("Very large (neg)", "Large (neg)", "Medium (neg)", "Small (neg)", "Very small (neg)", "Tiny (neg)", "Tiny (pos)", "Very small (pos)", "Small (pos)", "Medium (pos)", "Large (pos)", "Very large (pos)"),
                           breaks = c(-1, -0.4, -0.3, -0.2, -0.1, -0.05, 0, 0.05, 0.1, 0.2, 0.3, 0.4, 1),
                           include.lowest = TRUE, right = FALSE)
levels_count <- length(levels(df_diff$effsize_cut))

# Create the colors scale coresponded to the levels in cut
color_scales_fn <- colorRampPalette(c("red", "#007FFF"))
manual_color <- color_scales_fn(levels_count)
names(manual_color) <- levels(df_diff$effsize_cut)
# Create the sizes scale 
manual_size <-  c(3, 2.5, 2, 1.5, 1, 0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3)
names(manual_size) <- levels(df_diff$effsize_cut)

# New facet label names for dose variable
par_labels <- c(par="MA order",lag="Lagged prediction horizon")
#par_labels <- c(lag="Lagged prediction horizon",mdl="Model size",par="MA order")



bp_parameters <- df_diff %>%
  filter(lag.y>=3) %>%
  pivot_longer(c(lag.y, mdl.y, par.y)) %>%
  group_by(parameter) %>%
  filter(parameter!="mdl", name==paste0(parameter,".y")) %>%
  #filter(name==paste0(parameter,".y")) %>%
  ggplot(aes(x = factor(value,levels=c(unique(value))))) +
  facet_wrap(~factor(parameter, levels=c('par', 'lag'),labels=par_labels), scales = "free_x",ncol=1) +
  #facet_wrap(~factor(parameter, levels=c('lag', 'mdl', 'par'),labels=par_labels), scales = "free_x",ncol=1) +
  geom_boxplot(aes(y=diff_detlag_batch,fill=effsize_cut)) +
  geom_hline(linetype = "dashed",yintercept=0, size = 0.5, color = "gray28") +
  scale_fill_manual(values=manual_color,name="Effect size",limits = force) +
  labs(y = "Detection lag differences")+
  scale_x_discrete(name = "Parameters",labels=c("27"=expression(bold("27")), parse=TRUE))+
  theme_bw()+  
  geom_text(aes(y=max(diff_detlag_batch,na.rm=TRUE),label = p.adj.signif), position = position_dodge(width=0.9),vjust=-0.7)+
  scale_y_continuous(expand = expansion(mult = c(0.2, 0.5)))

sp_parameters <- df_diff %>%
  filter(lag.y>=3) %>%
  pivot_longer(c(lag.y, mdl.y, par.y)) %>%
  group_by(parameter) %>%
  filter(parameter!="mdl", name==paste0(parameter,".y")) %>%
  #filter(name==paste0(parameter,".y")) %>%
  ggplot(aes(x = factor(value,levels=c(unique(value))))) +
  geom_hline(linetype = "dashed",yintercept=0, size = 0.5, color = "gray28") +
  geom_path(aes(y=effsize.r_rank_biserial, group = 1),color="grey",size=1) +
  geom_point(aes(y=effsize.r_rank_biserial, fill=effsize_cut, color=effsize_cut, size=effsize_cut),shape=21) +
  facet_wrap(~factor(parameter, levels=c('par', 'lag'),labels=par_labels), scales = "free_x",ncol=1) +
  #facet_wrap(~factor(parameter, levels=c('lag', 'mdl', 'par'),labels=par_labels), scales = "free_x",ncol=1) +
  scale_size_manual(values = manual_size,name="Effect size",limits = force) +
  scale_color_manual(values = manual_color,name="Effect size",limits = force) +
  scale_fill_manual(values=manual_color,name="Effect size",limits = force) +
  scale_y_continuous(position = "right", name = "Effect size",expand = expansion(mult = c(0.2, 0.5)))+ 
  scale_x_discrete(name = "Parameters",labels=c("27"=expression(bold("27")), parse=TRUE))+
  theme_bw()+  
  geom_text(aes(y=max(effsize.r_rank_biserial),label = p.adj.signif), position = position_dodge(width=0.9),vjust=-0.7)

ggarrange(bp_parameters, sp_parameters,common.legend = FALSE, legend = "bottom")
ggsave("sensibility_parameters.pdf", width = 23, height = 13, units = "cm")



## Memory ----------------------------------------------------------------------

exp_mem_stats <- list()
for(m in unique(exp_res_df_mem$mem)){
  df_baseline <- exp_res_df_mem  %>%
    filter(exp %in% c("baseline") & var == "high" & mem==m)
  df_herald <- exp_res_df_mem  %>%
    filter(exp %in% c("D") & var == "high" & lag==27 & mdl==27 & par==27 & mem==m) %>%
    distinct(curve,var,seed,mem,lag,mdl,par, .keep_all= TRUE)
  df_res <- rbind(df_baseline,df_herald)  
  df_res$exp <- factor(df_res$exp, levels=c("baseline","D"))
  stats_signif <-  df_res %>%
    wilcox_test(detlag_batch ~ exp, paired = TRUE, alternative = "greater")  %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")
  stats_effsize <- df_res %>%
    group_map(~data.frame(mem=m,
                          effsize=interpret(rank_biserial(.x$detlag_batch, factor(.x$exp,levels=c("baseline",c("D"))), paired = TRUE, alternative = "greater"), rules = "funder2019")
    ),.keep=TRUE) %>%
    bind_rows()
  exp_mem_stats[[as.character(m)]] <- merge(stats_signif,stats_effsize, sort = F)
}
exp_mem_stats_df <- do.call(rbind,exp_mem_stats)
#View(exp_mem_stats_df)

# Getting data.frame of lag differences ........................................

df_diff <- data.frame()

for(m in unique(exp_res_df_mem$mem)){
  df_baseline <- exp_res_df_mem  %>%
    filter(exp %in% c("baseline") & var == "high" & mem==m)
  df_herald <- exp_res_df_mem  %>%
    filter(exp %in% c("D") & var == "high" & lag==27 & mdl==27 & par==27 & mem==m) %>%
    distinct(curve,var,seed,mem,lag,mdl,par, .keep_all= TRUE)
  df_res <- inner_join(df_baseline,df_herald, by = c("curve","var","seed","mem")) %>% 
    mutate(diff_detlag_batch = detlag_batch.x - detlag_batch.y,
           diff_detlag_obs = detlag_obs.x - detlag_obs.y,
           diff_detlag_time = detlag_time.x - detlag_time.y)
  df_diff <- rbind(df_diff,df_res)
}
df_diff <- inner_join(df_diff,exp_mem_stats_df, by = join_by(mem == mem))
#View(df_diff)


# Plot .........................................................................
# create the discrete variables using cut
df_diff$effsize_cut <- cut(df_diff$effsize.r_rank_biserial, 
                           labels= c("Very large (neg)", "Large (neg)", "Medium (neg)", "Small (neg)", "Very small (neg)", "Tiny (neg)", "Tiny (pos)", "Very small (pos)", "Small (pos)", "Medium (pos)", "Large (pos)", "Very large (pos)"),
                           breaks = c(-1, -0.4, -0.3, -0.2, -0.1, -0.05, 0, 0.05, 0.1, 0.2, 0.3, 0.4, 1),
                           include.lowest = TRUE, right = FALSE)
levels_count <- length(levels(df_diff$effsize_cut))

# Create the colors scale coresponded to the levels in cut
color_scales_fn <- colorRampPalette(c("red", "#007FFF"))
manual_color <- color_scales_fn(levels_count)
names(manual_color) <- levels(df_diff$effsize_cut)
# Create the sizes scale 
manual_size <-  c(3, 2.5, 2, 1.5, 1, 0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3)
names(manual_size) <- levels(df_diff$effsize_cut)

bp_memory <- df_diff %>%
  ggplot(aes(x = factor(mem))) +
  geom_boxplot(aes(y=diff_detlag_batch,fill=effsize_cut)) +
  geom_hline(linetype = "dashed",yintercept=0, size = 0.5, color = "gray28") +
  scale_fill_manual(values=manual_color,name="Effect size",limits = force) +
  labs(y = "Detection lag differences")+
  scale_x_discrete(name = "Memory sizes",labels=c("27"=expression(bold("27")), parse=TRUE))+
  theme_bw()+  
  geom_text(aes(y=max(diff_detlag_batch,na.rm=TRUE),label = p.adj.signif), position = position_dodge(width=0.9),vjust=-0.7)+
  scale_y_continuous(expand = expansion(mult = c(0.2, 0.2)))

sp_memory <- df_diff %>%
  ggplot(aes(x = factor(mem))) +
  geom_hline(linetype = "dashed",yintercept=0, size = 0.5, color = "gray28") +
  geom_path(aes(y=effsize.r_rank_biserial, group = 1),color="grey",size=1) +
  geom_point(aes(y=effsize.r_rank_biserial, fill=effsize_cut, color=effsize_cut, size=effsize_cut),shape=21) +
  scale_size_manual(values = manual_size,name="Effect size",limits = force) +
  scale_color_manual(values = manual_color,name="Effect size",limits = force) +
  scale_fill_manual(values=manual_color,name="Effect size",limits = force) +
  scale_y_continuous(position = "right", name = "Effect size",expand = expansion(mult = c(0.2, 0.2)))+ 
  scale_x_discrete(name = "Memory sizes",labels=c("27"=expression(bold("27")), parse=TRUE))+
  theme_bw()+  
  geom_text(aes(y=max(effsize.r_rank_biserial),label = p.adj.signif), position = position_dodge(width=0.9),vjust=-0.7)

ggarrange(bp_memory, sp_memory,common.legend = TRUE, legend = "bottom")



## Number of events + memory ---------------------------------------------------

exp_mem_stats_lags <- list()
for(nevt in 1:length(unique(exp_res_df_mem_lags$event_idx))){
  df_baseline <- exp_res_df_mem_lags  %>%
    filter(exp %in% c("baseline") & var == "high" & event_idx %in% unique(event_idx)[1:nevt])%>%
    distinct(curve,var,seed,mem,lag,mdl,par,event_idx, .keep_all= TRUE)
  df_herald <- exp_res_df_mem_lags  %>%
    filter(exp %in% c("D") & var == "high" & lag==27 & mdl==27 & par==27 & event_idx %in% unique(event_idx)[1:nevt]) %>%
    distinct(curve,var,seed,mem,lag,mdl,par,event_idx, .keep_all= TRUE)
  df_res <- rbind(df_baseline,df_herald)  
  df_res$exp <- factor(df_res$exp, levels=c("baseline","D"))
  stats_signif <-  df_res %>%
    group_by(curve,var,seed,mem,event_idx) %>%
    filter(n() > 1 & !any(is.na(detlag_batch))) %>%
    group_by(mem) %>%
    wilcox_test(detlag_batch ~ exp, paired = TRUE, alternative = "greater")  %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")
  stats_effsize <- df_res %>%
    group_by(curve,var,seed,mem,event_idx) %>%
    filter(n() > 1 & !any(is.na(detlag_batch))) %>%
    group_by(mem) %>%
    group_map(~data.frame(mem=unique(.x$mem), n_evt= nevt,
                          effsize=interpret(rank_biserial(.x$detlag_batch, factor(.x$exp,levels=c("baseline",c("D"))), paired = TRUE, alternative = "greater"), rules = "funder2019")
    ),.keep=TRUE) %>%
    bind_rows()

  exp_mem_stats_lags[[nevt]] <- merge(stats_signif,stats_effsize, sort = F)
}
exp_mem_stat_lags_df <- do.call(rbind,exp_mem_stats_lags)
#View(exp_mem_stat_lags_df)


# Getting data.frame of lag differences ........................................

df_diff <- data.frame()

for(nevt in 1:length(unique(exp_res_df_mem_lags$event_idx))){
  df_baseline <- exp_res_df_mem_lags  %>%
    filter(exp %in% c("baseline") & var == "high" & event_idx %in% unique(event_idx)[1:nevt])%>%
    distinct(curve,var,seed,mem,lag,mdl,par,event_idx, .keep_all= TRUE)
  df_herald <- exp_res_df_mem_lags  %>%
    filter(exp %in% c("D") & var == "high" & lag==27 & mdl==27 & par==27 & event_idx %in% unique(event_idx)[1:nevt]) %>%
    distinct(curve,var,seed,mem,lag,mdl,par,event_idx, .keep_all= TRUE)
  
  df_res <- inner_join(df_baseline,df_herald, by = c("curve","var","seed","mem","event_idx")) %>% 
    mutate(n_evt=nevt,
           diff_detlag_batch = detlag_batch.x - detlag_batch.y,
           diff_detlag_obs = detlag_obs.x - detlag_obs.y,
           diff_detlag_time = detlag_time.x - detlag_time.y)
  df_diff <- rbind(df_diff,df_res)
}
df_diff <- inner_join(df_diff,exp_mem_stat_lags_df, by = join_by(mem == mem, n_evt==n_evt))
#View(df_diff)


# Plot .........................................................................
# create the discrete variables using cut
df_diff$effsize_cut <- cut(df_diff$effsize.r_rank_biserial, 
                           labels= c("Very large (neg)", "Large (neg)", "Medium (neg)", "Small (neg)", "Very small (neg)", "Tiny (neg)", "Tiny (pos)", "Very small (pos)", "Small (pos)", "Medium (pos)", "Large (pos)", "Very large (pos)"),
                           breaks = c(-1, -0.4, -0.3, -0.2, -0.1, -0.05, 0, 0.05, 0.1, 0.2, 0.3, 0.4, 1),
                           include.lowest = TRUE, right = FALSE)
levels_count <- length(levels(df_diff$effsize_cut))

# Create the colors scale coresponded to the levels in cut
color_scales_fn <- colorRampPalette(c("red", "#007FFF"))
manual_color <- color_scales_fn(levels_count)
names(manual_color) <- levels(df_diff$effsize_cut)
# Create the sizes scale 
manual_size <-  c(3, 2.5, 2, 1.5, 1, 0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3)
names(manual_size) <- levels(df_diff$effsize_cut)

# New facet label names for dose variable
par_labels <- c("1"="One event","2"="Two events","3"="Three events")

bp_nevt_memory <- df_diff %>%
  ggplot(aes(x = factor(mem))) +
  geom_boxplot(aes(y=diff_detlag_batch,fill=effsize_cut)) +
  facet_wrap(~factor(n_evt, levels=c(1, 2, 3),labels=par_labels), scales = "free_x",ncol=1) +
  geom_hline(linetype = "dashed",yintercept=0, size = 0.5, color = "gray28") +
  scale_fill_manual(values=manual_color,name="Effect size",limits = force) +
  labs(y = "Detection lag differences")+
  scale_x_discrete(name = "Memory sizes",labels=c("27"=expression(bold("27")), parse=TRUE))+
  theme_bw()+  
  geom_text(aes(y=max(diff_detlag_batch,na.rm=TRUE),label = p.adj.signif), position = position_dodge(width=0.9),vjust=-0.7)+
  scale_y_continuous(expand = expansion(mult = c(0.2, 0.5)))

sp_nevt_memory <- df_diff %>%
  arrange(mem) %>% 
  ggplot(aes(x = factor(mem))) +
  geom_hline(linetype = "dashed",yintercept=0, size = 0.5, color = "gray28") +
  geom_path(aes(y=effsize.r_rank_biserial, group = 1),color="grey",size=1) +
  geom_point(aes(y=effsize.r_rank_biserial, fill=effsize_cut, color=effsize_cut, size=effsize_cut),shape=21) +
  facet_wrap(~factor(n_evt, levels=c(1, 2, 3),labels=par_labels), scales = "free_x",ncol=1) +
  scale_size_manual(values = manual_size,name="Effect size",limits = force) +
  scale_color_manual(values = manual_color,name="Effect size",limits = force) +
  scale_fill_manual(values=manual_color,name="Effect size",limits = force) +
  scale_y_continuous(position = "right", name = "Effect size",expand = expansion(mult = c(0.2, 0.5)))+ 
  scale_x_discrete(name = "Memory sizes",labels=c("27"=expression(bold("27")), parse=TRUE))+
  theme_bw()+  
  geom_text(aes(y=max(effsize.r_rank_biserial),label = p.adj.signif), position = position_dodge(width=0.9),vjust=-0.7)

ggarrange(bp_nevt_memory, sp_nevt_memory,common.legend = TRUE, legend = "bottom")
ggsave("sensibility_memory.pdf", width = 20, height = 16, units = "cm")
