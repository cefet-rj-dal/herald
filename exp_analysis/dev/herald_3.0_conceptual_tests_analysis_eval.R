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

#Contains results for baseline and D for memory sensibility analysis
exp_res_df_mem_eval <- readRDS(paste0(getwd(),"/dev/plots/exp_res_df_0D_mem_eval.rds"))

exp_res_df_mem_eval$curve <- factor(exp_res_df_mem_eval$curve,levels=curves)
exp_res_df_mem_eval$var <- factor(exp_res_df_mem_eval$var, ordered = TRUE, levels=names(volatility))
exp_res_df_mem_eval$seed <- factor(exp_res_df_mem_eval$seed, levels=seeds)
exp_res_df_mem_eval$exp <- factor(exp_res_df_mem_eval$exp, ordered = TRUE, levels=c("baseline","D"))

#Contains example plots for baseline and D
exp_res_plots_eval <- readRDS(paste0(getwd(),"/dev/plots/exp_res_plots_eval.rds"))



# Baseline versus Herald =======================================================

## Statistical tests -----------------------------------------------------------

df_baseline <- exp_res_df_mem_eval  %>%
  filter(exp %in% c("baseline") & mem==27)%>%
  distinct(curve,var,seed,mem,lag,mdl,par, .keep_all= TRUE)
df_herald <- exp_res_df_mem_eval  %>%
  filter(exp %in% c("D") & mem==27 & lag==27 & mdl==27 & par==27) %>%
  distinct(curve,var,seed,mem,lag,mdl,par, .keep_all= TRUE)

df_res <- rbind(df_baseline,df_herald)
df_res$curve <- factor(df_res$curve,levels=curves)
df_res$var <- factor(df_res$var, ordered = TRUE,levels=names(volatility))
df_res$exp <- factor(df_res$exp,levels=c("baseline","D"))

cols_soft_less <- paste(c("TPs","TNs","precision","recall","F1"),"_soft",sep="")
cols_soft_greater <- paste(c("FPs","FNs"),"_soft",sep="")

stats_soft <- list()
for(col in cols_soft_less){
  stats_signif <- df_res  %>%
    wilcox_test(as.formula(paste0(col," ~ exp")), paired = TRUE, alternative = "less")  %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")%>% 
    add_xy_position(x = "exp")
  stats_effsize <- df_res %>%
    group_map(~data.frame(effsize=interpret(rank_biserial(.x[[col]], factor(.x$exp,levels=c("baseline","D")), paired = TRUE, alternative = "less"), rules = "funder2019")
    ),.keep=TRUE) %>%
    bind_rows()
    stats_soft[[col]] <- merge(stats_signif,stats_effsize, sort = F)
}
for(col in cols_soft_greater){
  stats_signif <- df_res  %>%
    wilcox_test(as.formula(paste0(col," ~ exp")), paired = TRUE, alternative = "greater")  %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")%>% 
    add_xy_position(x = "exp")
  stats_effsize <- df_res %>%
    group_map(~data.frame(effsize=interpret(rank_biserial(.x[[col]], factor(.x$exp,levels=c("baseline","D")), paired = TRUE, alternative = "greater"), rules = "funder2019")
    ),.keep=TRUE) %>%
    bind_rows()
    stats_soft[[col]] <- merge(stats_signif,stats_effsize, sort = F)
}
herald_test_total <- do.call(rbind,stats_soft)
#View(herald_test_total)

 
stats_soft <- list()
for(col in cols_soft_less){
  stats_signif <-  df_res %>%
    group_by(curve,var,seed,mem) %>%
    filter(n() > 1) %>%
    group_by(var) %>%
    wilcox_test(as.formula(paste0(col," ~ exp")), paired = TRUE, alternative = "less")  %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")%>% 
    add_xy_position(x = "exp")
  stats_effsize <- df_res %>%
    group_by(curve,var,seed,mem) %>%
    filter(n() > 1) %>%
    group_by(var) %>%
    group_map(~data.frame(var=unique(.x$var), 
                          effsize=interpret(rank_biserial(.x[[col]], factor(.x$exp,levels=c("baseline","D")), paired = TRUE, alternative = "less"), rules = "funder2019")
    ),.keep=TRUE) %>%
    bind_rows()
    stats_soft[[col]] <- merge(stats_signif,stats_effsize, sort = F)
}
for(col in cols_soft_greater){
  stats_signif <-  df_res %>%
    group_by(curve,var,seed,mem) %>%
    filter(n() > 1) %>%
    group_by(var) %>%
    wilcox_test(as.formula(paste0(col," ~ exp")), paired = TRUE, alternative = "greater")  %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")%>% 
    add_xy_position(x = "exp")
  stats_effsize <- df_res %>%
    group_by(curve,var,seed,mem) %>%
    filter(n() > 1) %>%
    group_by(var) %>%
    group_map(~data.frame(var=unique(.x$var), 
                          effsize=interpret(rank_biserial(.x[[col]], factor(.x$exp,levels=c("baseline","D")), paired = TRUE, alternative = "greater"), rules = "funder2019")
    ),.keep=TRUE) %>%
    bind_rows()
    stats_soft[[col]] <- merge(stats_signif,stats_effsize, sort = F)
}
herald_test <- do.call(rbind,stats_soft)
#View(herald_test)


# Plot .........................................................................
# Box plot

# New facet label names for dose variable
par_labels_metrics <- c(recall_soft="SoftED recall",precision_soft="SoftED precision",F1_soft="SoftED F1")
par_labels_confmatrix <- c(TPs_soft="SoftED TP",FNs_soft="SoftED FN",FPs_soft="SoftED FP",TNs_soft="SoftED TN")

herald_test_total$name <- factor(herald_test_total$.y.,levels=names(c(par_labels_metrics,par_labels_confmatrix)))
herald_test$name <- factor(herald_test$.y.,levels=names(c(par_labels_metrics,par_labels_confmatrix)))
                                 
bp <- df_res %>%
  pivot_longer(c(unique(herald_test_total$name))) %>%
  filter(name %in% names(par_labels_metrics)) %>%
  ggplot(aes(x = factor(exp,levels=c("baseline","D")),y=value, group==exp)) +
  geom_boxplot(aes(fill=exp), show.legend = FALSE) +
  facet_wrap(~factor(name,levels=names(par_labels_metrics),labels=par_labels_metrics), scales = "free_y",ncol=3) +
  scale_fill_manual(values=c("white","#007FFF")) +
  labs(y = "Detection lags")+
  scale_x_discrete(name="Method",labels = c('EDfR','Herald'))+
  theme_bw()+ 
  stat_pvalue_manual(herald_test_total %>% filter(name %in% names(par_labels_metrics)), 
                     label = "Wilcoxon, p = {p.adj} {p.adj.signif},\n Effect = {round(-effsize.r_rank_biserial,2)} ({effsize.Interpretation})", tip.length = 0.01)+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.5)))

bp_var <- df_res %>%
  pivot_longer(c(unique(herald_test_total$name))) %>%
  filter(name %in% names(par_labels_metrics)) %>%
  ggplot(aes(x = factor(exp,levels=c("baseline","D")),y=value, group==exp)) +
  geom_boxplot(aes(fill=exp), show.legend = FALSE) +
  facet_grid(factor(var, levels=c('high', 'medium', 'low'), labels=c('High volatility', 'Medium volatility', 'Low volatility'))
             ~factor(name,levels=names(par_labels_metrics),labels=par_labels_metrics), scales = "free") +
  scale_fill_manual(values=c("white","#007FFF")) +
  labs(y = "Detection lags")+
  scale_x_discrete(name="Method",labels = c('EDfR','Herald'))+
  theme_bw()+ 
  stat_pvalue_manual(herald_test %>% filter(name %in% names(par_labels_metrics)),
                     label = "Wilcoxon, p = {p.adj} {p.adj.signif},\n Effect = {round(-effsize.r_rank_biserial,2)} ({effsize.Interpretation})", tip.length = 0.01)+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.6)))

ggarrange(bp, bp_var, nrow=2, heights = c(1, 2) ,common.legend = TRUE, legend = "bottom")
ggsave("herald_softed_stats.pdf", width = 20, height = 18, units = "cm")


bp <- df_res %>%
  pivot_longer(c(unique(herald_test_total$name))) %>%
  filter(name %in% names(par_labels_confmatrix)) %>%
  ggplot(aes(x = factor(exp,levels=c("baseline","D")),y=value, group==exp)) +
  geom_boxplot(aes(fill=exp), show.legend = FALSE) +
  facet_wrap(~factor(name,levels=names(par_labels_confmatrix),labels=par_labels_confmatrix), scales = "free_y",ncol=1) +
  scale_fill_manual(values=c("white","#007FFF")) +
  labs(y = "Detection lags")+
  scale_x_discrete(name="Method",labels = c('EDfR','Herald'))+
  theme_bw()+ 
  stat_pvalue_manual(herald_test_total %>% filter(name %in% names(par_labels_confmatrix)), 
                     label = "Wilcoxon, p = {p.adj} {p.adj.signif},\n Effect = {round(-effsize.r_rank_biserial,2)} ({effsize.Interpretation})", tip.length = 0.01)+
  scale_y_continuous(expand = expansion(mult = c(0.1, 1.2)))

bp_var <- df_res %>%
  pivot_longer(c(unique(herald_test_total$name))) %>%
  filter(name %in% names(par_labels_confmatrix)) %>%
  ggplot(aes(x = factor(exp,levels=c("baseline","D")),y=value, group==exp)) +
  geom_boxplot(aes(fill=exp), show.legend = FALSE) +
  facet_grid(factor(name,levels=names(par_labels_confmatrix),labels=par_labels_confmatrix)
             ~factor(var, levels=c('high', 'medium', 'low'), labels=c('High volatility', 'Medium volatility', 'Low volatility')), scales = "free") +
  scale_fill_manual(values=c("white","#007FFF")) +
  labs(y = "Detection lags")+
  scale_x_discrete(name="Method",labels = c('EDfR','Herald'))+
  theme_bw()+ 
  stat_pvalue_manual(herald_test %>% filter(name %in% names(par_labels_confmatrix)),
                     label = "Wilcoxon, p = {p.adj} {p.adj.signif},\n Effect = {round(-effsize.r_rank_biserial,2)} ({effsize.Interpretation})", tip.length = 0.01)+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.8)))

#ggarrange(bp, bp_var, widths = c(1, 2) ,common.legend = TRUE, legend = "bottom")



## Example plots ---------------------------------------------------------------

df_baseline <- exp_res_df_mem_eval  %>%
  filter(exp %in% c("baseline") & var == "high" & seed == 123 & mem==27) %>%
  distinct(curve,var,seed,mem,lag,mdl,par, .keep_all= TRUE)
df_herald <- exp_res_df_mem_eval  %>%
  filter(exp %in% c("D") & var == "high" & seed == 123 & mem==27 & lag==27 & mdl==27 & par==27) %>%
  distinct(curve,var,seed,mem,lag,mdl,par, .keep_all= TRUE)
df_res <- rbind(df_baseline,df_herald)

list_plots <- lapply(exp_res_plots_eval, function(curve) curve[c("baseline","D")])
labels <- names(unlist(list_plots,recursive=FALSE))
labels <- gsub("^.*\\.", "", labels)
labels <- gsub("D", "Herald", labels)
labels <- gsub("baseline", "EDfR", labels)

list_plots <- lapply(names(list_plots), function(crv) 
  lapply(names(list_plots[[crv]]), function(method) 
    list_plots[[crv]][[method]] + ggtitle(paste0("SoftED recall: ",round(df_res %>% filter(curve==crv,exp==method) %>% select(recall_soft),2),
                         ", SoftED precision: ",round(df_res %>% filter(curve==crv,exp==method) %>% select(precision_soft),2),
                         ", SoftED F1: ",round(df_res %>% filter(curve==crv,exp==method) %>% select(F1_soft),2))) 
    )
)
list_plots <- unlist(list_plots,recursive=FALSE)
list_plots <- lapply(1:length(list_plots),function(p) list_plots[[p]] + ggtitle(paste0(labels[p],"\n(",list_plots[[p]]$labels$title,")")))

#ggarrange(plotlist=list_plots, nrow=7,ncol=2)
plots <- list_plots[c(1,2,13,14)]
plots[[3]] <- plots[[3]]+ ggtitle(gsub("^.*\\\n", "", plots[[3]]$labels$title))
plots[[4]] <- plots[[4]]+ ggtitle(gsub("^.*\\\n", "", plots[[4]]$labels$title))

ggarrange(plotlist=plots, nrow=2, ncol=2, heights=c(1.1, 1))
ggsave("plots_eval.pdf", width = 28, height = 12, units = "cm")


# Sensibility analysis =========================================================

## Parameters ------------------------------------------------------------------

df_baseline <- exp_res_df_mem_eval  %>%
  filter(exp %in% c("baseline") & var == "high" & mem==27)%>%
  distinct(curve,var,seed,mem,lag,mdl,par, .keep_all= TRUE)
df_herald <- exp_res_df_mem_eval  %>%
  filter(exp %in% c("D") & var == "high" & mem==27) %>%
  distinct(curve,var,seed,mem,lag,mdl,par, .keep_all= TRUE)

par_comb <- list(lag=unique(df_herald$lag),mdl=unique(df_herald$mdl),par=unique(df_herald$par))

cols_soft_less <- paste(c("recall","precision","F1"),"_soft",sep="")

exp_par_stats <- list()
for(col in cols_soft_less){
  exp_par_stats[[col]] <- list()
  for(name_par in names(par_comb)){
    exp_par_stats[[col]][[name_par]] <- list()
    for(p in par_comb[[name_par]]){
      l <- ifelse(name_par=="lag",p,27)
      m <- ifelse(name_par=="mdl",p,27)
      o <- ifelse(name_par=="par",p,27)
      
      df_herald <- exp_res_df_mem_eval  %>%
        filter(exp %in% c("D") & var == "high" & mem==27 & lag==l & mdl==m & par==o) %>%
        distinct(curve,var,seed,mem,lag,mdl,par, .keep_all= TRUE)
      df_res <- rbind(df_baseline,df_herald)
      df_res$exp <- factor(df_res$exp, levels=c("baseline","D"))
      stats_signif <- df_res  %>%
        wilcox_test(as.formula(paste0(col," ~ exp")), paired = TRUE, alternative = "less")  %>%
        adjust_pvalue(method = "bonferroni") %>%
        add_significance("p.adj")%>% 
        add_xy_position(x = "exp")
      stats_effsize <- df_res %>%
        group_map(~data.frame(metric=col,parameter=name_par,lag=l, mdl=m, par=o,
                              effsize=interpret(rank_biserial(.x[[col]], factor(.x$exp,levels=c("baseline","D")), paired = TRUE, alternative = "less"), rules = "funder2019")
        ),.keep=TRUE) %>%
        bind_rows()
      exp_par_stats[[col]][[name_par]][[as.character(p)]] <- merge(stats_signif,stats_effsize, sort = F)
    }
    exp_par_stats[[col]][[name_par]] <- do.call(rbind,exp_par_stats[[col]][[name_par]])
  }
  exp_par_stats[[col]] <- do.call(rbind,exp_par_stats[[col]])
}
exp_par_stats_df <- do.call(rbind,exp_par_stats)
#View(exp_par_stats_df)


# Getting data.frame of lag differences ........................................

df_diff <- data.frame()
df_baseline <- exp_res_df_mem_eval  %>%
  filter(exp %in% c("baseline") & var == "high" & mem==27)%>%
  distinct(curve,var,seed,mem,lag,mdl,par, .keep_all= TRUE)
df_herald <- exp_res_df_mem_eval  %>%
  filter(exp %in% c("D") & var == "high" & mem==27) %>%
  distinct(curve,var,seed,mem,lag,mdl,par, .keep_all= TRUE)

par_comb <- list(lag=unique(df_herald$lag),mdl=unique(df_herald$mdl),par=unique(df_herald$par))

cols_soft_less <- paste(c("recall","precision","F1"),"_soft",sep="")

for(col in cols_soft_less){
  for(name_par in names(par_comb)){
    for(p in par_comb[[name_par]]){
      l <- ifelse(name_par=="lag",p,27)
      m <- ifelse(name_par=="mdl",p,27)
      o <- ifelse(name_par=="par",p,27)
      
      df_herald <- exp_res_df_mem_eval  %>%
        filter(exp %in% c("D") & var == "high" & mem==27 & lag==l & mdl==m & par==o) %>%
        distinct(curve,var,seed,mem,lag,mdl,par, .keep_all= TRUE)
      
      df_res <- inner_join(df_baseline,df_herald, by = c("curve","var","seed","mem")) %>% 
        mutate(parameter = name_par, col=col, 
               diff_recall_soft = recall_soft.y - recall_soft.x,
               diff_precision_soft = precision_soft.y - precision_soft.x,
               diff_F1_soft = F1_soft.y - F1_soft.x)
      df_diff <- rbind(df_diff,df_res)
    }
  }
}

df_diff <- inner_join(df_diff,exp_par_stats_df, by = join_by(lag.y == lag, mdl.y == mdl, par.y == par, parameter, col==.y.))
#View(df_diff)

# Plot .........................................................................
# create the discrete variables using cut
df_diff$effsize_cut <- cut(df_diff$effsize.r_rank_biserial, 
                           labels= rev(c("Very large (neg)", "Large (neg)", "Medium (neg)", "Small (neg)", "Very small (neg)", "Tiny (neg)", "Tiny (pos)", "Very small (pos)", "Small (pos)", "Medium (pos)", "Large (pos)", "Very large (pos)")),
                           breaks = rev(c(-1, -0.4, -0.3, -0.2, -0.1, -0.05, 0, 0.05, 0.1, 0.2, 0.3, 0.4, 1)),
                           include.lowest = TRUE, right = FALSE)
levels_count <- length(levels(df_diff$effsize_cut))

# Create the colors scale coresponded to the levels in cut
color_scales_fn <- colorRampPalette(c("#007FFF","red"))
manual_color <- color_scales_fn(levels_count)
names(manual_color) <- levels(df_diff$effsize_cut)
# Create the sizes scale 
manual_size <-  c(3, 2.5, 2, 1.5, 1, 0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3)
names(manual_size) <- levels(df_diff$effsize_cut)

# New facet label names
par_labels <- c(par="MA order",lag="Lagged prediction horizon")
#par_labels <- c(lag="Lagged prediction horizon",mdl="Model size",par="MA order")
par_labels_metrics <- c(recall_soft="SoftED recall",precision_soft="SoftED precision",F1_soft="SoftED F1")


bp_recall <- df_diff %>%
  filter(lag.y>=3,col=="recall_soft") %>%
  pivot_longer(c(lag.y, mdl.y, par.y)) %>%
  group_by(parameter) %>%
  filter(parameter!="mdl", name==paste0(parameter,".y")) %>%
  #filter(name==paste0(parameter,".y")) %>%
  ggplot(aes(x = factor(value,levels=c(unique(value))))) +
  facet_grid("SoftED recall" ~factor(parameter, levels=c('par', 'lag'),labels=par_labels), scales = "free_x") +
  #facet_grid("SoftED recall" ~factor(parameter, levels=c('lag', 'mdl', 'par'),labels=par_labels), scales = "free_x") +
  geom_boxplot(aes(y=diff_recall_soft,fill=effsize_cut)) +
  geom_hline(linetype = "dashed",yintercept=0, size = 0.5, color = "gray28") +
  scale_fill_manual(values=manual_color,name="Effect size",limits = force) +
  labs(y = "Detection lag differences")+
  scale_x_discrete(name = "Parameters",labels=c("27"=expression(bold("27")), parse=TRUE))+
  theme_bw()+  
  geom_text(aes(y=max(diff_recall_soft,na.rm=TRUE),label = p.adj.signif), position = position_dodge(width=0.9),vjust=-0.7)+
  scale_y_continuous(expand = expansion(mult = c(0.2, 0.5)),labels = function(x) format(x, nsmall = 2))

bp_precision <- df_diff %>%
  filter(lag.y>=3,col=="precision_soft") %>%
  pivot_longer(c(lag.y, mdl.y, par.y)) %>%
  group_by(parameter) %>%
  filter(parameter!="mdl", name==paste0(parameter,".y")) %>%
  #filter(name==paste0(parameter,".y")) %>%
  ggplot(aes(x = factor(value,levels=c(unique(value))))) +
  facet_grid("SoftED precision" ~factor(parameter, levels=c('par', 'lag'),labels=par_labels), scales = "free_x") +
  #facet_grid("SoftED precision" ~factor(parameter, levels=c('lag', 'mdl', 'par'),labels=par_labels), scales = "free_x") +
  geom_boxplot(aes(y=diff_precision_soft,fill=effsize_cut)) +
  geom_hline(linetype = "dashed",yintercept=0, size = 0.5, color = "gray28") +
  scale_fill_manual(values=manual_color,name="Effect size",limits = force) +
  labs(y = "Detection lag differences")+
  scale_x_discrete(name = "Parameters",labels=c("27"=expression(bold("27")), parse=TRUE))+
  theme_bw()+  
  geom_text(aes(y=max(diff_precision_soft,na.rm=TRUE),label = p.adj.signif), position = position_dodge(width=0.9),vjust=-0.7)+
  scale_y_continuous(expand = expansion(mult = c(0.2, 0.5)),labels = function(x) format(x, nsmall = 2))

bp_F1 <- df_diff %>%
  filter(lag.y>=3,col=="F1_soft") %>%
  pivot_longer(c(lag.y, mdl.y, par.y)) %>%
  group_by(parameter) %>%
  filter(parameter!="mdl", name==paste0(parameter,".y")) %>%
  #filter(name==paste0(parameter,".y")) %>%
  ggplot(aes(x = factor(value,levels=c(unique(value))))) +
  facet_grid("SoftED F1" ~factor(parameter, levels=c('par', 'lag'),labels=par_labels), scales = "free_x") +
  #facet_grid("SoftED F1" ~factor(parameter, levels=c('lag', 'mdl', 'par'),labels=par_labels), scales = "free_x") +
  geom_boxplot(aes(y=diff_F1_soft,fill=effsize_cut)) +
  geom_hline(linetype = "dashed",yintercept=0, size = 0.5, color = "gray28") +
  scale_fill_manual(values=manual_color,name="Effect size",limits = force) +
  labs(y = "Detection lag differences")+
  scale_x_discrete(name = "Parameters",labels=c("27"=expression(bold("27")), parse=TRUE))+
  theme_bw()+  
  geom_text(aes(y=max(diff_F1_soft,na.rm=TRUE),label = p.adj.signif), position = position_dodge(width=0.9),vjust=-0.7)+
  scale_y_continuous(expand = expansion(mult = c(0.2, 0.5)),labels = function(x) format(x, nsmall = 2))

remove_x <- theme(
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.title.x = element_blank()
)
require(grid)
bp <- ggarrange(bp_recall+ ylab(""), bp_precision+ ylab(""), bp_F1+ ylab(""),common.legend = F, nrow=3,legend = "bottom")
bp <- bp %>% annotate_figure(left = textGrob("Detection lag differences", rot = 90, vjust = 1, gp = gpar(cex = 1.0)))

sp_recall <- df_diff %>%
  filter(lag.y>=3,col=="recall_soft") %>%
  pivot_longer(c(lag.y, mdl.y, par.y)) %>%
  group_by(parameter) %>%
  filter(parameter!="mdl", name==paste0(parameter,".y")) %>%
  #filter(name==paste0(parameter,".y")) %>%
  ggplot(aes(x = factor(value,levels=c(unique(value))))) +
  geom_hline(linetype = "dashed",yintercept=0, size = 0.5, color = "gray28") +
  geom_path(aes(y=-effsize.r_rank_biserial, group = 1),color="grey",size=1) +
  geom_point(aes(y=-effsize.r_rank_biserial, fill=effsize_cut, color=effsize_cut, size=effsize_cut),shape=21) +
  facet_grid("SoftED recall" ~factor(parameter, levels=c('par', 'lag'),labels=par_labels), scales = "free_x", switch = "y") +
  #facet_grid("SoftED recall" ~factor(parameter, levels=c('lag', 'mdl', 'par'),labels=par_labels), scales = "free_x") +
  scale_size_manual(values = manual_size,name="Effect size",limits = force) +
  scale_color_manual(values = manual_color,name="Effect size",limits = force) +
  scale_fill_manual(values=manual_color,name="Effect size",limits = force) +
  scale_y_continuous(position = "right",expand = expansion(mult = c(0.2, 0.5)))+ 
  labs(y = "Effect size")+
  scale_x_discrete(name = "Parameters",labels=c("27"=expression(bold("27")), parse=TRUE))+
  theme_bw()+  
  geom_text(aes(y=max(-effsize.r_rank_biserial),label = p.adj.signif), position = position_dodge(width=0.9),vjust=-0.7)+
  theme(plot.margin=margin(5.5,5.5,5.5,0, "pt"))

sp_precision <- df_diff %>%
  filter(lag.y>=3,col=="precision_soft") %>%
  pivot_longer(c(lag.y, mdl.y, par.y)) %>%
  group_by(parameter) %>%
  filter(parameter!="mdl", name==paste0(parameter,".y")) %>%
  #filter(name==paste0(parameter,".y")) %>%
  ggplot(aes(x = factor(value,levels=c(unique(value))))) +
  geom_hline(linetype = "dashed",yintercept=0, size = 0.5, color = "gray28") +
  geom_path(aes(y=-effsize.r_rank_biserial, group = 1),color="grey",size=1) +
  geom_point(aes(y=-effsize.r_rank_biserial, fill=effsize_cut, color=effsize_cut, size=effsize_cut),shape=21) +
  facet_grid("SoftED precision" ~factor(parameter, levels=c('par', 'lag'),labels=par_labels), scales = "free_x", switch = "y") +
  #facet_grid("SoftED precision" ~factor(parameter, levels=c('lag', 'mdl', 'par'),labels=par_labels), scales = "free_x") +
  scale_size_manual(values = manual_size,name="Effect size",limits = force) +
  scale_color_manual(values = manual_color,name="Effect size",limits = force) +
  scale_fill_manual(values=manual_color,name="Effect size",limits = force) +
  scale_y_continuous(position = "right",expand = expansion(mult = c(0.2, 0.5)))+ 
  labs(y = "Effect size")+
  scale_x_discrete(name = "Parameters",labels=c("27"=expression(bold("27")), parse=TRUE))+
  theme_bw()+  
  geom_text(aes(y=max(-effsize.r_rank_biserial),label = p.adj.signif), position = position_dodge(width=0.9),vjust=-0.7)+
  theme(plot.margin=margin(5.5,5.5,5.5,0, "pt"))

sp_F1 <- df_diff %>%
  filter(lag.y>=3,col=="F1_soft") %>%
  pivot_longer(c(lag.y, mdl.y, par.y)) %>%
  group_by(parameter) %>%
  filter(parameter!="mdl", name==paste0(parameter,".y")) %>%
  #filter(name==paste0(parameter,".y")) %>%
  ggplot(aes(x = factor(value,levels=c(unique(value))))) +
  geom_hline(linetype = "dashed",yintercept=0, size = 0.5, color = "gray28") +
  geom_path(aes(y=-effsize.r_rank_biserial, group = 1),color="grey",size=1) +
  geom_point(aes(y=-effsize.r_rank_biserial, fill=effsize_cut, color=effsize_cut, size=effsize_cut),shape=21) +
  facet_grid("SoftED F1" ~factor(parameter, levels=c('par', 'lag'),labels=par_labels), scales = "free_x", switch = "y") +
  #facet_grid("SoftED F1" ~factor(parameter, levels=c('lag', 'mdl', 'par'),labels=par_labels), scales = "free_x") +
  scale_size_manual(values = manual_size,name="Effect size",limits = force) + 
  scale_color_manual(values = manual_color,name="Effect size",limits = force) +
  scale_fill_manual(values=manual_color,name="Effect size",limits = force) +
  scale_y_continuous(position = "right",expand = expansion(mult = c(0.2, 0.5)))+ 
  labs(y = "Effect size")+
  scale_x_discrete(name = "Parameters",labels=c("27"=expression(bold("27")), parse=TRUE))+
  theme_bw()+  
  geom_text(aes(y=max(-effsize.r_rank_biserial),label = p.adj.signif), position = position_dodge(width=0.9),vjust=-0.7)+
  theme(plot.margin=margin(5.5,5.5,5.5,0, "pt"))

remove_x <- theme(
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.title.x = element_blank()
)
require(grid)
sp <- ggarrange(sp_recall+ylab("")+theme(strip.text.y = element_blank()), 
                sp_precision+ylab("")+theme(strip.text.y = element_blank()),
                sp_F1+ylab("")+theme(strip.text.y = element_blank()),common.legend = FALSE, nrow=3,legend = "bottom")
sp <- sp %>% annotate_figure(right = textGrob("Effect size", rot = 270, vjust = 1, gp = gpar(cex = 1.0)))


ggarrange(bp, 
          ggplot()+theme_void(), 
          sp, common.legend = TRUE, legend = "bottom",nrow=1,widths = c(1, -0.0145, 1))
ggsave("sensibility_parameters_eval.pdf", width = 25, height = 20, units = "cm")


## Memory ----------------------------------------------------------------------

cols_soft_less <- paste(c("recall","precision","F1"),"_soft",sep="")
exp_mem_stats <- list()
for(col in cols_soft_less){
  exp_mem_stats[[col]] <- list()
  for(m in unique(exp_res_df_mem_eval$mem)){
    
    df_baseline <- exp_res_df_mem_eval  %>%
      filter(exp %in% c("baseline") & var == "high" & mem==m)%>%
      distinct(curve,var,seed,mem,lag,mdl,par, .keep_all= TRUE)
    df_herald <- exp_res_df_mem_eval  %>%
      filter(exp %in% c("D") & var == "high" & lag==27 & mdl==27 & par==27 & mem==m) %>%
      distinct(curve,var,seed,mem,lag,mdl,par, .keep_all= TRUE)
    df_res <- rbind(df_baseline,df_herald)  
    df_res$exp <- factor(df_res$exp, levels=c("baseline","D"))
    
    stats_signif <- df_res  %>%
      wilcox_test(as.formula(paste0(col," ~ exp")), paired = TRUE, alternative = "less")  %>%
      adjust_pvalue(method = "bonferroni") %>%
      add_significance("p.adj")
    stats_effsize <- df_res %>%
      group_map(~data.frame(metric=col,mem=m,
                            effsize=interpret(rank_biserial(.x[[col]], factor(.x$exp,levels=c("baseline","D")), paired = TRUE, alternative = "less"), rules = "funder2019")
      ),.keep=TRUE) %>%
      bind_rows()
    exp_mem_stats[[col]][[as.character(m)]] <- merge(stats_signif,stats_effsize, sort = F)
  }
  exp_mem_stats[[col]] <- do.call(rbind,exp_mem_stats[[col]])
}
exp_mem_stats_df <- do.call(rbind,exp_mem_stats)
#View(exp_mem_stats_df)


# Getting data.frame of lag differences ........................................

df_diff <- data.frame()

for(col in cols_soft_less){
  for(m in unique(exp_res_df_mem_eval$mem)){
    df_baseline <- exp_res_df_mem_eval  %>%
      filter(exp %in% c("baseline") & var == "high" & mem==m)%>%
      distinct(curve,var,seed,mem,lag,mdl,par, .keep_all= TRUE)
    df_herald <- exp_res_df_mem_eval  %>%
      filter(exp %in% c("D") & var == "high" & lag==27 & mdl==27 & par==27 & mem==m) %>%
      distinct(curve,var,seed,mem,lag,mdl,par, .keep_all= TRUE)
    
    df_res <- inner_join(df_baseline,df_herald, by = c("curve","var","seed","mem")) %>% 
      mutate(mem = m, col=col, 
             diff_recall_soft = recall_soft.y - recall_soft.x,
             diff_precision_soft = precision_soft.y - precision_soft.x,
             diff_F1_soft = F1_soft.y - F1_soft.x)
    df_diff <- rbind(df_diff,df_res)
  }
}
df_diff <- inner_join(df_diff,exp_mem_stats_df, by = join_by(mem == mem, col==.y.))
#View(df_diff)


# Plot .........................................................................
# create the discrete variables using cut
df_diff$effsize_cut <- cut(df_diff$effsize.r_rank_biserial, 
                           labels= rev(c("Very large (neg)", "Large (neg)", "Medium (neg)", "Small (neg)", "Very small (neg)", "Tiny (neg)", "Tiny (pos)", "Very small (pos)", "Small (pos)", "Medium (pos)", "Large (pos)", "Very large (pos)")),
                           breaks = rev(c(-1, -0.4, -0.3, -0.2, -0.1, -0.05, 0, 0.05, 0.1, 0.2, 0.3, 0.4, 1)),
                           include.lowest = TRUE, right = FALSE)
levels_count <- length(levels(df_diff$effsize_cut))

# Create the colors scale coresponded to the levels in cut
color_scales_fn <- colorRampPalette(c("#007FFF","red"))
manual_color <- color_scales_fn(levels_count)
names(manual_color) <- levels(df_diff$effsize_cut)
# Create the sizes scale 
manual_size <-  c(3, 2.5, 2, 1.5, 1, 0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3)
names(manual_size) <- levels(df_diff$effsize_cut)

par_labels_metrics <- c(diff_recall_soft="SoftED recall",diff_precision_soft="SoftED precision",diff_F1_soft="SoftED F1")

bp_memory <- df_diff %>%
  pivot_longer(c(diff_recall_soft, diff_precision_soft, diff_F1_soft)) %>%
  group_by(col) %>%
  filter(name==paste0("diff_",col)) %>%
  ggplot(aes(x = factor(mem,levels=c(unique(mem))))) +
  facet_wrap(~factor(name, levels=c('diff_recall_soft', 'diff_precision_soft', 'diff_F1_soft'),labels=par_labels_metrics), scales = "free_x",ncol=1) +
  geom_boxplot(aes(y=value,fill=effsize_cut)) +
  geom_hline(linetype = "dashed",yintercept=0, size = 0.5, color = "gray28") +
  scale_fill_manual(values=manual_color,name="Effect size",limits = force) +
  labs(y = "Detection lag differences")+
  scale_x_discrete(name = "Memory sizes",labels=c("27"=expression(bold("27")), parse=TRUE))+
  theme_bw()+  
  geom_text(aes(y=max(value,na.rm=TRUE),label = p.adj.signif), position = position_dodge(width=0.9),vjust=-0.7)+
  scale_y_continuous(expand = expansion(mult = c(0.2, 0.5)))

sp_memory <- df_diff %>%
  pivot_longer(c(diff_recall_soft, diff_precision_soft, diff_F1_soft)) %>%
  group_by(col) %>%
  filter(name==paste0("diff_",col)) %>%
  ggplot(aes(x = factor(mem,levels=c(unique(mem))))) +
  geom_hline(linetype = "dashed",yintercept=0, size = 0.5, color = "gray28") +
  geom_path(aes(y=-effsize.r_rank_biserial, group = 1),color="grey",size=1) +
  geom_point(aes(y=-effsize.r_rank_biserial, fill=effsize_cut, color=effsize_cut, size=effsize_cut),shape=21) +
  facet_wrap(~factor(name, levels=c('diff_recall_soft', 'diff_precision_soft', 'diff_F1_soft'),labels=par_labels_metrics), scales = "free_x",ncol=1) +
  scale_size_manual(values = manual_size,name="Effect size",limits = force) +
  scale_color_manual(values = manual_color,name="Effect size",limits = force) +
  scale_fill_manual(values=manual_color,name="Effect size",limits = force) +
  scale_y_continuous(position = "right", name = "Effect size",expand = expansion(mult = c(0.2, 0.5)))+ 
  scale_x_discrete(name = "Memory sizes",labels=c("27"=expression(bold("27")), parse=TRUE))+
  theme_bw()+  
  geom_text(aes(y=max(-effsize.r_rank_biserial),label = p.adj.signif), position = position_dodge(width=0.9),vjust=-0.7)

ggarrange(bp_memory, sp_memory,common.legend = FALSE, legend = "bottom")
ggsave("sensibility_memory_eval.pdf", width = 25, height = 15, units = "cm")
