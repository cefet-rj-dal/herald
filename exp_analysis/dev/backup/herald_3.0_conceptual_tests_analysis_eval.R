library(ggplot2)
library(tidyverse)
library(rstatix)
#install.packages("rcompanion") # for wilcoxonPairedRC
library(rcompanion)
#install.packages("effectsize") # for rank_biserial
library(effectsize)
#https://easystats.github.io/effectsize/reference/interpret_r.html
#https://easystats.github.io/effectsize/articles/interpret.html
#https://easystats.github.io/effectsize/reference/rank_biserial.html

# Gathering of results =========================================================

curves <- c("linear","log","sqrt","poly2","poly3","exp","sigmoid")
volatility <- c(low=0.025,medium=0.05,high=0.1) #volatility levels
seeds <- as.character(seq(from=123,by=1,length.out=30))

# Concatenate results from all series folders...................................
exp_res_df_lags <- data.frame()
for(curve in curves){
  for(var in names(volatility)){
    for(seed in seeds){
      data_name <- paste0(curve,"_",var,"_",seed)
      if(data_name=="linear_medium_146") stop("Let's break out!")
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
      #if(data_name=="linear_medium_146") stop("Let's break out!")
      df <- readRDS(paste0(getwd(),"/dev/plots/",data_name,"/","res_exp_",data_name,"_eval.rds"))
      exp_res_df_eval <- rbind(exp_res_df_eval, cbind(curve=curve,var=var,seed=seed,df))
    }
  }
}
saveRDS(exp_res_df_eval,file=paste0(getwd(),"/dev/plots/exp_res_df_eval.rds"))


# OR

#Concatenate results from all series from result dataframes.....................
#Concatenate results for D and baseline for memory sensibility analysis
exp_res_df_lags <- readRDS(paste0(getwd(),"/dev/plots/exp_res_df_0D_mem_lags.rds"))
exp_res_df_eval <- readRDS(paste0(getwd(),"/dev/plots/exp_res_df_0D_mem_eval.rds"))



# Best results =================================================================

# Best lags --------------------------------------------------------------------
best_res_df_lags <- exp_res_df_lags %>%
  filter(mem ==27) %>%
  group_by(curve,var,seed,exp,event_idx) %>%
  filter(detlag_batch == min(detlag_batch, na.rm = TRUE)) %>%
  slice(1) %>%
  ungroup()
best_res_df_lags$curve <- factor(best_res_df_lags$curve,levels=curves)
best_res_df_lags$var <- factor(best_res_df_lags$var, ordered = TRUE,levels=names(volatility))
best_res_df_lags$exp <- factor(best_res_df_lags$exp,levels=c("baseline","A","B","C","D"))

saveRDS(best_res_df_lags,file=paste0(getwd(),"/dev/plots/best_res_df_lags.rds"))

# Detection plots
best_res_plots_lags <- lapply(1:nrow(best_res_df_lags), function(r) {
  b <- best_res_df_lags[r,]
  data_name <- paste0(b$curve,"_",b$var,"_",b$seed)
  det <- readRDS(paste0(getwd(),"/dev/plots/",data_name,"/","det_res_exp_",data_name,".rds"))
  det[[as.character(b$exp)]]$detection[[b$id]]$detection$plot
})
saveRDS(best_res_plots_lags,file=paste0(getwd(),"/dev/plots/best_res_plots_lags.rds"))


#Best results by memory size
best_res_df_lags_mem <- exp_res_df_lags %>%
  group_by(curve,var,seed,exp,mem,event_idx) %>%
  filter(detlag_batch == min(detlag_batch, na.rm = TRUE)) %>%
  slice(1) %>%
  ungroup()
best_res_df_lags_mem$curve <- factor(best_res_df_lags_mem$curve,levels=curves)
best_res_df_lags_mem$var <- factor(best_res_df_lags_mem$var, ordered = TRUE,levels=names(volatility))
best_res_df_lags_mem$exp <- factor(best_res_df_lags_mem$exp,levels=c("baseline","A","B","C","D"))

saveRDS(best_res_df_lags_mem,file=paste0(getwd(),"/dev/plots/best_res_df_lags_mem.rds"))


# Detection plots
best_res_plots_lags_mem <- lapply(1:nrow(best_res_df_lags_mem), function(r) {
  b <- best_res_df_lags_mem[r,]
  data_name <- paste0(b$curve,"_",b$var,"_",b$seed)
  det <- readRDS(paste0(getwd(),"/dev/plots/",data_name,"/","det_res_exp_",data_name,".rds"))
  det[[as.character(b$exp)]]$detection[[b$id]]$detection$plot
})
saveRDS(best_res_plots_lags_mem,file=paste0(getwd(),"/dev/plots/best_res_plots_lags_mem.rds"))


# Best eval --------------------------------------------------------------------
best_res_df_eval <- exp_res_df_eval %>%
  filter(mem ==27) %>%
  group_by(curve,var,seed,exp) %>%
  filter(F1_soft == min(F1_soft, na.rm = TRUE)) %>%
  slice(1) %>%
  ungroup()
best_res_df_eval$curve <- factor(best_res_df_eval$curve,levels=curves)
best_res_df_eval$var <- factor(best_res_df_eval$var, ordered = TRUE,levels=names(volatility))
best_res_df_eval$exp <- factor(best_res_df_eval$exp,levels=c("baseline","A","B","C","D"))

saveRDS(best_res_df_eval,file=paste0(getwd(),"/dev/plots/best_res_df_eval.rds"))

# Detection plots
best_res_plots_eval <- lapply(1:nrow(best_res_df_eval), function(r) {
  b <- best_res_df_eval[r,]
  data_name <- paste0(b$curve,"_",b$var,"_",b$seed)
  det <- readRDS(paste0(getwd(),"/dev/plots/",data_name,"/","det_res_exp_",data_name,".rds"))
  det[[as.character(b$exp)]]$detection[[b$id]]$detection$plot
})
saveRDS(best_res_plots_eval,file=paste0(getwd(),"/dev/plots/best_res_plots_eval.rds"))


#Best results by memory size
best_res_df_eval_mem <- exp_res_df_eval %>%
  group_by(curve,var,seed,exp,mem) %>%
  filter(F1_soft == min(F1_soft, na.rm = TRUE)) %>%
  slice(1) %>%
  ungroup()
best_res_df_eval_mem$curve <- factor(best_res_df_eval_mem$curve,levels=curves)
best_res_df_eval_mem$var <- factor(best_res_df_eval_mem$var, ordered = TRUE,levels=names(volatility))
best_res_df_eval_mem$exp <- factor(best_res_df_eval_mem$exp,levels=c("baseline","A","B","C","D"))

saveRDS(best_res_df_eval_mem,file=paste0(getwd(),"/dev/plots/best_res_df_eval_mem.rds"))


# Detection plots
best_res_plots_eval_mem <- lapply(1:nrow(best_res_df_eval_mem), function(r) {
  b <- best_res_df_eval_mem[r,]
  data_name <- paste0(b$curve,"_",b$var,"_",b$seed)
  det <- readRDS(paste0(getwd(),"/dev/plots/",data_name,"/","det_res_exp_",data_name,".rds"))
  det[[as.character(b$exp)]]$detection[[b$id]]$detection$plot
})
saveRDS(best_res_plots_eval_mem,file=paste0(getwd(),"/dev/plots/best_res_plots_eval_mem.rds"))



# Statistical tests ============================================================

# Multiple lags ----------------------------------------------------------------

# Tests by volatility
best_res_df_lags  %>%
  filter(exp %in% c("baseline","D")) %>%
  group_by(curve,var,seed,event_idx) %>%
  filter(n() > 1) %>%
  group_by(var) %>%
  wilcox_test(detlag_batch ~ exp, paired = TRUE, alternative = "greater")  %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")

# Tests by volatility with Wilcoxon effect size
best_res_df_lags  %>%
  filter(exp %in% c("baseline","D")) %>%
  group_by(curve,var,seed,event_idx) %>%
  filter(n() > 1) %>%
  group_by(var) %>%
  wilcox_effsize(detlag_batch ~ exp, paired = TRUE, alternative = "greater")


# Tests by volatility and events
best_res_df_lags  %>%
  filter(exp %in% c("baseline","D")) %>%
  group_by(curve,var,seed,event_idx) %>%
  filter(n() > 1) %>%
  group_by(var,event_idx) %>%
  wilcox_test(detlag_batch ~ exp, paired = TRUE, alternative = "greater")  %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")

# Tests by volatility with Wilcoxon effect size
best_res_df_lags  %>%
  filter(exp %in% c("baseline","D")) %>%
  group_by(curve,var,seed,event_idx) %>%
  filter(n() > 1) %>%
  group_by(var,event_idx) %>%
  wilcox_effsize(detlag_batch ~ exp, paired = TRUE, alternative = "greater")


# Multiple lags (obs) -------------------------------------------------------------

# Tests by volatility
best_res_df_lags  %>%
  filter(exp %in% c("baseline","D")) %>%
  group_by(curve,var,seed,event_idx) %>%
  filter(n() > 1) %>%
  group_by(var) %>%
  wilcox_test(detlag_obs ~ exp, paired = TRUE, alternative = "greater")  %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")

# Tests by volatility with Wilcoxon effect size
best_res_df_lags  %>%
  filter(exp %in% c("baseline","D")) %>%
  group_by(curve,var,seed,event_idx) %>%
  filter(n() > 1) %>%
  group_by(var) %>%
  wilcox_effsize(detlag_obs ~ exp, paired = TRUE, alternative = "greater")


# Tests by volatility and events
best_res_df_lags  %>%
  filter(exp %in% c("baseline","D")) %>%
  group_by(curve,var,seed,event_idx) %>%
  filter(n() > 1) %>%
  group_by(var,event_idx) %>%
  wilcox_test(detlag_obs ~ exp, paired = TRUE, alternative = "greater")  %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")


# Tests by volatility with Wilcoxon effect size
best_res_df_lags  %>%
  filter(exp %in% c("baseline","D")) %>%
  group_by(curve,var,seed,event_idx) %>%
  filter(n() > 1) %>%
  group_by(var,event_idx) %>%
  wilcox_effsize(detlag_obs ~ exp, paired = TRUE, alternative = "greater")



# SoftED evaluation ------------------------------------------------------------

# Tests by volatility
best_res_df_eval  %>%
  filter(exp %in% c("baseline","D")) %>%
  group_by(curve,var,seed) %>%
  filter(n() > 1) %>%
  group_by(var) %>%
  wilcox_test(F1_soft ~ exp, paired = TRUE, alternative = "greater")  %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")


# Tests by volatility with Wilcoxon effect size
best_res_df_lags  %>%
  filter(exp %in% c("baseline","D")) %>%
  group_by(curve,var,seed,event_idx) %>%
  filter(n() > 1) %>%
  group_by(var) %>%
  wilcox_effsize(detlag_batch ~ exp, paired = TRUE, alternative = "greater")


# Tests by volatility and events
best_res_df_lags  %>%
  filter(exp %in% c("baseline","D")) %>%
  group_by(curve,var,seed,event_idx) %>%
  filter(n() > 1) %>%
  group_by(var,event_idx) %>%
  wilcox_test(detlag_batch ~ exp, paired = TRUE, alternative = "greater")  %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")

# Tests by volatility with Wilcoxon effect size
best_res_df_lags  %>%
  filter(exp %in% c("baseline","D")) %>%
  group_by(curve,var,seed,event_idx) %>%
  filter(n() > 1) %>%
  group_by(var,event_idx) %>%
  wilcox_effsize(detlag_batch ~ exp, paired = TRUE, alternative = "greater")


# baseline versus ablation ............................................................

ablation_tests_signif <- best_res_df  %>%
  filter(exp %in% c("baseline","D")) %>%
  group_by(curve,var,seed) %>%
  filter(n() > 1) %>%
  group_by(var) %>%
  wilcox_test(detlag_batch ~ exp, paired = TRUE, alternative = "greater")  %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")

#wilcoxonPairedRC: matched-pairs rank biserial correlation coefficient
#(gives negative effect sizes when appropriate)
ablation_tests_effsize <- best_res_df  %>%
  filter(exp %in% c("baseline","D")) %>%
  group_by(curve,var,seed) %>%
  filter(n() > 1) %>%
  group_by(var) %>%
  group_map(~data.frame(var=unique(.x$var), 
                        effsize=interpret(rank_biserial(.x$detlag_batch, factor(.x$exp,levels=c("baseline",a)), paired = TRUE, alternative = "greater"), rules = "funder2019")
  ),.keep=TRUE) %>%
  bind_rows()
  
ablation_tests <- merge(ablation_tests_signif,ablation_tests_effsize, sort = F)

ablation_tests



# TODO
# Plots ------------------------------------------------------------------------
# Create the ggplot with faceted boxplots
exp_res_df %>%
  filter(exp %in% c("A","B","C","D","E","F")) %>%
  ggplot(aes(x = exp, y = detlag_batch)) +
  geom_boxplot(fill="deepskyblue") +
  stat_summary(fun = "mean", geom = "point", shape = 17, size = 3, color = "red") +
  stat_summary(fun = "mean", geom = "line", aes(group = 1), size = 1, color = "red") +
  facet_wrap(~ factor(var, levels=c('low', 'medium', 'high')), scales = "free_x", ncol = 1) +
  labs(title = "Boxplots of detection lags (batch) grouped by data volatility",
       x = "Herald settings",
       y = "Detection lag")+
  theme_bw()

# Create the ggplot with faceted boxplots
exp_res_df %>%
  filter(exp %in% c("A","B","C","D","E","F")) %>%
  ggplot(aes(x = exp, y = detlag_batch)) +
  geom_boxplot(fill="deepskyblue") +
  stat_summary(fun = "mean", geom = "point", shape = 17, size = 3, color = "red") +
  stat_summary(fun = "mean", geom = "line", aes(group = 1), size = 1, color = "red") +
  facet_wrap(~ factor(var, levels=c('low', 'medium', 'high')), scales = "free_x", ncol = 1) +
  labs(title = "Boxplots of detection lags (obs) grouped by data volatility",
       x = "Herald settings",
       y = "Detection lag")+
  theme_bw()


# Create the ggplot with faceted barplots
best_res_df %>%
  filter(exp %in% c("A","B","C","D","E","F")) %>%
  group_by(curve,var,seed) %>%
  arrange(detlag_batch) %>%
  slice_head(n = 3) %>%
  group_by(var, exp) %>%
  summarise(frequency = n()) %>%
  mutate(perc = frequency/sum(frequency)) %>%
  ggplot(aes(x = exp, y = perc)) +
  geom_bar(stat = "identity",fill="deepskyblue") +
  facet_wrap(~ factor(var, levels=c('low', 'medium', 'high')), ncol = 1) +
  labs(title = "Frequency of settings among the top 3 detection lags (batch)",
       x = "Herald settings",
       y = "Frequency")+
  theme_bw()

best_res_df %>%
  filter(exp %in% c("A","B","C","D","E","F")) %>%
  group_by(curve,var,seed) %>%
  arrange(detlag_batch) %>%
  slice_head(n = 3) %>%
  group_by(var, exp) %>%
  summarise(frequency = n()) %>%
  mutate(perc = frequency/sum(frequency)) %>%
  ggplot(aes(x = exp, y = perc)) +
  geom_bar(stat = "identity",fill="deepskyblue") +
  facet_wrap(~ factor(var, levels=c('low', 'medium', 'high')), ncol = 1) +
  labs(title = "Frequency of settings among the top 3 detection lags (obs)",
       x = "Herald settings",
       y = "Frequency")+
  theme_bw()

# Create the ggplot with faceted barplots
best_res_df %>%
  filter(exp %in% c("A","B","C","D","E","F")) %>%
  group_by(curve,var,seed) %>%
  arrange(detlag_batch) %>%
  slice_head(n = 1) %>%
  group_by(var, exp) %>%
  summarise(frequency = n()) %>%
  mutate(perc = frequency/sum(frequency)) %>%
  ggplot(aes(x = exp, y = perc)) +
  geom_bar(stat = "identity",fill="deepskyblue") +
  facet_wrap(~ factor(var, levels=c('low', 'medium', 'high')), ncol = 1) +
  labs(title = "Frequency of settings among the top 1 detection lags (batch)",
       x = "Herald settings",
       y = "Frequency")+
  theme_bw()

best_res_df %>%
  filter(exp %in% c("A","B","C","D","E","F")) %>%
  group_by(curve,var,seed) %>%
  arrange(detlag_obs) %>%
  slice_head(n = 1) %>%
  group_by(var, exp) %>%
  summarise(frequency = n()) %>%
  mutate(perc = frequency/sum(frequency)) %>%
  ggplot(aes(x = exp, y = perc)) +
  geom_bar(stat = "identity",fill="deepskyblue") +
  facet_wrap(~ factor(var, levels=c('low', 'medium', 'high')), ncol = 1) +
  labs(title = "Frequency of settings among the top 1 detection lags (obs)",
       x = "Herald settings",
       y = "Frequency")+
  theme_bw()


# Create the ggplot with faceted barplots
best_res_df %>%
  filter(exp %in% c("baseline","F")) %>%
  group_by(curve,var,seed) %>%
  arrange(detlag_batch) %>%
  slice_head(n = 1) %>%
  group_by(var, exp) %>%
  summarise(frequency = n()) %>%
  mutate(perc = frequency/sum(frequency)) %>%
  ggplot(aes(x = exp, y = perc)) +
  geom_bar(stat = "identity",fill="deepskyblue") +
  facet_wrap(~ factor(var, levels=c('low', 'medium', 'high')), ncol = 1) +
  labs(title = "Frequency of settings among the top 1 detection lags (batch)",
       x = "Herald settings",
       y = "Frequency")+
  theme_bw()

best_res_df %>%
  filter(exp %in% c("baseline","F")) %>%
  group_by(curve,var,seed) %>%
  arrange(detlag_obs) %>%
  slice_head(n = 1) %>%
  group_by(var, exp) %>%
  summarise(frequency = n()) %>%
  mutate(perc = frequency/sum(frequency)) %>%
  ggplot(aes(x = exp, y = perc)) +
  geom_bar(stat = "identity",fill="deepskyblue") +
  facet_wrap(~ factor(var, levels=c('low', 'medium', 'high')), ncol = 1) +
  labs(title = "Frequency of settings among the top 1 detection lags (obs)",
       x = "Herald settings",
       y = "Frequency")+
  theme_bw()


# Create the ggplot with faceted barplots
best_res_df %>%
  filter(exp %in% c("baseline", "F")) %>%
  group_by(curve, var, seed) %>%
  summarise(diff_detlag = detlag_batch[exp == "baseline"] - detlag_batch[exp == "F"])  %>%
  mutate(comparison = ifelse(diff_detlag > 0, "F < Baseline", ifelse(diff_detlag == 0, "F = Baseline", "F > Baseline"))) %>%
  group_by(var, comparison) %>%
  summarise(frequency = n()) %>%
  mutate(perc = frequency / sum(frequency)) %>%
  ggplot(aes(x = factor(var, levels=c('low', 'medium', 'high')), y = perc, fill = comparison)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Comparison of detection lags (batch) between F and Baseline",
       x = "Volatility",
       y = "Percentage") +
  scale_fill_manual(values = c("lightgreen", "deepskyblue", "tomato"), guide = guide_legend(title = "Lag comparison")) +
  theme_bw()


best_res_df %>%
  filter(exp %in% c("baseline", "F")) %>%
  group_by(curve, var, seed) %>%
  summarise(diff_detlag = detlag_obs[exp == "baseline"] - detlag_obs[exp == "F"])  %>%
  mutate(comparison = ifelse(diff_detlag > 0, "F < Baseline", ifelse(diff_detlag == 0, "F = Baseline", "F > Baseline"))) %>%
  group_by(var, comparison) %>%
  summarise(frequency = n()) %>%
  mutate(perc = frequency / sum(frequency)) %>%
  ggplot(aes(x = factor(var, levels=c('low', 'medium', 'high')), y = perc, fill = comparison)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Comparison of detection lags (obs) between F and Baseline",
       x = "Volatility",
       y = "Percentage") +
  scale_fill_manual(values = c("lightgreen", "deepskyblue", "tomato"), guide = guide_legend(title = "Lag comparison")) +
  theme_bw()