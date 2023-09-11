library(ggplot2)
library(tidyverse)
library(rstatix)
#install.packages("rcompanion") # for wilcoxonPairedRC
library(rcompanion)

# Gathering of results =========================================================

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


# OR

#Concatenate results from all series from result dataframes.....................
exp_res_df_0F <- readRDS(paste0(getwd(),"/dev/plots/exp_res_df_0F.rds"))
exp_res_df_AE <- readRDS(paste0(getwd(),"/dev/plots/exp_res_df_AE.rds"))
exp_res_df <- rbind(exp_res_df_0F,exp_res_df_AE)
saveRDS(exp_res_df,file=paste0(getwd(),"/dev/plots/exp_res_df.rds"))

exp_res_df$curve <- factor(exp_res_df$curve,levels=curves)
exp_res_df$var <- factor(exp_res_df$var, ordered = TRUE,levels=names(volatility))
exp_res_df$exp <- factor(exp_res_df$exp,levels=c("baseline","A","B","C","D","E","F"))

#Concatenate results for D and baseline for memory sensibility analysis
exp_res_df_0D_mem <- readRDS(paste0(getwd(),"/dev/plots/exp_res_df_0D_mem.rds"))



# Best results =================================================================

# Best lags --------------------------------------------------------------------
best_res_df <- exp_res_df %>%
  filter(mem ==27) %>%
  group_by(curve,var,seed,exp) %>%
  filter(detlag_batch == min(detlag_batch, na.rm = TRUE)) %>%
  slice(1) %>%
  ungroup()
best_res_df$curve <- factor(best_res_df$curve,levels=curves)
best_res_df$var <- factor(best_res_df$var, ordered = TRUE,levels=names(volatility))
best_res_df$exp <- factor(best_res_df$exp,levels=c("baseline","A","B","C","D","E","F"))

saveRDS(best_res_df,file=paste0(getwd(),"/dev/plots/best_res_df.rds"))

# Detection plots
best_res_plots <- lapply(1:nrow(best_res_df), function(r) {
  b <- best_res_df[r,]
  data_name <- paste0(b$curve,"_",b$var,"_",b$seed)
  det <- readRDS(paste0(getwd(),"/dev/plots/",data_name,"/","det_res_exp_",data_name,".rds"))
  det[[as.character(b$exp)]]$detection[[b$id]]$detection$plot
})
saveRDS(best_res_plots,file=paste0(getwd(),"/dev/plots/best_res_plots.rds"))


# Best lags by memory size -----------------------------------------------------
#Best results by memory size
best_res_df_mem <- exp_res_df_0D_mem %>%
  group_by(curve,var,seed,exp,mem) %>%
  filter(detlag_batch == min(detlag_batch, na.rm = TRUE)) %>%
  slice(1) %>%
  ungroup()
best_res_df$curve <- factor(best_res_df$curve,levels=curves)
best_res_df$var <- factor(best_res_df$var, ordered = TRUE,levels=names(volatility))
best_res_df$exp <- factor(best_res_df$exp,levels=c("baseline","A","B","C","D","E","F"))

saveRDS(best_res_df_mem,file=paste0(getwd(),"/dev/plots/best_res_df_mem.rds"))


# Detection plots
best_res_plots_mem <- lapply(1:nrow(best_res_df_mem), function(r) {
  b <- best_res_df_mem[r,]
  data_name <- paste0(b$curve,"_",b$var,"_",b$seed)
  det <- readRDS(paste0(getwd(),"/dev/plots/",data_name,"/","det_res_exp_",data_name,".rds"))
  det[[as.character(b$exp)]]$detection[[b$id]]$detection$plot
})
saveRDS(best_res_plots_mem,file=paste0(getwd(),"/dev/plots/best_res_plots.rds"))




# Statistical tests ============================================================

# baseline versus D ............................................................

#General results test
baseline <- best_res_df %>%
  filter(exp %in% c("baseline"))
expD <- best_res_df %>%
  filter(exp %in% c("D"))
wilcox.test(baseline$detlag_batch, expD$detlag_batch, paired = TRUE, alternative = "greater")

# Tests by volatility
best_res_df  %>%
  filter(exp %in% c("baseline","D")) %>%
  group_by(var) %>%
  wilcox_test(detlag_batch ~ exp, paired = TRUE, alternative = "greater")  %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")

# Tests by volatility with Wilcoxon effect size
best_res_df  %>%
  filter(exp %in% c("baseline","D")) %>%
  group_by(var) %>%
  wilcox_effsize(detlag_batch ~ exp, paired = TRUE, alternative = "greater")


# baseline versus ablation ............................................................

#source function from rstatix
get_wilcox_effsize_magnitude <- function(d){
  magnitude.levels = c(0.3, 0.5, Inf)
  magnitude = c("small","moderate","large")
  d.index <- findInterval(abs(d), magnitude.levels)+1
  magnitude <- factor(magnitude[d.index], levels = magnitude, ordered = TRUE)
  magnitude
}

ablation_tests <- list()
for(a in c("A","B","C","D")){
  ablation_tests_signif <- best_res_df  %>%
    filter(exp %in% c("baseline",a)) %>%
    group_by(curve,var,seed) %>%
    filter(n() > 1) %>%
    group_by(var) %>%
    wilcox_test(detlag_batch ~ exp, paired = TRUE, alternative = "greater")  %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")
  
  #wilcoxonPairedRC: matched-pairs rank biserial correlation coefficient
  #(gives negative effect sizes when appropriate)
  ablation_tests_effsize <- best_res_df  %>%
    filter(exp %in% c("baseline",a)) %>%
    group_by(curve,var,seed) %>%
    filter(n() > 1) %>%
    group_by(var) %>%
    group_map(~data.frame(var=unique(.x$var), 
                          effsize=wilcoxonPairedRC(.x$detlag_batch, factor(.x$exp,levels=c("baseline",a)), zero.method="Wilcoxon")
    ),.keep=TRUE) %>%
    bind_rows() %>%
    mutate(magnitude = get_wilcox_effsize_magnitude(effsize))
  
  ablation_tests[[a]] <- merge(ablation_tests_signif,ablation_tests_effsize, sort = F)
}

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