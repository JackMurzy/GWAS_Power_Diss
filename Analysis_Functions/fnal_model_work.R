# Complete workflow for predicting AUC which we want to be maximal
library(caret)
library(car)
library(ggplot2)
library(ggsignif)
library(ggnewscale)
library(tidyverse)
library(gridExtra)
library(grid)
library(cowplot)
library(ggrepel)

directory <- "~/Desktop/University/Biomedicine/Y3/Dissertation/PLINK_Tutorials/Simulate_Data_Test"
setwd(directory)
# Get plots for the individual results

# Read in all the data
nindv_data<- Indv_Rdata$summary_data %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  mutate(lcAUC = ifelse(lcAUC < 0.5,0.5,lcAUC),
         ucAUC = ifelse(ucAUC > 1,1,ucAUC),
         lcbFstat = ifelse(lcbFstat < 0,0,lcbFstat),
         ucbFstat = ifelse(ucbFstat > 1,1,ucbFstat),
         lcfFstat = ifelse(lcfFstat < 0,0,lcfFstat),
         ucfFstat = ifelse(ucfFstat > 1,1,ucfFstat))
nSNP_data <- nSNP_Rdata$summary_data %>%
  filter(n_snps != 5000) %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  mutate(lcAUC = ifelse(lcAUC < 0.5,0.5,lcAUC),
         ucAUC = ifelse(ucAUC > 1,1,ucAUC),
         lcbFstat = ifelse(lcbFstat < 0,0,lcbFstat),
         ucbFstat = ifelse(ucbFstat > 1,1,ucbFstat),
         lcfFstat = ifelse(lcfFstat < 0,0,lcfFstat),
         ucfFstat = ifelse(ucfFstat > 1,1,ucfFstat))
causal_snp_data <- ncSNP_Rdata$summary_data %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  mutate(lcAUC = ifelse(lc_AUC_val < 0.5,0.5,lc_AUC_val),
         ucAUC = ifelse(uc_AUC_val > 1,1,uc_AUC_val),
         lcbFstat = ifelse(lc_Bonf_Fval < 0,0,lc_Bonf_Fval),
         ucbFstat = ifelse(uc_Bonf_Fval > 1,1,uc_Bonf_Fval),
         lcfFstat = ifelse(lc_FDR_Fval < 0,0,lc_FDR_Fval),
         ucfFstat = ifelse(uc_FDR_Fval > 1,1,uc_FDR_Fval))
herit_data <- sHerit_Rdata$summary_data %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  mutate(lcAUC = ifelse(lc_AUC_val < 0.5,0.5,lc_AUC_val),
         ucAUC = ifelse(uc_AUC_val > 1,1,uc_AUC_val),
         lcbFstat = ifelse(lc_Bonf_Fval < 0,0,lc_Bonf_Fval),
         ucbFstat = ifelse(uc_Bonf_Fval > 1,1,uc_Bonf_Fval),
         lcfFstat = ifelse(lc_FDR_Fval < 0,0,lc_FDR_Fval),
         ucfFstat = ifelse(uc_FDR_Fval > 1,1,uc_FDR_Fval))
cMAF_data <- cMAF_Rdata$summary_data %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  mutate(lcAUC = ifelse(lcAUC < 0.5,0.5,lcAUC),
         ucAUC = ifelse(ucAUC > 1,1,ucAUC),
         lcbFstat = ifelse(lcbFstat < 0,0,lcbFstat),
         ucbFstat = ifelse(ucbFstat > 1,1,ucbFstat),
         lcfFstat = ifelse(lcfFstat < 0,0,lcfFstat),
         ucfFstat = ifelse(ucfFstat > 1,1,ucfFstat))
cascntrl_data <- CCRatio_Rdata$summary_data %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  mutate(lcAUC = ifelse(lc_AUC_val < 0.5,0.5,lc_AUC_val),
         ucAUC = ifelse(uc_AUC_val > 1,1,uc_AUC_val),
         lcbFstat = ifelse(lc_Bonf_Fval < 0,0,lc_Bonf_Fval),
         ucbFstat = ifelse(uc_Bonf_Fval > 1,1,uc_Bonf_Fval),
         lcfFstat = ifelse(lc_FDR_Fval < 0,0,lc_FDR_Fval),
         ucfFstat = ifelse(uc_FDR_Fval > 1,1,uc_FDR_Fval))


# Get min and max AUC and Fvals
all_auc_fvals <- nindv_data %>%
  select(lcbFstat,ucbFstat,lcfFstat,ucfFstat,lcAUC,ucAUC) %>%
  mutate(class = "indv") %>%
  rbind(cMAF_data %>%
          select(lcbFstat,ucbFstat,lcfFstat,ucfFstat,lcAUC,ucAUC) %>%
          mutate(class = "cmaf")) %>%
  rbind(nSNP_data %>%
          select(lcbFstat,ucbFstat,lcfFstat,ucfFstat,lcAUC,ucAUC) %>%
          mutate(class = "snp")) %>%
  rbind(cascntrl_data %>%
          select(lcbFstat,ucbFstat,lcfFstat,ucfFstat,lcAUC,ucAUC) %>%
          mutate(class = "ccrat")) %>%
  rbind(herit_data %>%
          select(lcbFstat,ucbFstat,lcfFstat,ucfFstat,lcAUC,ucAUC) %>%
          mutate(class = "herit")) %>%
  rbind(causal_snp_data %>%
          select(lcbFstat,ucbFstat,lcfFstat,ucfFstat,lcAUC,ucAUC) %>%
          mutate(class = "csnp")) %>%
  mutate(lcFstat = ifelse(lcbFstat < lcfFstat,lcbFstat,lcfFstat),
         ucFstat = ifelse(ucbFstat < ucfFstat,ucbFstat,ucfFstat)) 
total_summarised_auc_fval <- all_auc_fvals %>%
  dplyr::summarise(across(everything(), list(min = min, max = max), na.rm = TRUE)) %>%
  select(lcAUC_min,ucAUC_max,lcFstat_min,ucFstat_max) %>%
  mutate(lcAUC_min = ifelse(lcAUC_min < 0,0,lcAUC_min),
         ucAUC_max = ifelse(ucAUC_max > 1,1,ucAUC_max),
         lcFstat_min = ifelse(lcFstat_min < 0,0,lcFstat_min),
         ucFstat_max = ifelse(ucFstat_max > 1,1,ucFstat_max),
         lcAUC_min_round = floor(lcAUC_min*10)/10,
         ucAUC_max_round = ceiling(ucAUC_max*10)/10,
         lcFstat_min_round = floor(lcFstat_min*10)/10,
         ucFstat_max_round = ceiling(ucFstat_max*10)/10)

min_auc <- total_summarised_auc_fval$lcAUC_min_round
max_auc <- total_summarised_auc_fval$ucAUC_max_round
min_fval <- total_summarised_auc_fval$lcFstat_min_round
max_fval <- total_summarised_auc_fval$ucFstat_max_round

# Individual data
indv_auc_plot <- ggplot(nindv_data, aes(x = n_individuals)) +
  geom_line(aes(y = mAUC, color = "AUC")) +
  geom_point(aes(y = mAUC, color = "AUC")) +
  geom_ribbon(aes(ymin = lcAUC, ymax = ucAUC, fill = "AUC"), alpha = 0.2) +
  scale_color_manual(name = str_wrap("Mean AUC ± 95% CI", width = 9), values = c("AUC" = "blue")) +
  scale_fill_manual(name = str_wrap("Mean AUC ± 95% CI", width = 9), values = c("AUC" = "blue")) +
  labs(title = "AUC vs. nIndividuals", x = "nIndividuals", y = "AUC") +
  theme_minimal() +
  theme(axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.5,"in")
  ) + 
  ylim(min_auc, max_auc)

indv_fstat_plot <- ggplot(nindv_data, aes(x = n_individuals)) +
  geom_line(aes(y = mbFstat, color = "Bonferroni")) +
  geom_point(aes(y = mbFstat, color = "Bonferroni")) +
  geom_ribbon(aes(ymin = lcbFstat, ymax = ucbFstat, fill = "Bonferroni"), alpha = 0.2) +
  geom_line(aes(y = mfFstat, color = "FDR")) +
  geom_point(aes(y = mfFstat, color = "FDR")) +
  geom_ribbon(aes(ymin = lcfFstat, ymax = ucfFstat, fill = "FDR"), alpha = 0.2) +
  scale_color_manual(name = str_wrap("Mean F-value ± 95% CI", width = 12), values = c("Bonferroni" = "blue", "FDR" = "red")) +
  scale_fill_manual(name = str_wrap("Mean F-value ± 95% CI", width = 12), values = c("Bonferroni" = "blue", "FDR" = "red")) +
  labs(title = "F-value at Bonferroni and FDR thresholds",
       x = "nIndividuals",
       y = "F-value") +
  theme_minimal() + 
  theme(axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.5,"in")
  ) + 
  ylim(min_fval, max_fval)


# Number of SNPs
snp_auc_plot <- ggplot(nSNP_data, aes(x = n_snps)) +
  geom_line(aes(y = mAUC, color = "AUC")) +
  geom_point(aes(y = mAUC, color = "AUC")) +
  geom_ribbon(aes(ymin = lcAUC, ymax = ucAUC, fill = "AUC"), alpha = 0.2) +
  scale_color_manual(name = "Mean AUC ± 95% CI", values = c("AUC" = "blue")) +
  scale_fill_manual(name = "Mean AUC ± 95% CI", values = c("AUC" = "blue")) +
  labs(title = "AUC vs. nSNPs", x = "nSNPs", y = "AUC") +
  theme_minimal() +
  theme(axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.5,"in")
  ) + 
  ylim(min_auc, max_auc)
snp_fstat_plot <- ggplot(nSNP_data, aes(x = n_snps)) +
  geom_line(aes(y = mbFstat, color = "Bonferroni")) +
  geom_point(aes(y = mbFstat, color = "Bonferroni")) +
  geom_ribbon(aes(ymin = lcbFstat, ymax = ucbFstat, fill = "Bonferroni"), alpha = 0.2) +
  geom_line(aes(y = mfFstat, color = "FDR")) +
  geom_point(aes(y = mfFstat, color = "FDR")) +
  geom_ribbon(aes(ymin = lcfFstat, ymax = ucfFstat, fill = "FDR"), alpha = 0.2) +
  scale_color_manual(name = "Mean F-value ± 95% CI", values = c("Bonferroni" = "blue", "FDR" = "red")) +
  scale_fill_manual(name = "Mean F-value ± 95% CI", values = c("Bonferroni" = "blue", "FDR" = "red")) +
  labs(title = "F-value at Bonferroni and FDR thresholds",
       x = "nSNPs",
       y = "F-value") +
  theme_minimal() + 
  theme(axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.5,"in")
  ) + 
  ylim(min_fval, max_fval)


# Number of Causal SNPs
causal_snp_auc_plot <- ggplot(causal_snp_data, aes(x = ncSNP)) +
  geom_line(aes(y = m_AUC_val, color = "AUC")) +
  geom_point(aes(y = m_AUC_val, color = "AUC")) +
  geom_ribbon(aes(ymin = lcAUC, ymax = ucAUC, fill = "AUC"), alpha = 0.2) +
  scale_color_manual(name = "Mean AUC ± 95% CI", values = c("AUC" = "blue")) +
  scale_fill_manual(name = "Mean AUC ± 95% CI", values = c("AUC" = "blue")) +
  labs(title = "AUC vs. nCausal SNPs", x = "nCausal SNPs", y = "AUC") +
  theme_minimal() +
  theme(axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.5,"in")
  ) + 
  ylim(min_auc, max_auc)
causal_snp_fstat_plot <- ggplot(causal_snp_data, aes(x = ncSNP)) +
  geom_line(aes(y = m_Bonf_Fval, color = "Bonferroni")) +
  geom_point(aes(y = m_Bonf_Fval, color = "Bonferroni")) +
  geom_ribbon(aes(ymin = lcbFstat, ymax = ucbFstat, fill = "Bonferroni"), alpha = 0.2) +
  geom_line(aes(y = m_FDR_Fval, color = "FDR")) +
  geom_point(aes(y = m_FDR_Fval, color = "FDR")) +
  geom_ribbon(aes(ymin = lcfFstat, ymax = ucfFstat, fill = "FDR"), alpha = 0.2) +
  scale_color_manual(name = "Mean F-value ± 95% CI", values = c("Bonferroni" = "blue", "FDR" = "red")) +
  scale_fill_manual(name = "Mean F-value ± 95% CI", values = c("Bonferroni" = "blue", "FDR" = "red")) +
  labs(title = "F-value at Bonferroni and FDR thresholds",
       x = "nCausal SNPs",
       y = "F-value") +
  theme_minimal() + 
  theme(axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.5,"in")
  ) + 
  ylim(min_fval, max_fval)

# Causal MAF
cmaf_auc_plot <- ggplot(cMAF_data, aes(x = causal_maf)) +
  geom_line(aes(y = mAUC, color = "AUC")) +
  geom_point(aes(y = mAUC, color = "AUC")) +
  geom_ribbon(aes(ymin = lcAUC, ymax = ucAUC, fill = "AUC"), alpha = 0.2) +
  scale_color_manual(name = "Mean AUC ± 95% CI", values = c("AUC" = "blue")) +
  scale_fill_manual(name = "Mean AUC ± 95% CI", values = c("AUC" = "blue")) +
  labs(title = "AUC vs. Causal MAF", x = "Causal MAF", y = "AUC") +
  theme_minimal() +
  theme(axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.5,"in")
  ) + 
  ylim(min_auc, max_auc)
cMAF_fstat_plot <- ggplot(cMAF_data, aes(x = causal_maf)) +
  geom_line(aes(y = mbFstat, color = "Bonferroni")) +
  geom_point(aes(y = mbFstat, color = "Bonferroni")) +
  geom_ribbon(aes(ymin = lcbFstat, ymax = ucbFstat, fill = "Bonferroni"), alpha = 0.2) +
  geom_line(aes(y = mfFstat, color = "FDR")) +
  geom_point(aes(y = mfFstat, color = "FDR")) +
  geom_ribbon(aes(ymin = lcfFstat, ymax = ucfFstat, fill = "FDR"), alpha = 0.2) +
  scale_color_manual(name = "Mean F-value ± 95% CI", values = c("Bonferroni" = "blue", "FDR" = "red")) +
  scale_fill_manual(name = "Mean F-value ± 95% CI", values = c("Bonferroni" = "blue", "FDR" = "red")) +
  labs(title = "F-value at Bonferroni and FDR thresholds",
       x = "Causal MAF",
       y = "F-value") +
  theme_minimal() + 
  theme(axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.5,"in")
  ) + 
  ylim(min_fval, max_fval)

# Heritability
herit_auc_plot <- ggplot(herit_data, aes(x = sHerit)) +
  geom_line(aes(y = m_AUC_val, color = "AUC")) +
  geom_point(aes(y = m_AUC_val, color = "AUC")) +
  geom_ribbon(aes(ymin = lcAUC, ymax = ucAUC, fill = "AUC"), alpha = 0.2) +
  scale_color_manual(name = "Mean AUC ± 95% CI", values = c("AUC" = "blue")) +
  scale_fill_manual(name = "Mean AUC ± 95% CI", values = c("AUC" = "blue")) +
  labs(title = "AUC vs. Trait Heritability", x = "Trait Heritability", y = "AUC") +
  theme_minimal() +
  theme(axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.5,"in")
  ) + 
  ylim(min_auc, max_auc)
herit_fstat_plot <- ggplot(herit_data, aes(x = sHerit)) +
  geom_line(aes(y = m_Bonf_Fval, color = "Bonferroni")) +
  geom_point(aes(y = m_Bonf_Fval, color = "Bonferroni")) +
  geom_ribbon(aes(ymin = lcbFstat, ymax = ucbFstat, fill = "Bonferroni"), alpha = 0.2) +
  geom_line(aes(y = m_FDR_Fval, color = "FDR")) +
  geom_point(aes(y = m_FDR_Fval, color = "FDR")) +
  geom_ribbon(aes(ymin = lcfFstat, ymax = ucfFstat, fill = "FDR"), alpha = 0.2) +
  scale_color_manual(name = "Mean F-value ± 95% CI", values = c("Bonferroni" = "blue", "FDR" = "red")) +
  scale_fill_manual(name = "Mean F-value ± 95% CI", values = c("Bonferroni" = "blue", "FDR" = "red")) +
  labs(title = "F-value at Bonferroni and FDR thresholds",
       x = "Trait Heritability",
       y = "F-value") +
  theme_minimal() + 
  theme(axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.5,"in")
  ) + 
  ylim(min_fval, max_fval)

# Case Control Ratio
ccratio_auc_plot <- ggplot(cascntrl_data, aes(x = cc_ratio)) +
  geom_line(aes(y = m_AUC_val, color = "AUC")) +
  geom_point(aes(y = m_AUC_val, color = "AUC")) +
  geom_ribbon(aes(ymin = lcAUC, ymax = ucAUC, fill = "AUC"), alpha = 0.2) +
  scale_color_manual(name = "Mean AUC ± 95% CI", values = c("AUC" = "blue")) +
  scale_fill_manual(name = "Mean AUC ± 95% CI", values = c("AUC" = "blue")) +
  labs(title = "AUC vs. Case/Control Ratio", x = "Case/Control Ratio", y = "AUC") +
  theme_minimal() +
  theme(axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.5,"in")
  ) + 
  ylim(min_auc, max_auc)
ccratio_fstat_plot <- ggplot(cascntrl_data, aes(x = cc_ratio)) +
  geom_line(aes(y = m_Bonf_Fval, color = "Bonferroni")) +
  geom_point(aes(y = m_Bonf_Fval, color = "Bonferroni")) +
  geom_ribbon(aes(ymin = lcbFstat, ymax = ucbFstat, fill = "Bonferroni"), alpha = 0.2) +
  geom_line(aes(y = m_FDR_Fval, color = "FDR")) +
  geom_point(aes(y = m_FDR_Fval, color = "FDR")) +
  geom_ribbon(aes(ymin = lcfFstat, ymax = ucfFstat, fill = "FDR"), alpha = 0.2) +
  scale_color_manual(name = "Mean F-value ± 95% CI", values = c("Bonferroni" = "blue", "FDR" = "red")) +
  scale_fill_manual(name = "Mean F-value ± 95% CI", values = c("Bonferroni" = "blue", "FDR" = "red")) +
  labs(title = "Bonferroni and FDR F-statistics with Confidence Intervals",
       x = "Case/Control Ratio",
       y = "F-value") +
  theme_minimal() + 
  theme(axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.5,"in")
  ) + 
  ylim(min_fval, max_fval)

# Join the AUC plots
indv_auc_plot <- indv_auc_plot + theme(legend.position = "none", plot.title = element_blank())
causal_snp_auc_plot <- causal_snp_auc_plot + theme(legend.position = "none", plot.title = element_blank())
herit_auc_plot <- herit_auc_plot + theme(legend.position = "none", plot.title = element_blank())
ccratio_auc_plot <- ccratio_auc_plot + theme(legend.position = "none", plot.title = element_blank())
cmaf_auc_plot <- cmaf_auc_plot + theme(legend.position = "none", plot.title = element_blank())
snp_auc_plot <- snp_auc_plot + theme(legend.position = "none", plot.title = element_blank())
legend <- get_legend(
  indv_auc_plot + 
    theme(legend.position = "right", 
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10),
          legend.key.size = unit(0.4, "in"))
)
combined_auc_plot <- plot_grid(
  ggdraw() + draw_label("AUC Modulation", 
                        fontface = 'bold', x = 0.5, hjust = 0.5, size = 16),
  plot_grid(
    plot_grid(indv_auc_plot, causal_snp_auc_plot, herit_auc_plot,
              ccratio_auc_plot, cmaf_auc_plot, snp_auc_plot,
              ncol = 3, align = "hv"),
    legend, 
    ncol = 2, 
    rel_widths = c(3, 0.4)
  ),
  ncol = 1,
  rel_heights = c(0.1, 1)
)
ggsave(filename = paste0(directory,"/",response,"_model_auc.png"), plot = combined_auc_plot, device = "png", bg = 'white')


# Join the Fval plots
indv_fstat_plot <- indv_fstat_plot + theme(legend.position = "none", plot.title = element_blank())
causal_snp_fstat_plot <- causal_snp_fstat_plot + theme(legend.position = "none", plot.title = element_blank())
herit_fstat_plot <- herit_fstat_plot + theme(legend.position = "none", plot.title = element_blank())
ccratio_fstat_plot <- ccratio_fstat_plot + theme(legend.position = "none", plot.title = element_blank())
cMAF_fstat_plot <- cMAF_fstat_plot + theme(legend.position = "none", plot.title = element_blank())
snp_fstat_plot <- snp_fstat_plot + theme(legend.position = "none", plot.title = element_blank())

# Extract the legend from one of the plots
legend <- get_legend(
  indv_fstat_plot + 
    theme(legend.position = "right", 
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10),
          legend.key.size = unit(0.4, "in"))
)

# Add the common title at the top
combined_fval_plot <- plot_grid(
  ggdraw() + draw_label("F-value Modulation", 
                        fontface = 'bold', x = 0.5, hjust = 0.5, size = 16),
  plot_grid(
    plot_grid(indv_fstat_plot, causal_snp_fstat_plot, herit_fstat_plot,
              ccratio_fstat_plot, cMAF_fstat_plot, snp_fstat_plot,
              ncol = 3, align = "hv"),
    legend, 
    ncol = 2, 
    rel_widths = c(3, 0.4)
  ),
  ncol = 1,
  rel_heights = c(0.1, 1)
)
ggsave(filename = paste0(directory,"/",response,"_model_fval.png"), plot = combined_fval_plot, device = "png", bg = 'white')

# Combine both figures
label_A <- ggdraw() + draw_label("A", fontface = 'bold', x = 0, hjust = 0, size = 20)
label_B <- ggdraw() + draw_label("B", fontface = 'bold', x = 0, hjust = 0, size = 20)
combined_plot <- plot_grid(
  plot_grid(label_A, combined_auc_plot, ncol = 1, rel_heights = c(0.05, 1)),
  plot_grid(label_B, combined_fval_plot, ncol = 1, rel_heights = c(0.05, 1)),
  ncol = 1,
  rel_heights = c(1, 1)
)
ggsave(filename = paste0(directory,"/",response,"_model_aucfval.png"), plot = combined_plot, device = "png", bg = 'white')

# Build the mixed results
assoc_stats_list <- list()
auc_list <- list()
for(data in list(mixed_Rdata,mixed_Rdata2,mixed_Rdata3,mixed_Rdata4,mixed_Rdata5)){
  for (name in names(data)) {
    if (startsWith(name, "simgene_")) {
      
      # Extract the relevant values from the name using regular expressions
      n_individuals <- as.numeric(str_extract(name, "(?<=simgene_)\\d+"))
      n_snps <- as.numeric(str_extract(name, "(?<=inds_)(\\d+e[+-]?\\d+|\\d+)"))
      causal_snps <- as.numeric(str_extract(name, "(?<=snps_)\\d+"))
      cc_ratio <- as.numeric(str_extract(name, "(?<=causal_)\\d+\\.\\d+"))
      causal_maf <- as.numeric(str_extract(name, "(?<=ccratio_)\\d+\\.\\d+"))
      heritability <- as.numeric(str_extract(name, "(?<=maf_)\\d+\\.\\d+"))
      stats_track <- data.frame(
        n_individuals = n_individuals,
        n_snps = n_snps,
        causal_snps = causal_snps,
        cc_ratio = cc_ratio,
        causal_maf = causal_maf,
        heritability = heritability
      )
      
      # Extract the assoc_stats data frame
      assoc_stats <- data[[name]]$assoc_stats
      # Extract AUC
      auc <- as.numeric(data[[name]]$auc)
      
      
      # Add these as new columns to the assoc_stats data frame
      assoc_stats_list[[name]] <- cbind(assoc_stats,stats_track)
      auc_list[[name]] <- cbind(auc,stats_track)
    }
  }
}
mixed_data_F_stats <- bind_rows(assoc_stats_list) %>% filter(!is.na(F_statistic))
mixed_data_AUC_stats <- bind_rows(auc_list)

# Define the data and response we are interested in
data <- mixed_data_AUC_stats
response <- "auc"

# Scale the data
interest_vals <- data[[response]] 
total_data <- data %>%
  dplyr::mutate(
    n_individuals_scaled = c(scale(n_individuals)),
    n_snps_scaled = c(scale(n_snps)),
    causal_snps_scaled = c(scale(causal_snps)),
    cc_ratio_scaled = c(scale(cc_ratio)),
    causal_maf_scaled = c(scale(causal_maf)),
    heritability_scaled = c(scale(heritability)),
    response = interest_vals,
    response_scaled = c(scale(response))
  ) %>%
  dplyr::select(response, n_individuals,
                n_snps,causal_snps,cc_ratio, causal_maf, 
                heritability,response_scaled, n_individuals_scaled,
                n_snps_scaled, causal_snps_scaled,
                cc_ratio_scaled,causal_maf_scaled, heritability_scaled)

# Generate a plot for all terms
# Function to calculate and extract R-squared from the linear model
get_r_squared <- function(data, x, y) {
  model <- lm(as.formula(paste(y, "~", x)), data = data)
  summary(model)$r.squared
}

create_plot <- function(data, x, y, title, xlab, y_label) {
  # Perform a Pearson Correlation test
  cor_test_result <- cor.test(data[[x]], data[[y]], method = "pearson")
  r_squared <- round(cor_test_result$estimate^2, 2)
  p_value <- format(cor_test_result$p.value, scientific = TRUE, digits = 2)
  xlab_with_r2_p <- paste0(xlab, " (R² = ", r_squared, ", p = ", p_value, ")")
  ggplot(data, aes_string(x = x, y = y)) + 
    geom_point() + 
    geom_smooth(method = "lm", col = "red") + 
    labs(title = title) + 
    xlab(xlab_with_r2_p) + 
    theme_minimal() + 
    theme(axis.title.y = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          plot.title = element_text(hjust = 0.5, size = 16),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12)) + 
    ylab(y_label)
}


n_ind_plot <- create_plot(total_data, "n_individuals_scaled", "response_scaled","AUC vs nIndividuals" ,"nIndividuals", "AUC")
n_snp_plot <- create_plot(total_data, "n_snps_scaled", "response_scaled", "AUC vs nSNPs", "nSNPs", "AUC")
n_causal_plot <- create_plot(total_data, "causal_snps_scaled", "response_scaled", "AUC vs nCausal SNPs", "nCausal SNPs", "AUC")
n_case_plot <- create_plot(total_data, "cc_ratio_scaled", "response_scaled", "AUC vs nCases", "nCases", "AUC")
s_maf_plot <- create_plot(total_data, "causal_maf_scaled", "response_scaled", "AUC vs Causal MAF", "Causal MAF", "AUC")
s_herit_plot <- create_plot(total_data, "heritability_scaled", "response_scaled", "AUC vs Heritability", "Heritability", "AUC")

s_response_plot <- ggplot(total_data, aes(x = response)) +
  geom_histogram(binwidth = 0.05, fill = 'blue', color = 'black', alpha = 0.7) +
  labs(title = paste0("Histogram of ",toupper(response)), x = toupper(response), y = "Frequency") +
  theme_minimal() + 
  theme(axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)) + 
  xlab(paste0(toupper(response)))

layout_matrix <- rbind(c(1, 2),c(3, 4),c(5, 6),c(7, 7),c(7, 7))
heights <- c(1, 1, 1, 0.9, 0.9)
title_grob <- textGrob("Scaled Associations", gp = gpar(fontsize = 18, fontface = "bold"))
initial_sum_plots <- grid.arrange(
  n_ind_plot, n_snp_plot, n_case_plot, n_causal_plot, s_maf_plot, s_herit_plot, s_response_plot,
  layout_matrix = layout_matrix, top = title_grob, heights = heights
)
ggsave(filename = paste0(directory,"/",response,"_allmixed_assocs.png"), plot = initial_sum_plots, device = "png", bg = 'white')

# Select just the scaled correlations
data_for_correlation <- total_data %>%
  dplyr::select(ends_with("_scaled"))

# Create a correlation plot
correlation_matrix <- cor(data_for_correlation, method = 'pearson')
correlation_data <- as.data.frame(correlation_matrix)
correlation_data$row <- rownames(correlation_data)
correlation_long <- pivot_longer(correlation_data, cols = -row, names_to = "variable", values_to = "value") %>%
  arrange(row) %>%
  mutate(
    row = str_replace(row, "_scaled", ""),
    variable = str_replace(variable, "_scaled", "")
  )
correlation_long$row <- factor(correlation_long$row, levels = rev(unique(correlation_long$row)))
label_mapping <- c(
  "response" = "AUC",
  "n_snps" = "nSNPs",
  "n_individuals" = "nIndividuals",
  "heritability" = "Heritability",
  "cc_ratio" = "nCases",
  "causal_snps" = "nCausal SNPs",
  "causal_maf" = "Causal MAF"
)

# Apply the mapping to relabel rows and columns
correlation_long <- correlation_long %>%
  mutate(
    row = label_mapping[row],
    variable = label_mapping[variable]
  )

correlation_long %>%
  filter(variable == "AUC")

# Create the heatmap with circles
correlation_heatmap <- ggplot(correlation_long, aes(x = row, y = as.factor(variable), fill = value, size = abs(value))) +
  geom_point(shape = 21, color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-1, 1), space = "Lab",
                       name = str_wrap("Pearson Correlation",width = 10)) +
  scale_size(range = c(8, 20), guide = "none") +
  labs(title = "Response/Predictor Correlation Heatmap") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                   size = 12, hjust = 1),
        axis.text.y = element_text(size = 12),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
ggsave(filename = paste0(directory,"/",response,"_cor_heatmap.png"), plot = correlation_heatmap, device = "png", bg = 'white')


# Combine with R^2 values
n_ind_plot <- n_ind_plot + theme(plot.title = element_blank())
n_snp_plot <- n_snp_plot + theme(plot.title = element_blank())
n_case_plot <- n_case_plot + theme(plot.title = element_blank())
n_causal_plot <- n_causal_plot + theme(plot.title = element_blank())
s_maf_plot <- s_maf_plot + theme(plot.title = element_blank())
s_herit_plot <- s_herit_plot + theme(plot.title = element_blank())
correlation_heatmap <- correlation_heatmap + theme(plot.title = element_blank())

# Create a title plot for the first section
title_plot1 <- ggdraw() + 
  draw_label("Individual Variable Linear Regression", fontface = 'bold', size = 14, x = 0.5, hjust = 0.5)

# Create a title plot for the second section
title_plot2 <- ggdraw() + 
  draw_label("Pearson Correlation Heatmap", fontface = 'bold', size = 14, x = 0.5, hjust = 0.5)

# Create labels
label_A <- ggdraw() + draw_label("A", fontface = 'bold', x = 0, hjust = 0, size = 20)
label_B <- ggdraw() + draw_label("B", fontface = 'bold', x = 0, hjust = 0, size = 20)

# Combine the six individual plots into a 3x2 grid
plot_rsqrs <- plot_grid(
  n_ind_plot, n_snp_plot, n_case_plot,
  n_causal_plot, s_maf_plot, s_herit_plot,
  ncol = 3, nrow = 2
)
mixed_initial_plot <- plot_grid(
  plot_grid(label_A, title_plot1, plot_rsqrs, ncol = 1, rel_heights = c(0.1, 0.1, 1)),
  plot_grid(label_B, title_plot2, correlation_heatmap, ncol = 1, rel_heights = c(0.1, 0.1, 1)),
  ncol = 1,
  rel_heights = c(0.8, 1) # Adjust heights as necessary
)
ggsave(filename = paste0(directory,"/",response,"_allmixed_assocs.png"), plot = mixed_initial_plot, device = "png", bg = 'white')

# Join two datasets
auc_fval_stats <- plot_grid(
  mixed_initial_plot_auc,
  patchwork::plot_spacer(),
  mixed_initial_plot_Fval,
  rel_widths = c(1,0.02,1),
  ncol = 3
)
ggsave(filename = paste0("~/Desktop/University/Biomedicine/Y3/Dissertation/PLINK_Tutorials/Simulate_Data_Tes/auc_fval_allmixed_assocs.png"), 
       plot = auc_fval_stats, device = "png", bg = 'white',
       width = 20, height = 15, units = "in")



# Model data

# Partition data
set.seed(1)
trainIndex <- createDataPartition(data_for_correlation$response_scaled, p = 0.8, 
                                  list = FALSE, 
                                  times = 1)

# Split the data into training and testing sets
data_train <- data_for_correlation[trainIndex, ]
data_test <- data_for_correlation[-trainIndex, ]

# Perform cross validation of each model
train_control <- trainControl(method="cv", number=10)

# Train each model using cross-validation and collect stats
# Linear model
set.seed(1)
lm_cv <- train(response_scaled ~ ., data=data_train, method="lm", trControl=train_control)

# Polynomial model
degree_indv<- length(unique(data_train$n_individuals_scaled))
degree_snps<- length(unique(data_train$n_snps_scaled))
degree_csnps<- length(unique(data_train$causal_snps_scaled))
degree_cindv<- length(unique(data_train$cc_ratio_scaled))
degree_cmaf<- length(unique(data_train$causal_maf_scaled))
degree_herit<- length(unique(data_train$heritability_scaled))

poly_formula_str <- paste("response_scaled ~",
                          paste("poly(n_individuals_scaled,", degree_indv-1, ")", sep = ""),
                          paste("poly(n_snps_scaled,", degree_snps-1, ")", sep = ""),
                          paste("poly(causal_snps_scaled,", degree_csnps-1, ")", sep = ""),
                          paste("poly(cc_ratio_scaled,", degree_cindv-1, ")", sep = ""),
                          paste("poly(causal_maf_scaled,", degree_cmaf-1, ")", sep = ""),
                          paste("poly(heritability_scaled,", degree_herit-1, ")", sep = ""),
                          sep = " + ")
poly_formula <- as.formula(poly_formula_str)
set.seed(1)
poly_cv <- train(poly_formula, data=data_train, method="lm", trControl=train_control)

# Regression tree model
tree_grid <- expand.grid(
  cp = seq(0.01, 0.1, by = 0.01)
)
set.seed(1)
tree_cv <- train(response_scaled ~ ., data = data_train, method = "rpart", trControl=train_control,
                 tuneGrid = tree_grid)

# Random forest model
rf_grid <- expand.grid(
  mtry = c(2, 4, 6)
)
set.seed(1)
rf_cv <- train(response_scaled ~ ., data=data_train, method="rf", trControl=train_control,
               tuneGrid = rf_grid)

# Gradient boosting regression
gbm_grid <- expand.grid(
  interaction.depth = c(1,2,3,4),
  n.trees = c(50, 100, 200, 400, 600),
  shrinkage = c(0.01,0.1),
  n.minobsinnode = 10
)
set.seed(1)
gbm_cv <- train(response_scaled ~ ., data=data_train, method="gbm", trControl=train_control, verbose=FALSE,
                tuneGrid = gbm_grid,
                distribution = "laplace")


# Collect model stats
resamples_list <- resamples(list(
  lm = lm_cv,
  poly = poly_cv,
  tree = tree_cv,
  rf = rf_cv,
  gbm = gbm_cv
))
resample_summary <- summary(resamples_list)

## Check for normality of differences using a shapiro wilks test
check_normality <- function(diffs) {
  qqnorm(diffs)
  qqline(diffs, col = "red")
  return(shapiro.test(diffs))
}

lm_rmse <- resamples_list$values$`lm~RMSE`
poly_rmse <- resamples_list$values$`poly~RMSE`
tree_rmse <- resamples_list$values$`tree~RMSE`
rf_rmse <- resamples_list$values$`rf~RMSE`
gbm_rmse <- resamples_list$values$`gbm~RMSE`

rmse_diffs <- list(
  lm_poly = lm_rmse - poly_rmse,
  lm_tree = lm_rmse - tree_rmse,
  lm_rf = lm_rmse - rf_rmse,
  lm_gbm = lm_rmse - gbm_rmse,
  poly_tree = poly_rmse - tree_rmse,
  poly_rf = poly_rmse - rf_rmse,
  poly_gbm = poly_rmse - gbm_rmse,
  tree_rf = tree_rmse - rf_rmse,
  tree_gbm = tree_rmse - gbm_rmse,
  rf_gbm = rf_rmse - gbm_rmse
)


rmse_normality <- lapply(rmse_diffs, check_normality)
for(test in names(rmse_normality)){
  if(rmse_normality[[test]]$p.value < 0.05){
    print(paste0("Significant difference in ",test, " p value = ",rmse_normality[[test]]$p.value))
  }
}

# Build a dataframe of the W and P statistics for normality
rmse_results <- data.frame(predictor_1 = character(),
                           predictor_2 = character(),
                           RMSE_W_value = numeric(),
                           RMSE_p_value = numeric(),
                           stringsAsFactors = FALSE)
for (name in names(rmse_normality)) {
  parts <- unlist(strsplit(name, "_"))
  W_value <- rmse_normality[[name]]$statistic
  p_value <- rmse_normality[[name]]$p.value
  rmse_results <- rbind(rmse_results, data.frame(predictor_1 = parts[1],
                                                 predictor_2 = parts[2],
                                                 RMSE_W_value = W_value,
                                                 RMSE_p_value = p_value,
                                                 stringsAsFactors = FALSE))
}

lm_rsqr <- resamples_list$values$`lm~Rsquared`
poly_rsqr <- resamples_list$values$`poly~Rsquared`
tree_rsqr <- resamples_list$values$`tree~Rsquared`
rf_rsqr <- resamples_list$values$`rf~Rsquared`
gbm_rsqr <- resamples_list$values$`gbm~Rsquared`

rsqu_diffs <- list(
  lm_poly = lm_rsqr - poly_rsqr,
  lm_tree = lm_rsqr - tree_rsqr,
  lm_rf = lm_rsqr - rf_rsqr,
  lm_gbm = lm_rsqr - gbm_rsqr,
  poly_tree = poly_rsqr - tree_rsqr,
  poly_rf = poly_rsqr - rf_rsqr,
  poly_gbm = poly_rsqr - gbm_rsqr,
  tree_rf = tree_rsqr - rf_rsqr,
  tree_gbm = tree_rsqr - gbm_rsqr,
  rf_gbm = rf_rsqr - gbm_rsqr
)

rsqur_normality <- lapply(rsqu_diffs, check_normality)
for(test in names(rsqur_normality)){
  if(rsqur_normality[[test]]$p.value < 0.05){
    print(paste0("Significant difference in ",test, " p value = ",rsqur_normality[[test]]$p.value))
  }
}

# Build a dataframe of the W and P statistics for normality
rsqr_results <- data.frame(predictor_1 = character(),
                           predictor_2 = character(),
                           Rsq_W_value = numeric(),
                           Rsq_p_value = numeric(),
                           stringsAsFactors = FALSE)
for (name in names(rmse_normality)) {
  parts <- unlist(strsplit(name, "_"))
  W_value <- rmse_normality[[name]]$statistic
  p_value <- rmse_normality[[name]]$p.value
  rsqr_results <- rbind(rsqr_results, data.frame(predictor_1 = parts[1],
                                                 predictor_2 = parts[2],
                                                 Rsq_W_value = W_value,
                                                 Rsq_p_value = p_value,
                                                 stringsAsFactors = FALSE))
}

# Perform pariwise t-test comparisons of resampling results
diffs <- diff(resamples_list)

# By default summary 
significant_differences <- summary(diffs, adjust = "bonferroni")


# Extract RMSE and R-squared values and calculate mean and 95% CI
summary_stats <- list()
for(measure in c("RMSE","Rsquared")){
  # Get confidence intervals
  model_summary <- resamples_list$values %>%
    select(contains(paste0(measure))) %>%
    pivot_longer(everything(), names_to = "Model", values_to = "Value") %>%
    mutate(Model = gsub(paste0("~",measure), "", Model)) %>%
    group_by(Model) %>%
    summarize(
      Mean = mean(Value),
      Lower_CI = mean(Value) - qt(1 - 0.05 / 2, n() - 1) * sd(Value) / sqrt(n()),
      Upper_CI = mean(Value) + qt(1 - 0.05 / 2, n() - 1) * sd(Value) / sqrt(n()),
      SE = sd(Value) / sqrt(n())
    ) %>% {
      if(measure == "RMSE"){
        arrange(.,Mean) %>%
          mutate(rank = row_number())
      } else {
        arrange(.,desc(Mean)) %>%
          mutate(rank = row_number())
      }
    }
  model_summary$Model <- factor(model_summary$Model, levels = model_summary$Model)
  summary_stats[[paste0(measure)]] <- model_summary
  model_rank <- model_summary %>% select(Model,Upper_CI,rank)
  # Get significance values 
  model_matrix <- significant_differences$table[[paste0(measure)]]
  model_matrix[upper.tri(model_matrix)] <- ""
  signif_comparisons <- model_matrix %>%
    as.data.frame() %>%
    rownames_to_column(var = "Model1") %>%
    pivot_longer(cols = -Model1, names_to = "Model2", values_to = "p_value") %>%
    filter(p_value != "") %>%
    mutate(p_value = as.numeric(p_value)) %>%
    mutate(significance = case_when(
      p_value < 0.0001 ~ "****",
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*"
    )) %>%
    filter(!if_any(everything(), is.na)) %>%
    select(Model1, Model2, p_value, significance) %>%
    left_join(model_rank, by = c("Model1" = "Model")) %>%
    left_join(model_rank, by = c("Model2" = "Model")) %>%
    mutate(priCI = ifelse(Upper_CI.x > Upper_CI.y, Upper_CI.x, Upper_CI.y),
           secCI = ifelse(Upper_CI.x < Upper_CI.y, Upper_CI.x, Upper_CI.y)) %>%
    arrange(priCI) %>%
    group_by(secCI) %>%
    mutate(secrank = cur_group_id()) %>%
    ungroup() %>%
    arrange(priCI,secCI) 
  # Get coordinates for sig lines
  comparisons_list <- list()
  annotations <- c()
  y_positions <- c()
  base_y_position <- max(model_summary$Upper_CI) * 1
  max_ci_track <- c()
  prev_pos <- NULL
  for (i in 1:nrow(signif_comparisons)) {
    comparisons_list[[i]] <- c(signif_comparisons$Model1[i], signif_comparisons$Model2[i])
    annotations[i] <- signif_comparisons$significance[i]
    priCI <- signif_comparisons[i,]$priCI
    secrank <- signif_comparisons[i,]$secrank
    predict_pos <- priCI + 0.04 * secrank
    if(secrank != 1){
      if(prev_pos > priCI){
        predict_pos <- prev_pos + 0.06
      }
    } else {
      if(i != 1){
        if(prev_pos + 0.06 > priCI){
          predict_pos <- prev_pos + 0.06
        } else {
          predict_pos <- predict_pos - 0.02 * secrank # Initial elevation above the plot
        }
      } else {
        predict_pos <- predict_pos - 0.02 * secrank # Initial elevation above the plot
      }
    }
    prev_pos <- predict_pos
    
    y_positions[i] <- predict_pos
  }
  summary_stats[[paste0(measure)]] <- list(model_summary = model_summary,
                                           model_pvals = signif_comparisons,
                                           comparisons_list = comparisons_list,
                                           annotations = annotations,
                                           y_positions = y_positions)
}


# Plot RMSE with error bars and significance annotations
rmse_summary <- summary_stats$RMSE$model_summary
rmse_lables <- rmse_summary %>% 
  mutate(Model = case_when(
    str_detect(Model, "lm") ~ "Linear",
    str_detect(Model, "poly") ~ "Polynomial",
    str_detect(Model, "tree") ~ "Regression Tree",
    str_detect(Model, "rf") ~ "Random Forest",
    str_detect(Model, "gbm") ~ "Gradient Boosting",
    TRUE ~ Model
  )) %>% select(Model) %>% as.matrix() %>% as.character()
y_breaks <- pretty(c(0, max(rmse_summary$Upper_CI) * 1.01), n = 5)
y_labels <- c(0,0.2,0.4,0.6,0.8,1.0,1.2,1.4)[1:length(y_breaks)]
rmse_sig_plot <- ggplot(rmse_summary, aes(x = Model, y = Mean)) +
  geom_bar(stat = "identity", fill = "grey", color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.2) +
  labs(title = "AUC Residual Error", y = "Mean RMSE", x = "Model") +
  theme_minimal() +
  geom_signif(
    comparisons = summary_stats$RMSE$comparisons_list,
    annotations = summary_stats$RMSE$annotations,
    y_position = summary_stats$RMSE$y_positions,
    tip_length = 0.03,
    textsize = 4,
    vjust = 0.4
  ) + 
  scale_y_continuous(
    breaks = y_breaks, 
    labels = y_labels,
    expand = expansion(mult = c(0, 0.2))
  ) +
  scale_x_discrete(labels = rmse_lables) + 
  theme(
    plot.margin = unit(c(3, 1, 1, 1), "lines"),
    panel.grid.major.y = element_line(color = "grey", size = 0.5),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y = element_text(hjust = 0.25, vjust = 2.5, size = 14),
    axis.title.x = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12)
  )

# Plot Rsquared with error bars and significance annotations
rsquare_summary <- summary_stats$Rsquared$model_summary
rsqr_lables <- rsquare_summary %>% 
  mutate(Model = case_when(
    str_detect(Model, "lm") ~ "Linear",
    str_detect(Model, "poly") ~ "Polynomial",
    str_detect(Model, "tree") ~ "Regression Tree",
    str_detect(Model, "rf") ~ "Random Forest",
    str_detect(Model, "gbm") ~ "Gradient Boosting",
    TRUE ~ Model
  )) %>% select(Model) %>% as.matrix() %>% as.character()
y_breaks <- pretty(c(0, max(rsquare_summary$Upper_CI) * 1.01), n = 5)
y_labels <- c(0,0.2,0.4,0.6,0.8,1.0,1.2,1.4)[1:length(y_breaks)]
rsquare_sig_plot <- ggplot(rsquare_summary, aes(x = Model, y = Mean)) + 
  geom_bar(stat = "identity", fill = "grey", color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.2) +
  labs(title = "AUC Varience Explained", y = "Mean Rsquared", x = "Model") +
  theme_minimal() +
  geom_signif(
    comparisons = summary_stats$Rsquared$comparisons_list,
    annotations = summary_stats$Rsquared$annotations,
    y_position = summary_stats$Rsquared$y_positions,
    tip_length = 0.03,
    textsize = 4,
    vjust = 0.4
  ) + 
  scale_y_continuous(
    breaks = y_breaks, 
    labels = y_labels,
    expand = expansion(mult = c(0, 0.2))
  ) +
  scale_x_discrete(labels = rsqr_lables) + 
  theme(
    plot.margin = unit(c(3, 1, 1, 1), "lines"),
    panel.grid.major.y = element_line(color = "grey", size = 0.5),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y = element_text(hjust = 0.25,vjust = 2.5, size = 14),
    axis.title.x = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12)
  )

residual_rsqu_plot <- plot_grid(
  rmse_sig_plot,
  rsquare_sig_plot,
  ncol = 2,
  rel_heights = c(1,1),
  labels = c("A","B"),
  label_size = 20
)
ggsave(filename = paste0(directory,"/",response,"_model_rmsersq.png"), plot = residual_rsqu_plot, device = "png", bg = 'white')


# Get the importance of each model
standardise_col <- function(df, column_name) {
  total_sum <- sum(df[[column_name]], na.rm = TRUE)
  df[[column_name]] <- df[[column_name]] / total_sum
  return(df)
}

# Linear model
varImp_lm <- as.data.frame(varImp(lm_cv)$importance)
varImp_lm$Predictor <- rownames(varImp_lm)
colnames(varImp_lm)[1] <- "lm_importance"
varImp_lm <- standardise_col(varImp_lm, "lm_importance")

# Polynomial Model
poly_coefs <- coef(poly_cv$finalModel)
poly_terms <- names(poly_coefs)
poly_predictors <- gsub("`?poly\\(([^,]+),.*", "\\1", poly_terms)
varImp_poly <- data.frame(Predictor = poly_predictors, poly_importance = abs(poly_coefs)) %>%
  group_by(Predictor) %>%
  summarise(poly_importance = sum(poly_importance)) %>%
  as.data.frame()
varImp_poly <- standardise_col(varImp_poly, "poly_importance")

# Regression Tree
varImp_tree <- as.data.frame(varImp(tree_cv)$importance)
varImp_tree$Predictor <- rownames(varImp_tree)
colnames(varImp_tree)[1] <- "tree_importance"
varImp_tree <- standardise_col(varImp_tree, "tree_importance")

# Random Forest
varImp_rf <- as.data.frame(varImp(rf_cv)$importance)
varImp_rf$Predictor <- rownames(varImp_rf)
colnames(varImp_rf)[1] <- "rf_importance"
varImp_rf <- standardise_col(varImp_rf, "rf_importance")

# Gradient Boosting Model
varImp_gbm <- summary(gbm_cv$finalModel)
varImp_gbm <- as.data.frame(varImp_gbm)
colnames(varImp_gbm) <- c("Predictor", "gbm_importance")
varImp_gbm <- standardise_col(varImp_gbm, "gbm_importance")

importance_df <- reduce(list(varImp_lm, 
                             varImp_poly, 
                             varImp_tree, 
                             varImp_rf, 
                             varImp_gbm), 
                        function(x, y) merge(x, y, by = "Predictor", all = TRUE)) %>%
  filter(!is.na(gbm_importance)) %>%
  mutate(Predictor = gsub("_scaled", "", Predictor))

importance_long <- importance_df %>% 
  mutate(total_imp = rowSums(select(., ends_with("_importance")))) %>%
  arrange(desc(total_imp)) %>%
  select(-total_imp) %>%
  pivot_longer(., cols = ends_with("importance"), names_to = "Model", values_to = "Importance") %>%
  mutate(Model = case_when(
    str_detect(Model, "lm_") ~ "Linear",
    str_detect(Model, "poly_") ~ "Polynomial",
    str_detect(Model, "tree_") ~ "Regression Tree",
    str_detect(Model, "rf_") ~ "Random Forest",
    str_detect(Model, "gbm_") ~ "Gradient Boosting",
    TRUE ~ Model
  ))

# Get the RMSE of model
heatmap_rmse <- rmse_summary %>% 
  select(Model,"Importance" = Mean, "RMSE_rank" = rank) %>%
  mutate(Predictor = "RMSE",
         Model = case_when(
           str_detect(Model, "lm") ~ "Linear",
           str_detect(Model, "poly") ~ "Polynomial",
           str_detect(Model, "tree") ~ "Regression Tree",
           str_detect(Model, "rf") ~ "Random Forest",
           str_detect(Model, "gbm") ~ "Gradient Boosting",
           TRUE ~ Model
         )) %>%
  bind_rows(importance_long) %>%
  mutate(Type = if_else(Predictor == "RMSE", "RMSE", "Importance"))

heatmap_rmse_rsq <- rsquare_summary %>% 
  select(Model,"Importance" = Mean, "RSQ_rank" = rank) %>%
  mutate(Predictor = "RSQ",
         Model = case_when(
           str_detect(Model, "lm") ~ "Linear",
           str_detect(Model, "poly") ~ "Polynomial",
           str_detect(Model, "tree") ~ "Regression Tree",
           str_detect(Model, "rf") ~ "Random Forest",
           str_detect(Model, "gbm") ~ "Gradient Boosting",
           TRUE ~ Model
         )) %>%
  bind_rows(heatmap_rmse) %>%
  mutate(Type = if_else(Predictor == "RSQ", "RSQ", Type)) 
vertical_line_pos <- max(which(heatmap_rmse_rsq$Predictor != "RMSE" & 
                                 heatmap_rmse_rsq$Predictor != "RSQ"))

# Sort models by RMSE to ensure the lowest RMSE is at the top
rmse_order <- heatmap_rmse %>%
  arrange(desc(RMSE_rank)) %>%
  pull(Model) %>%
  unique()

heatmap_rmse_rsq$Model <- factor(heatmap_rmse_rsq$Model, 
                                 levels = rmse_order)
# Ensure the column factor levels are ordered by importance
col_order <- c(unique(importance_long$Predictor),"RMSE","RSQ")
heatmap_rmse_rsq$Predictor <- factor(heatmap_rmse_rsq$Predictor, 
                                     levels = col_order)

# Plot the heatmap with the RMSE column
importance_heat <- ggplot() +
  geom_tile(data = heatmap_rmse_rsq[heatmap_rmse_rsq$Type == "Importance",], aes(x = Predictor, y = Model, fill = Importance)) +
  scale_fill_gradient(low = "white", high = "darkred", name = "Importance") +
  new_scale("fill") + 
  geom_tile(data = heatmap_rmse_rsq[heatmap_rmse_rsq$Type == "RMSE",], aes(x = Predictor, y = Model, fill = Importance)) +
  scale_fill_gradient(low = "darkblue", high = "lightblue", name = "RMSE") +
  new_scale("fill") + 
  geom_tile(data = heatmap_rmse_rsq[heatmap_rmse_rsq$Type == "RSQ",], aes(x = Predictor, y = Model, fill = Importance)) +
  scale_fill_gradient(low = "lightgreen", high = "darkgreen", name = "RSQ") +
  theme_minimal() +
  geom_vline(xintercept = 6.55, color = "white", linewidth = 6) +
  geom_vline(xintercept = 7.5, color = "white", linewidth = 6) +
  geom_vline(xintercept = 8.45, color = "white", linewidth = 6) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14, hjust = 0.45),
        axis.title.y = element_text(size = 14),
        legend.title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold")) +
  labs(title = "Heatmap of Model Predictor Importance, RMSE, and RSQ",
       x = "Predictor",
       y = "Model") + 
  scale_x_discrete(drop = FALSE)
ggsave(filename = paste0(directory,"/",response,"_model_impheat.png"), plot = importance_heat, device = "png", bg = 'white')


# Use each model to predict the response based off the training data
predictions_lm <- predict(lm_cv, newdata = data_test)
predictions_poly <- predict(poly_cv, newdata = data_test)
predictions_tree <- predict(tree_cv, newdata = data_test)
predictions_rf <- predict(rf_cv, newdata = data_test)
predictions_gbm <- predict(gbm_cv, newdata = data_test)

# Rank models based on RMSE
rmse_lm <- data.frame(model = "lm", RMSE = sqrt(mean((data_test$response_scaled - predictions_lm)^2)))
rmse_poly <- data.frame(model = "pm", RMSE = sqrt(mean((data_test$response_scaled - predictions_poly)^2)))
rmse_tree <- data.frame(model = "tm", RMSE = sqrt(mean((data_test$response_scaled - predictions_tree)^2)))
rmse_rf <- data.frame(model = "rfm", RMSE = sqrt(mean((data_test$response_scaled - predictions_rf)^2)))
rmse_gbm <- data.frame(model = "gbm", RMSE = sqrt(mean((data_test$response_scaled - predictions_gbm)^2)))
ordered_models <- rbind(rmse_lm,rmse_poly,rmse_tree,rmse_rf,rmse_gbm) %>%
  arrange(RMSE) %>%
  mutate(rank = row_number(),
         colour = c("#00FF00","#CCFF00","#FFCC33","#FF9900","#FF0000"))

# Visually plot these models
plot_data <- data_test
plot_data$pred_lm <- predictions_lm
plot_data$pred_poly <- predictions_poly
plot_data$pred_tree <- predictions_tree
plot_data$pred_rf <- predictions_rf
plot_data$pred_gbm <- predictions_gbm
plot_data$resid_lm <- plot_data$response_scaled - plot_data$pred_lm
plot_data$resid_poly <- plot_data$response_scaled - plot_data$pred_poly
plot_data$resid_tree <- plot_data$response_scaled - plot_data$pred_tree
plot_data$resid_rf <- plot_data$response_scaled - plot_data$pred_rf
plot_data$resid_gbm <- plot_data$response_scaled - plot_data$pred_gbm

min_X <- floor(min(plot_data$response_scaled))
max_X <- ceiling(max(plot_data$response_scaled))
min_y <- floor(min(plot_data %>% select(starts_with("pred"))))
max_y <- ceiling(max(plot_data %>% select(starts_with("pred"))))

rsq_lm <- rsquare_summary %>% filter(Model == "lm") %>% select(Mean) %>% as.numeric()
plot_lm <- ggplot(plot_data, aes(x = response_scaled, y = pred_lm)) +
  geom_point(color = ordered_models %>% filter(model == "lm") %>% select(colour) %>% as.character(), alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(title = "Linear Predictions",
       x = paste0("Actual ",toupper(response)),
       y = paste0("Predicted ",toupper(response))) +
  theme_minimal() + 
  xlim(min_X, max_X) + 
  ylim(min_y, max_y) + 
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  ) + 
  annotate("text", x = 3, y = -1.7, label = paste("R²:", round(rsq_lm, 4)),
           hjust = 1, vjust = 1, size = 4, color = "black", fontface = "bold")
rsq_pm <- rsquare_summary %>% filter(Model == "poly") %>% select(Mean) %>% as.numeric()
plot_poly <- ggplot(plot_data, aes(x = response_scaled, y = pred_poly)) +
  geom_point(color = ordered_models %>% filter(model == "pm") %>% select(colour) %>% as.character(), alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(title = "Polynomial Predictions",
       x = paste0("Actual ",toupper(response)),
       y = paste0("Predicted ",toupper(response))) +
  theme_minimal() + 
  xlim(min_X, max_X) + 
  ylim(min_y, max_y) + 
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    # axis.title.x = element_text(size = 14),
    # axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  ) + 
  annotate("text", x = 3, y = -1.7, label = paste("R²:", round(rsq_pm, 4)),
           hjust = 1, vjust = 1, size = 4, color = "black", fontface = "bold")
rsq_tm <- rsquare_summary %>% filter(Model == "tree") %>% select(Mean) %>% as.numeric()
plot_tree <- ggplot(plot_data, aes(x = response_scaled, y = pred_tree)) +
  geom_point(color = ordered_models %>% filter(model == "tm") %>% select(colour) %>% as.character(), alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(title = "Regression Tree Predictions",
       x = paste0("Actual ",toupper(response)),
       y = paste0("Predicted ",toupper(response))) +
  theme_minimal() + 
  xlim(min_X, max_X) + 
  ylim(min_y, max_y) + 
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    # axis.title.x = element_text(size = 14),
    # axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  ) + 
  annotate("text", x = 3, y = -1.7, label = paste("R²:", round(rsq_tm, 4)),
           hjust = 1, vjust = 1, size = 4, color = "black", fontface = "bold")
rsq_rf <- rsquare_summary %>% filter(Model == "rf") %>% select(Mean) %>% as.numeric()
plot_rf <- ggplot(plot_data, aes(x = response_scaled, y = pred_rf)) +
  geom_point(color = ordered_models %>% filter(model == "rfm") %>% select(colour) %>% as.character(), alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(title = "Random Forest Predictions",
       x = paste0("Actual ",toupper(response)),
       y = paste0("Predicted ",toupper(response))) +
  theme_minimal() + 
  xlim(min_X, max_X) + 
  ylim(min_y, max_y) + 
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    # axis.title.x = element_text(size = 14),
    # axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  ) + 
  annotate("text", x = 3, y = -1.7, label = paste("R²:", round(rsq_rf, 4)),
           hjust = 1, vjust = 1, size = 4, color = "black", fontface = "bold")
rsq_gbm <- rsquare_summary %>% filter(Model == "gbm") %>% select(Mean) %>% as.numeric()
plot_gbm <- ggplot(plot_data, aes(x = response_scaled, y = pred_gbm)) +
  geom_point(color = ordered_models %>% filter(model == "gbm") %>% select(colour) %>% as.character(), alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(title = "Gradient Boosting Predictions",
       x = paste0("Actual ",toupper(response)),
       y = paste0("Predicted ",toupper(response))) +
  theme_minimal() + 
  xlim(min_X, max_X) + 
  ylim(min_y, max_y) + 
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    # axis.title.x = element_text(size = 14),
    # axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  ) + 
  annotate("text", x = 3, y = -1.7, label = paste("R²:", round(rsq_gbm, 4)),
           hjust = 1, vjust = 1, size = 4, color = "black", fontface = "bold")

# Plot residuals 
min_X <- floor(min(plot_data %>% select(starts_with("pred"))))
max_X <- ceiling(max(plot_data %>% select(starts_with("pred"))))
min_y <- floor(min(plot_data %>% select(starts_with("resid"))))
max_y <- ceiling(max(plot_data %>% select(starts_with("resid"))))
predict_lm_rmse <- ordered_models %>% filter(model == "lm") %>% select(RMSE) %>% as.numeric()
resid_plot_lm <- ggplot(plot_data, aes(x = pred_lm, y = resid_lm)) +
  geom_point(color = ordered_models %>% filter(model == "lm") %>% select(colour) %>% as.character(), alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Linear Residuals",
       x = paste0("Predicted ",toupper(response)),
       y = "Residuals") +
  theme_minimal() + 
  ylim(min_y, max_y) + 
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  ) + 
  annotate("text", x = 3, y = -2.5, label = paste("RMSE:", round(predict_lm_rmse, 4)),
           hjust = 1, vjust = 1, size = 4, color = "black", fontface = "bold")
predict_pm_rmse <- ordered_models %>% filter(model == "pm") %>% select(RMSE) %>% as.numeric()
resid_plot_poly <- ggplot(plot_data, aes(x = pred_poly, y = resid_poly)) +
  geom_point(color = ordered_models %>% filter(model == "pm") %>% select(colour) %>% as.character(), alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Polynomial Residuals",
       x = paste0("Predicted ",toupper(response)),
       y = "Residuals") +
  theme_minimal() + 
  ylim(min_y, max_y) + 
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    # axis.title.x = element_text(size = 14),
    # axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  ) + 
  annotate("text", x = 3, y = -2.5, label = paste("RMSE:", round(predict_pm_rmse, 4)),
           hjust = 1, vjust = 1, size = 4, color = "black", fontface = "bold")
predict_tm_rmse <- ordered_models %>% filter(model == "tm") %>% select(RMSE) %>% as.numeric()
resid_plot_tree <- ggplot(plot_data, aes(x = pred_tree, y = resid_tree)) +
  geom_point(color = ordered_models %>% filter(model == "tm") %>% select(colour) %>% as.character(), alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Regression Tree Residuals",
       x = paste0("Predicted ",toupper(response)),
       y = "Residuals") +
  theme_minimal() + 
  ylim(min_y, max_y) + 
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    # axis.title.x = element_text(size = 14),
    # axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  ) + 
  annotate("text", x = 3, y = -2.5, label = paste("RMSE:", round(predict_tm_rmse, 4)),
           hjust = 1, vjust = 1, size = 4, color = "black", fontface = "bold")
predict_rf_rmse <- ordered_models %>% filter(model == "rfm") %>% select(RMSE) %>% as.numeric()
resid_plot_rf <- ggplot(plot_data, aes(x = pred_rf, y = resid_rf)) +
  geom_point(color = ordered_models %>% filter(model == "rfm") %>% select(colour) %>% as.character(), alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Random Forest Residuals",
       x = paste0("Predicted ",toupper(response)),
       y = "Residuals") +
  theme_minimal() + 
  ylim(min_y, max_y) + 
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    # axis.title.x = element_text(size = 14),
    # axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  ) + 
  annotate("text", x = 3, y = -2.5, label = paste("RMSE:", round(predict_rf_rmse, 4)),
           hjust = 1, vjust = 1, size = 4, color = "black", fontface = "bold")
predict_gbm_rmse <- ordered_models %>% filter(model == "gbm") %>% select(RMSE) %>% as.numeric()
resid_plot_gbm <- ggplot(plot_data, aes(x = pred_gbm, y = resid_gbm)) +
  geom_point(color = ordered_models %>% filter(model == "gbm") %>% select(colour) %>% as.character(), alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Gradient Boosting Residuals",
       x = paste0("Predicted ",toupper(response)),
       y = "Residuals") +
  theme_minimal() + 
  ylim(min_y, max_y) + 
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    # axis.title.x = element_text(size = 14),
    # axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  ) + 
  annotate("text", x = 3, y = -2.5, label = paste("RMSE:", round(predict_gbm_rmse, 4)),
           hjust = 1, vjust = 1, size = 4, color = "black", fontface = "bold")


# Plot distribution of residuals
min_X <- floor(min(plot_data %>% select(starts_with("resid"))))
max_X <- ceiling(max(plot_data %>% select(starts_with("resid"))))
resid_dist_lm <- ggplot(plot_data, aes(x = resid_lm)) +
  geom_histogram(fill = ordered_models %>% filter(model == "lm") %>% select(colour) %>% as.character(), bins = 30, alpha = 0.5) +
  labs(title = "Linear Residuals Distribution",
       x = "Residuals",
       y = "Frequency") +
  theme_minimal() + 
  xlim(min_X, max_X) + 
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  )
resid_dist_poly <- ggplot(plot_data, aes(x = resid_poly)) +
  geom_histogram(fill = ordered_models %>% filter(model == "pm") %>% select(colour) %>% as.character(), bins = 30, alpha = 0.5) +
  labs(title = "Polynomial Residuals Distribution",
       x = "Residuals",
       y = "Frequency") +
  theme_minimal() + 
  xlim(min_X, max_X) + 
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    # axis.title.x = element_text(size = 14),
    # axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  )
resid_dist_tree <- ggplot(plot_data, aes(x = resid_tree)) +
  geom_histogram(fill = ordered_models %>% filter(model == "tm") %>% select(colour) %>% as.character(), bins = 30, alpha = 0.5) +
  labs(title = "Regression Tree Residuals Distribution",
       x = "Residuals",
       y = "Frequency") +
  theme_minimal() + 
  xlim(min_X, max_X) + 
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    # axis.title.x = element_text(size = 14),
    # axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  )
resid_dist_rf <- ggplot(plot_data, aes(x = resid_rf)) +
  geom_histogram(fill = ordered_models %>% filter(model == "rfm") %>% select(colour) %>% as.character(), bins = 30, alpha = 0.5) +
  labs(title = "Random Forest Residuals Distribution",
       x = "Residuals",
       y = "Frequency") +
  theme_minimal() + 
  xlim(min_X, max_X) + 
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    # axis.title.x = element_text(size = 14),
    # axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  )
resid_dist_gbm <- ggplot(plot_data, aes(x = resid_gbm)) +
  geom_histogram(fill = ordered_models %>% filter(model == "gbm") %>% select(colour) %>% as.character(), bins = 30, alpha = 0.5) +
  labs(title = "Gradient Boosting Residuals Distribution",
       x = "Residuals",
       y = "Frequency") +
  theme_minimal() + 
  xlim(min_X, max_X) + 
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    # axis.title.x = element_text(size = 14),
    # axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  )

lm_label <- textGrob("Linear", gp = gpar(fontsize = 16), vjust = 1)
pm_label <- textGrob("Polynomial", gp = gpar(fontsize = 16), vjust = 1)
tm_label <- textGrob("Regression Tree", gp = gpar(fontsize = 16), vjust = 1)
rf_label <- textGrob("Random Forest", gp = gpar(fontsize = 16), vjust = 1)
gbm_label <- textGrob("Gradient Boosting", gp = gpar(fontsize = 16), vjust = 1)
acc_predict_label <- textGrob(paste0("Actual vs Predicted"), rot = 270, gp = gpar(fontsize = 16))
residual_fit_label <- textGrob("Predicted vs Residuals", rot = 270, gp = gpar(fontsize = 16))
residual_dist_label <- textGrob("Residuals vs Frequency", rot = 270, gp = gpar(fontsize = 16))
main_title <-  textGrob(paste0("Model prediction of GWAS ",toupper(response)), gp = gpar(fontsize = 18, fontface = "bold"))

predict_plots <- grid.arrange(
  arrangeGrob(
    lm_label,textGrob(""), pm_label, tm_label, gbm_label, rf_label, textGrob(""), 
    ncol = 7, widths = c(1.1, 0.05, 1, 1, 1, 1, 0.1)
  ),
  arrangeGrob(
    plot_lm + theme(axis.title = element_blank(), plot.title = element_blank()), 
    textGrob(""),
    plot_poly + theme(axis.title = element_blank(), plot.title = element_blank()), 
    plot_tree + theme(axis.title = element_blank(), plot.title = element_blank()), 
    plot_gbm + theme(axis.title = element_blank(), plot.title = element_blank()), 
    plot_rf + theme(axis.title = element_blank(), plot.title = element_blank()),
    acc_predict_label,
    ncol = 7, widths = c(1.1, 0.05, 1, 1, 1, 1, 0.1)
  ),
  arrangeGrob(
    resid_plot_lm + theme(axis.title = element_blank(), plot.title = element_blank()), 
    textGrob(""),
    resid_plot_poly + theme(axis.title = element_blank(), plot.title = element_blank()), 
    resid_plot_tree + theme(axis.title = element_blank(), plot.title = element_blank()), 
    resid_plot_gbm + theme(axis.title = element_blank(), plot.title = element_blank()), 
    resid_plot_rf + theme(axis.title = element_blank(), plot.title = element_blank()), 
    residual_fit_label,
    ncol = 7, widths = c(1.1, 0.05, 1, 1, 1, 1, 0.1)
  ),
  arrangeGrob(
    resid_dist_lm + theme(axis.title = element_blank(), plot.title = element_blank()), 
    textGrob(""),
    resid_dist_poly + theme(axis.title = element_blank(), plot.title = element_blank()), 
    resid_dist_tree + theme(axis.title = element_blank(), plot.title = element_blank()), 
    resid_dist_gbm + theme(axis.title = element_blank(), plot.title = element_blank()), 
    resid_dist_rf + theme(axis.title = element_blank(), plot.title = element_blank()), 
    residual_dist_label,
    ncol = 7, widths = c(1.1, 0.05, 1, 1, 1, 1, 0.1)
  ),
  ncol = 1,
  heights = c(0.1, 1, 1, 1),
  top = main_title
)
ggsave(filename = paste0(directory,"/",response,"_model_predict.png"), plot = predict_plots, device = "png", bg = 'white')


# Test the best performing model minus some of the predictors

# Define the Random Forest training function
train_rf <- function(data, response, predictors) {
  set.seed(1)
  # Random forest model
  rf_grid <- expand.grid(
    mtry = c(2, 4, 6)
  )
  model <- train(
    as.formula(paste(response, "~", paste(predictors, collapse = "+"))),
    data = data,
    method = "rf",
    trControl = train_control,
    tuneGrid = rf_grid
  )
  return(model)
}
# Define the GB training function
train_gb <- function(data, response, predictors) {
  set.seed(1)
  # Gradient boosting model forest model
  gb_grid <- expand.grid(
    interaction.depth = c(1,2,3,4),
    n.trees = c(50, 100, 200, 400, 600),
    shrinkage = c(0.01,0.1),
    n.minobsinnode = 10
  )
  model <- train(
    as.formula(paste(response, "~", paste(predictors, collapse = "+"))),
    data = data,
    method = "gbm",
    trControl = train_control,
    tuneGrid = gb_grid,
    distribution = "laplace",
    verbose = FALSE
  )
}

standardise_col <- function(df, column_name) {
  total_sum <- sum(df[[column_name]], na.rm = TRUE)
  df[[column_name]] <- df[[column_name]] / total_sum
  df[is.na(df[[column_name]]), column_name] <- 0
  return(df)
}
# Get the full model stats
all_predictors <- c("n_individuals_scaled", 
                    "n_snps_scaled", 
                    "causal_snps_scaled", 
                    "cc_ratio_scaled", 
                    "causal_maf_scaled", 
                    "heritability_scaled")

# Get the weights for the full RF model
rf_full <- rf_cv
varImp_rf <- as.data.frame(varImp(rf_full)$importance)
rf_models_list <- list(full = rf_full)
varImp_rf$Predictor <- rownames(varImp_rf)
colnames(varImp_rf)[1] <- "full_scaled"
varImp_rf <- standardise_col(varImp_rf, "full_scaled")
varImp_rf_list <- list(rf_full = varImp_rf)
# Get the weights for the full GB model
gbm_full <- gbm_cv
gb_models_list <- list(full = gbm_full)
varImp_gb <- varImp_gbm %>%
  select("full_scaled" = gbm_importance, Predictor)
varImp_gb <- standardise_col(varImp_gb, "full_scaled")
varImp_gb_list <- list(gbm_full = varImp_gb)

# Train models excluding each predictor
for (predictor in all_predictors) {
  print(predictor)
  predictors_subset <- setdiff(all_predictors, predictor)
  model_name <- paste("no", predictor, sep = "_")
  
  model_rf <- train_rf(data_train, "response_scaled", predictors_subset)
  model_gb <- train_gb(data_train, "response_scaled", predictors_subset)
  
  rf_models_list[[model_name]]  <- model_rf
  gb_models_list[[model_name]] <- model_gb
  
  rf_varImp_tmp <- as.data.frame(varImp(model_rf)$importance)
  rf_varImp_tmp$Predictor <- rownames(rf_varImp_tmp)
  colnames(rf_varImp_tmp)[1] <- paste("no", predictor, sep = "_")
  rf_varImp_tmp <- standardise_col(rf_varImp_tmp, colnames(rf_varImp_tmp)[1])
  varImp_rf_list[[model_name]] <- rf_varImp_tmp
  
  gb_varImp_tmp  <- summary(model_gb$finalModel)
  gb_varImp_tmp <- as.data.frame(gb_varImp_tmp) %>%
    select(rel.inf, "Predictor" = var)
  colnames(gb_varImp_tmp)[1] <- paste("no", predictor, sep = "_")
  gb_varImp_tmp <- standardise_col(gb_varImp_tmp, colnames(gb_varImp_tmp)[1])
  varImp_gb_list[[model_name]] <- gb_varImp_tmp
}

# Perform resampling on all models
rf_resamples_list <- resamples(rf_models_list)
rf_resample_summary <- summary(rf_resamples_list)
gb_resamples_list <- resamples(gb_models_list)
gb_resample_summary <- summary(gb_resamples_list)

# Check for normality of differences
rf_noindvs_rmse <- rf_resamples_list$values$`no_n_individuals_scaled~RMSE`
rf_nosnps_rmse <- rf_resamples_list$values$`no_n_snps_scaled~RMSE`
rf_nocausnps_rmse <- rf_resamples_list$values$`no_causal_snps_scaled~RMSE`
rf_noccratio_rmse <- rf_resamples_list$values$`no_cc_ratio_scaled~RMSE`
rf_nomaf_rmse <- rf_resamples_list$values$`no_causal_maf_scaled~RMSE`
rf_noherit_rmse <- rf_resamples_list$values$`no_heritability_scaled~RMSE`

gb_noindvs_rmse <- gb_resamples_list$values$`no_n_individuals_scaled~RMSE`
gb_nosnps_rmse <- gb_resamples_list$values$`no_n_snps_scaled~RMSE`
gb_nocausnps_rmse <- gb_resamples_list$values$`no_causal_snps_scaled~RMSE`
gb_noccratio_rmse <- gb_resamples_list$values$`no_cc_ratio_scaled~RMSE`
gb_nomaf_rmse <- gb_resamples_list$values$`no_causal_maf_scaled~RMSE`
gb_noherit_rmse <- gb_resamples_list$values$`no_heritability_scaled~RMSE`

rf_rmse_diffs <- list(
  noindvs_nosnps = rf_noindvs_rmse - rf_nosnps_rmse,
  noindvs_nocausnps = rf_noindvs_rmse - rf_nocausnps_rmse,
  noindvs_noccratio = rf_noindvs_rmse - rf_noccratio_rmse,
  noindvs_nomaf = rf_noindvs_rmse - rf_nomaf_rmse,
  noindvs_noherit = rf_noindvs_rmse - rf_noherit_rmse,
  nosnps_nocausnps = rf_nosnps_rmse - rf_nocausnps_rmse,
  nosnps_noccratio = rf_nosnps_rmse - rf_noccratio_rmse,
  nosnps_nomaf = rf_nosnps_rmse - rf_nomaf_rmse,
  nosnps_noherit = rf_nosnps_rmse - rf_noherit_rmse,
  nocausnps_noccratio = rf_nocausnps_rmse - rf_noccratio_rmse,
  nocausnps_nomaf = rf_nocausnps_rmse - rf_nomaf_rmse,
  nocausnps_noherit = rf_nocausnps_rmse - rf_noherit_rmse,
  noccratio_nomaf = rf_noccratio_rmse - rf_nomaf_rmse,
  noccratio_noherit = rf_noccratio_rmse - rf_noherit_rmse,
  nomaf_noherit = rf_nomaf_rmse - rf_noherit_rmse
)

gb_rmse_diffs <- list(
  noindvs_nosnps = gb_noindvs_rmse - gb_nosnps_rmse,
  noindvs_nocausnps = gb_noindvs_rmse - gb_nocausnps_rmse,
  noindvs_noccratio = gb_noindvs_rmse - gb_noccratio_rmse,
  noindvs_nomaf = gb_noindvs_rmse - gb_nomaf_rmse,
  noindvs_noherit = gb_noindvs_rmse - gb_noherit_rmse,
  nosnps_nocausnps = gb_nosnps_rmse - gb_nocausnps_rmse,
  nosnps_noccratio = gb_nosnps_rmse - gb_noccratio_rmse,
  nosnps_nomaf = gb_nosnps_rmse - gb_nomaf_rmse,
  nosnps_noherit = gb_nosnps_rmse - gb_noherit_rmse,
  nocausnps_noccratio = gb_nocausnps_rmse - gb_noccratio_rmse,
  nocausnps_nomaf = gb_nocausnps_rmse - gb_nomaf_rmse,
  nocausnps_noherit = gb_nocausnps_rmse - gb_noherit_rmse,
  noccratio_nomaf = gb_noccratio_rmse - gb_nomaf_rmse,
  noccratio_noherit = gb_noccratio_rmse - gb_noherit_rmse,
  nomaf_noherit = gb_nomaf_rmse - gb_noherit_rmse
)

rf_rmse_normality <- lapply(rf_rmse_diffs, check_normality)
for(test in names(rf_rmse_normality)){
  if(rf_rmse_normality[[test]]$p.value < 0.05){
    print(paste0("Significant difference in ",test, " p value = ",rf_rmse_normality[[test]]$p.value))
  }
}
gb_rmse_normality <- lapply(gb_rmse_diffs, check_normality)
for(test in names(gb_rmse_normality)){
  if(gb_rmse_normality[[test]]$p.value < 0.05){
    print(paste0("Significant difference in ",test, " p value = ",gb_rmse_normality[[test]]$p.value))
  }
}

rf_noindvs_rsqr <- rf_resamples_list$values$`no_n_individuals_scaled~Rsquared`
rf_nosnps_rsqr <- rf_resamples_list$values$`no_n_snps_scaled~Rsquared`
rf_nocausnps_rsqr <- rf_resamples_list$values$`no_causal_snps_scaled~Rsquared`
rf_noccratio_rsqr <- rf_resamples_list$values$`no_cc_ratio_scaled~Rsquared`
rf_nomaf_rsqr <- rf_resamples_list$values$`no_causal_maf_scaled~Rsquared`
rf_noherit_rsqr <- rf_resamples_list$values$`no_heritability_scaled~Rsquared`

gb_noindvs_rsqr <- gb_resamples_list$values$`no_n_individuals_scaled~Rsquared`
gb_nosnps_rsqr <- gb_resamples_list$values$`no_n_snps_scaled~Rsquared`
gb_nocausnps_rsqr <- gb_resamples_list$values$`no_causal_snps_scaled~Rsquared`
gb_noccratio_rsqr <- gb_resamples_list$values$`no_cc_ratio_scaled~Rsquared`
gb_nomaf_rsqr <- gb_resamples_list$values$`no_causal_maf_scaled~Rsquared`
gb_noherit_rsqr <- gb_resamples_list$values$`no_heritability_scaled~Rsquared`

rf_rsqur_diffs <- list(
  noindvs_nosnps = rf_noindvs_rsqr - rf_nosnps_rsqr,
  noindvs_nocausnps = rf_noindvs_rsqr - rf_nocausnps_rsqr,
  noindvs_noccratio = rf_noindvs_rsqr - rf_noccratio_rsqr,
  noindvs_nomaf = rf_noindvs_rsqr - rf_nomaf_rsqr,
  noindvs_noherit = rf_noindvs_rsqr - rf_noherit_rsqr,
  nosnps_nocausnps = rf_nosnps_rsqr - rf_nocausnps_rsqr,
  nosnps_noccratio = rf_nosnps_rsqr - rf_noccratio_rsqr,
  nosnps_nomaf = rf_nosnps_rsqr - rf_nomaf_rsqr,
  nosnps_noherit = rf_nosnps_rsqr - rf_noherit_rsqr,
  nocausnps_noccratio = rf_nocausnps_rsqr - rf_noccratio_rsqr,
  nocausnps_nomaf = rf_nocausnps_rsqr - rf_nomaf_rsqr,
  nocausnps_noherit = rf_nocausnps_rsqr - rf_noherit_rsqr,
  noccratio_nomaf = rf_noccratio_rsqr - rf_nomaf_rsqr,
  noccratio_noherit = rf_noccratio_rsqr - rf_noherit_rsqr,
  nomaf_noherit = rf_nomaf_rsqr - rf_noherit_rsqr
)

gb_rsqur_diffs <- list(
  noindvs_nosnps = gb_noindvs_rsqr - gb_nosnps_rsqr,
  noindvs_nocausnps = gb_noindvs_rsqr - gb_nocausnps_rsqr,
  noindvs_noccratio = gb_noindvs_rsqr - gb_noccratio_rsqr,
  noindvs_nomaf = gb_noindvs_rsqr - gb_nomaf_rsqr,
  noindvs_noherit = gb_noindvs_rsqr - gb_noherit_rsqr,
  nosnps_nocausnps = gb_nosnps_rsqr - gb_nocausnps_rsqr,
  nosnps_noccratio = gb_nosnps_rsqr - gb_noccratio_rsqr,
  nosnps_nomaf = gb_nosnps_rsqr - gb_nomaf_rsqr,
  nosnps_noherit = gb_nosnps_rsqr - gb_noherit_rsqr,
  nocausnps_noccratio = gb_nocausnps_rsqr - gb_noccratio_rsqr,
  nocausnps_nomaf = gb_nocausnps_rsqr - gb_nomaf_rsqr,
  nocausnps_noherit = gb_nocausnps_rsqr - gb_noherit_rsqr,
  noccratio_nomaf = gb_noccratio_rsqr - gb_nomaf_rsqr,
  noccratio_noherit = gb_noccratio_rsqr - gb_noherit_rsqr,
  nomaf_noherit = gb_nomaf_rsqr - gb_noherit_rsqr
)

rf_rsqur_normality <- lapply(rf_rsqur_diffs, check_normality)
for(test in names(rf_rsqur_normality)){
  if(rf_rsqur_normality[[test]]$p.value < 0.05){
    print(paste0("Significant difference in ",test, " p value = ",rf_rsqur_normality[[test]]$p.value))
  }
}
gb_rsqur_normality <- lapply(gb_rsqur_diffs, check_normality)
for(test in names(gb_rsqur_normality)){
  if(gb_rsqur_normality[[test]]$p.value < 0.05){
    print(paste0("Significant difference in ",test, " p value = ",gb_rsqur_normality[[test]]$p.value))
  }
}

# Perform pairwise t-test comparisons of resampling results
rf_diffs <- diff(rf_resamples_list)
rf_significant_differences <- summary(rf_diffs, adjust = "bonferroni")
gb_diffs <- diff(gb_resamples_list)
gb_significant_differences <- summary(gb_diffs, adjust = "bonferroni")

# Extract RMSE and R-squared values and calculate mean and 95% CI
rf_summary_stats <- list()
gb_summary_stats <- list()
for(model in 1:2){
  if(model == 1){
    resamples_list <- rf_resamples_list
    significant_differences <- rf_significant_differences
  } else {
    resamples_list <- gb_resamples_list
    significant_differences <- gb_significant_differences
  }
  for(measure in c("RMSE","Rsquared")){
    # Get confidence intervals
    model_summary <- resamples_list$values %>%
      select(contains(paste0(measure))) %>%
      pivot_longer(everything(), names_to = "Model", values_to = "Value") %>%
      mutate(Model = gsub(paste0("~",measure), "", Model)) %>%
      group_by(Model) %>%
      summarize(
        Mean = mean(Value),
        Lower_CI = mean(Value) - qt(1 - 0.05 / 2, n() - 1) * sd(Value) / sqrt(n()),
        Upper_CI = mean(Value) + qt(1 - 0.05 / 2, n() - 1) * sd(Value) / sqrt(n())
      ) %>% {
        if(measure == "RMSE"){
          arrange(.,Mean) %>%
            mutate(rank = row_number())
        } else {
          arrange(.,desc(Mean)) %>%
            mutate(rank = row_number())
        }
      }
    model_summary$Model <- factor(model_summary$Model, levels = model_summary$Model)
    summary_stats[[paste0(measure)]] <- model_summary
    model_rank <- model_summary %>% select(Model,Upper_CI,rank)
    # Get significance values 
    model_matrix <- significant_differences$table[[paste0(measure)]]
    model_matrix[upper.tri(model_matrix)] <- ""
    signif_comparisons <- model_matrix %>%
      as.data.frame() %>%
      rownames_to_column(var = "Model1") %>%
      pivot_longer(cols = -Model1, names_to = "Model2", values_to = "p_value") %>%
      filter(p_value != "") %>%
      mutate(p_value = as.numeric(p_value)) %>%
      mutate(significance = case_when(
        p_value < 0.0001 ~ "****",
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*"
      )) %>%
      filter(!if_any(everything(), is.na)) %>%
      select(Model1, Model2, p_value, significance) %>%
      left_join(model_rank, by = c("Model1" = "Model")) %>%
      left_join(model_rank, by = c("Model2" = "Model")) %>%
      mutate(priCI = ifelse(Upper_CI.x > Upper_CI.y, Upper_CI.x, Upper_CI.y),
             secCI = ifelse(Upper_CI.x < Upper_CI.y, Upper_CI.x, Upper_CI.y)) %>%
      arrange(priCI) %>%
      group_by(secCI) %>%
      mutate(secrank = cur_group_id()) %>%
      ungroup() %>%
      arrange(priCI,secCI) 
    # Get coordinates for sig lines
    comparisons_list <- list()
    annotations <- c()
    y_positions <- c()
    base_y_position <- max(model_summary$Upper_CI) * 1
    max_ci_track <- c()
    prev_pos <- NULL
    for (i in 1:nrow(signif_comparisons)) {
      comparisons_list[[i]] <- c(signif_comparisons$Model1[i], signif_comparisons$Model2[i])
      annotations[i] <- signif_comparisons$significance[i]
      priCI <- signif_comparisons[i,]$priCI
      secrank <- signif_comparisons[i,]$secrank
      predict_pos <- priCI + 0.04 * secrank
      if(secrank != 1){
        if(prev_pos > priCI){
          predict_pos <- prev_pos + 0.08
        }
      } else {
        if(i != 1){
          if(prev_pos + 0.06 > priCI){
            predict_pos <- prev_pos + 0.08
          } else {
            predict_pos <- predict_pos - 0.02 * secrank # Initial elevation above the plot
          }
        } else {
          predict_pos <- predict_pos - 0.02 * secrank # Initial elevation above the plot
        }
      }
      prev_pos <- predict_pos
      
      y_positions[i] <- predict_pos
    }
    if(model == 1){
      rf_summary_stats[[paste0(measure)]] <- list(model_summary = model_summary,
                                                  model_pvals = signif_comparisons,
                                                  comparisons_list = comparisons_list,
                                                  annotations = annotations,
                                                  y_positions = y_positions)
    } else {
      gb_summary_stats[[paste0(measure)]] <- list(model_summary = model_summary,
                                                  model_pvals = signif_comparisons,
                                                  comparisons_list = comparisons_list,
                                                  annotations = annotations,
                                                  y_positions = y_positions)
    }
  }
}


# Plot RMSE with error bars and significance annotations
model_ggplots <- list()
for(model in c("RF","GB")){
  if(model == "RF"){
    rmse_summary <- rf_summary_stats$RMSE$model_summary
    summary_stats <- rf_summary_stats
    varImp_list <- varImp_rf_list
  } else {
    rmse_summary <- gb_summary_stats$RMSE$model_summary
    summary_stats <- gb_summary_stats
    varImp_list <- varImp_gb_list
  }
  rmse_lables <- rmse_summary %>% 
    mutate(Model = case_when(
      str_detect(Model, "full") ~ "Full Model",
      str_detect(Model, "n_snps") ~ "N.SNPs",
      str_detect(Model, "causal_maf") ~ "Causal MAF",
      str_detect(Model, "n_individual") ~ "N.Individuals",
      str_detect(Model, "cc_ratio") ~ "Case/Control Ratio",
      str_detect(Model, "heritability") ~ "Heritability",
      str_detect(Model, "causal_snps") ~ "N.Causal SNPs",
      TRUE ~ Model
    )) %>% select(Model) %>% as.matrix() %>% as.character()
  y_breaks <- pretty(c(0, max(rmse_summary$Upper_CI) * 1.01), n = 5)
  y_labels <- c(0,0.2,0.4,0.6,0.8,1.0,1.2,1.4)[1:length(y_breaks)]
  rmse_sig_plot <- ggplot(rmse_summary, aes(x = Model, y = Mean)) +
    geom_bar(stat = "identity", fill = "grey", color = "black", width = 0.7) +
    geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.2) +
    labs(title = "Missing Predictor Residual Error", y = "Mean RMSE", x = "Missing Predictor") +
    theme_minimal() +
    geom_signif(
      comparisons = summary_stats$RMSE$comparisons_list,
      annotations = summary_stats$RMSE$annotations,
      y_position = summary_stats$RMSE$y_positions,
      tip_length = 0.03,
      textsize = 4,
      vjust = 0.4
    ) + 
    scale_y_continuous(
      breaks = y_breaks, 
      labels = y_labels,
      expand = expansion(mult = c(0, 0.2))
    ) +
    scale_x_discrete(labels = rmse_lables) + 
    theme(
      plot.margin = unit(c(3, 1, 1, 1), "lines"),
      panel.grid.major.y = element_line(color = "grey", size = 0.5),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.title.y = element_text(hjust = 0.25, vjust = 2.5, size = 14),
      axis.title.x = element_text(size = 14),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12)
    )
  
  # Plot Rsquared with error bars and significance annotations
  rsquare_summary <- summary_stats$Rsquared$model_summary
  rsqr_lables <- rmse_lables
  y_breaks <- pretty(c(0, max(rsquare_summary$Upper_CI) * 1.01), n = 5)
  y_labels <- c(0,0.2,0.4,0.6,0.8,1.0,1.2,1.4)[1:length(y_breaks)]
  rsquare_sig_plot <- ggplot(rsquare_summary, aes(x = Model, y = Mean)) + 
    geom_bar(stat = "identity", fill = "grey", color = "black", width = 0.7) +
    geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.2) +
    labs(title = "Missing Predictor Varience Explained", y = "Mean Rsquared", x = "Missing Predictor") +
    theme_minimal() +
    geom_signif(
      comparisons = summary_stats$Rsquared$comparisons_list,
      annotations = summary_stats$Rsquared$annotations,
      y_position = summary_stats$Rsquared$y_positions,
      tip_length = 0.03,
      textsize = 4,
      vjust = 0.4
    ) + 
    scale_y_continuous(
      breaks = y_breaks, 
      labels = y_labels,
      expand = expansion(mult = c(0, 0.2))
    ) +
    scale_x_discrete(labels = rsqr_lables) + 
    theme(
      plot.margin = unit(c(3, 1, 1, 1), "lines"),
      panel.grid.major.y = element_line(color = "grey", size = 0.5),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.title.y = element_text(hjust = 0.25,vjust = 2.5, size = 14),
      axis.title.x = element_text(size = 14),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12)
    )
  
  importance_df <- reduce(varImp_list, function(x, y) merge(x, y, by = "Predictor", all = TRUE)) %>%
    mutate(Predictor = gsub("_scaled", "", Predictor))
  importance_df[is.na(importance_df)] <- 0
  
  importance_long <- importance_df %>% 
    mutate(total_imp = rowSums(select(., ends_with("_scaled")))) %>%
    arrange(desc(total_imp)) %>%
    select(-total_imp) %>%
    pivot_longer(., cols = ends_with("_scaled"), names_to = "Model", values_to = "Importance") %>%
    mutate(Model = case_when(
      str_detect(Model, "full") ~ "Full Model",
      str_detect(Model, "n_snps") ~ "N.SNPs",
      str_detect(Model, "causal_maf") ~ "Causal MAF",
      str_detect(Model, "n_individual") ~ "N.Individuals",
      str_detect(Model, "cc_ratio") ~ "Case/Control Ratio",
      str_detect(Model, "heritability") ~ "Heritability",
      str_detect(Model, "causal_snps") ~ "N.Causal SNPs",
      TRUE ~ Model
    ))
  
  heatmap_rmse <- rmse_summary %>% 
    select(Model,"Importance" = Mean, "RMSE_rank" = rank) %>%
    mutate(Predictor = "RMSE",
           Model = case_when(
             str_detect(Model, "full") ~ "Full Model",
             str_detect(Model, "n_snps") ~ "N.SNPs",
             str_detect(Model, "causal_maf") ~ "Causal MAF",
             str_detect(Model, "n_individual") ~ "N.Individuals",
             str_detect(Model, "cc_ratio") ~ "Case/Control Ratio",
             str_detect(Model, "heritability") ~ "Heritability",
             str_detect(Model, "causal_snps") ~ "N.Causal SNPs",
             TRUE ~ Model
           )) %>%
    bind_rows(importance_long) %>%
    mutate(Type = if_else(Predictor == "RMSE", "RMSE", "Importance"))
  
  heatmap_rmse_rsq <- rsquare_summary %>% 
    select(Model,"Importance" = Mean, "RSQ_rank" = rank) %>%
    mutate(Predictor = "RSQ",
           Model = case_when(
             str_detect(Model, "full") ~ "Full Model",
             str_detect(Model, "n_snps") ~ "N.SNPs",
             str_detect(Model, "causal_maf") ~ "Causal MAF",
             str_detect(Model, "n_individual") ~ "N.Individuals",
             str_detect(Model, "cc_ratio") ~ "Case/Control Ratio",
             str_detect(Model, "heritability") ~ "Heritability",
             str_detect(Model, "causal_snps") ~ "N.Causal SNPs",
             TRUE ~ Model
           )) %>%
    bind_rows(heatmap_rmse) %>%
    mutate(Type = if_else(Predictor == "RSQ", "RSQ", Type)) 
  vertical_line_pos <- max(which(heatmap_rmse_rsq$Predictor != "RMSE" & 
                                   heatmap_rmse_rsq$Predictor != "RSQ"))
  
  rmse_order <- heatmap_rmse %>%
    arrange(desc(RMSE_rank)) %>%
    pull(Model) %>%
    unique()
  
  heatmap_rmse_rsq$Model <- factor(heatmap_rmse_rsq$Model, 
                                   levels = rmse_order)
  # Ensure the column factor levels are ordered by importance
  col_order <- c(unique(importance_long$Predictor),"RMSE","RSQ")
  heatmap_rmse_rsq$Predictor <- factor(heatmap_rmse_rsq$Predictor, 
                                       levels = col_order)
  
  # Plot the heatmap with the RMSE column
  importance_heat <- ggplot() +
    geom_tile(data = heatmap_rmse_rsq[heatmap_rmse_rsq$Type == "Importance",], aes(x = Predictor, y = Model, fill = Importance)) +
    scale_fill_gradient(low = "white", high = "darkred", name = "Importance") +
    new_scale("fill") + 
    geom_tile(data = heatmap_rmse_rsq[heatmap_rmse_rsq$Type == "RMSE",], aes(x = Predictor, y = Model, fill = Importance)) +
    scale_fill_gradient(low = "darkblue", high = "lightblue", name = "RMSE") +
    new_scale("fill") + 
    geom_tile(data = heatmap_rmse_rsq[heatmap_rmse_rsq$Type == "RSQ",], aes(x = Predictor, y = Model, fill = Importance)) +
    scale_fill_gradient(low = "lightgreen", high = "darkgreen", name = "RSQ") +
    theme_minimal() +
    geom_vline(xintercept = 6.55, color = "white", linewidth = 6) +
    geom_vline(xintercept = 7.5, color = "white", linewidth = 6) +
    geom_vline(xintercept = 8.45, color = "white", linewidth = 6) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 14, hjust = 0.45),
          axis.title.y = element_text(size = 14),
          legend.title = element_text(size = 12),
          plot.title = element_text(hjust = 0.5, size = 16, face = "bold")) +
    labs(title = "Random Forest Predictor Importance, RMSE, and RSQ",
         x = "Predictor",
         y = "Missing Predictor Model") + 
    scale_x_discrete(drop = FALSE)
  
  create_pie_chart <- function(data, model_name, title_size = 14, show_legend = TRUE, legend_text_size = 12) {
    data_filtered <- data %>% 
      filter(Model == model_name & Type == "Importance") %>%
      arrange(Importance) %>%
      mutate(
        Importance = as.numeric(Importance),
        csum = cumsum(Importance),
        ypos = csum - 0.5 * Importance,
        label = ifelse(Importance > 0.05, round(Importance, 2),""),
        angle = 90 - 360 * (csum - 0.5 * Importance) / sum(Importance),
        angle = ifelse(angle < -90, angle + 180, angle)
      )
    
    p <- ggplot(data_filtered, aes(x = "", y = Importance, fill = Predictor)) +
      geom_bar(width = 1, stat = "identity") +
      coord_polar("y", start = 0) +
      theme_void() +
      labs(title = paste0("-",model_name)) +
      theme(plot.title = element_text(size = title_size, vjust = -2, hjust = 0.2)) + 
      scale_fill_brewer(palette = "Set1")
    
    if (!show_legend) {
      p <- p + theme(legend.position = "none")
    } else {
      p <- p + theme(legend.position = "left",
                     legend.text = element_text(size = legend_text_size),
                     legend.title = element_text(size = legend_text_size + 2))
    }
    
    p <- p + geom_text(aes(y = ypos, label = label, angle = angle), size = 6)
    
    return(p)
  }
  
  # List of models including the full model and models with missing predictors
  models <- c("Full Model", "N.SNPs", "Causal MAF", "N.Individuals", "Case/Control Ratio", "Heritability", "N.Causal SNPs")
  
  # Create pie charts for all models without legends
  pie_charts <- lapply(models, function(model) create_pie_chart(heatmap_rmse_rsq, model, 
                                                                title_size = 14, 
                                                                show_legend = FALSE))
  names(pie_charts) <- models
  
  # Create a larger pie chart for the full model with legend
  data_filtered <- heatmap_rmse_rsq %>% 
    filter(Model == "Full Model" & Type == "Importance") %>%
    arrange(Importance) %>%
    mutate(
      Importance = as.numeric(Importance),
      csum = cumsum(Importance),
      ypos = csum - 0.5 * Importance,
      label = paste0(Predictor, ": ", round(Importance, 2))
    )
  
  # Create the pie chart with conditional labeling
  full_model_pie <- ggplot(data_filtered, aes(x = "", y = Importance, fill = Predictor)) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start = 0) +
    theme_void() +
    labs(title = "Full Model") +
    theme(plot.title = element_text(size = 18, face = "bold", vjust = -20, hjust = 0.1),
          legend.position = "left",
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16)) +
    scale_fill_brewer(palette = "Set1") + 
    geom_text(aes(y = ypos, label = ifelse(Importance > 0.2, label, "")), size = 6) +  # Display labels inside for larger slices
    geom_label_repel(aes(y = ypos, label = ifelse(Importance <= 0.2, label, "")),  # Repel labels for smaller slices
                     size = 6, nudge_x = 0.9, nudge_y = 0.02, box.padding = 0.5,
                     show.legend = FALSE,
                     force = 80,
                     force_pull = 1.5,
                     direction = "both",
                     seed = 123) + 
    theme(legend.margin = margin(t = 0, r = -40, b = 0, l = 20))
  
  
  
  # Arrange the pie charts in a custom layout
  plot_L<- plot_grid(NULL,
                     NULL,
                     NULL,
                     pie_charts[["N.SNPs"]],
                     NULL,
                     NULL,
                     NULL,
                     pie_charts[["Causal MAF"]],
                     pie_charts[["N.Causal SNPs"]],
                     pie_charts[["Heritability"]],
                     pie_charts[["Case/Control Ratio"]],
                     pie_charts[["N.Individuals"]],
                     ncol = 4)
  
  if(model == "RF"){
    title <- "Random Forest"
  } else {
    title <- "Gradient Boosting"
  }
  fnal_pie_charts <- ggdraw() + 
    draw_text(paste0(title," Predictor Weights"),
              x = 0.2, 
              y = 0.98,
              size = 16,
              fontface = "bold") + 
    draw_plot(plot_L) + 
    draw_plot(full_model_pie,
              x = 0.77, 
              y = 1.11, 
              hjust = 1, 
              vjust = 1, 
              halign = 1, 
              valign = 1,
              width = 0.78)
  
  
  model_ggplots[[model]] <- list(rmse_plot = rmse_sig_plot,
                                 rsqr_plot = rsquare_sig_plot,
                                 importance_heat_plot = importance_heat,
                                 fnal_pie_charts = fnal_pie_charts)
}


rf_rmse_sig_plot <- model_ggplots$RF$rmse_plot
gb_rmse_sig_plot <- model_ggplots$GB$rmse_plot
rf_rmse_sig_plot_notitle <- rf_rmse_sig_plot + theme(plot.title = element_blank()) + 
  theme(plot.margin = margin(t = -60, r = 5, b = 0, l = 5))
rf_rsquare_sig_plot_notitle <- rf_rsquare_sig_plot + theme(plot.title = element_blank()) + 
  theme(plot.margin = margin(t = -60, r = 5, b = 0, l = 5))
gb_rmse_sig_plot_notitle <- gb_rmse_sig_plot + theme(plot.title = element_blank()) + 
  theme(plot.margin = margin(t = -60, r = 5, b = 0, l = 5))
gb_rsquare_sig_plot_notitle <- gb_rsquare_sig_plot + theme(plot.title = element_blank()) + 
  theme(plot.margin = margin(t = -60, r = 5, b = 0, l = 5))

rf_title <- ggdraw() + draw_label("Random Forest", fontface = 'bold', size = 14, x = 0.5, hjust = 0.5)
gb_title <- ggdraw() + draw_label("Gradient Boosting", fontface = 'bold', size = 14, x = 0.5, hjust = 0.5)
rf_row <- plot_grid(rf_rmse_sig_plot_notitle, rf_rsquare_sig_plot_notitle, ncol = 2)
gb_row <- plot_grid(gb_rmse_sig_plot_notitle, gb_rsquare_sig_plot_notitle, ncol = 2)

residual_rsqu_plot <- plot_grid(
  rf_title, rf_row,
  gb_title, gb_row,
  ncol = 1,
  rel_heights = c(0.15, 1, 0.15, 1),
  labels = c("A", "", "B", ""),
  label_size = 20
)
ggsave(filename = paste0(directory,"/",response,"_mergedpredict_rmsersq.png"), plot = residual_rsqu_plot, device = "png", bg = 'white')

rf_imp_pie_plot <- model_ggplots$RF$fnal_pie_charts
ggsave(filename = paste0(directory,"/",response,"_rfweights.png"), plot = rf_imp_pie_plot, device = "png", bg = 'white')

gb_imp_pie_plot <- model_ggplots$GB$fnal_pie_charts
ggsave(filename = paste0(directory,"/",response,"_gbweights.png"), plot = gb_imp_pie_plot, device = "png", bg = 'white')
