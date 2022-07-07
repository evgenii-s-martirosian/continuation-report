library(tidyverse)
library(ggplot2)
library(rstatix)
library(ggpubr)

setwd("C:/Users/44784/OneDrive - The University of Manchester/PhD_Genomics/iPSCs-stuff/DNA_methylation_paper/Entropy-DMPs")

Kabuki_Central_Cluster_Stats <- read_csv("Kabuki_Central_Cluster_Stats.csv")

Control_Central_Cluster_Stats <- read_csv("Control_Central_Cluster_Stats.csv")

Kabuki_Central_Cluster_Stats <- Kabuki_Central_Cluster_Stats %>%
  mutate(dataset = "Kabuki Experimental")

Control_Central_Cluster_Stats <- Control_Central_Cluster_Stats %>%
  mutate(dataset = "Control Experimental")

Central_Cluster_Stats_Combined <- rbind(Control_Central_Cluster_Stats, Kabuki_Central_Cluster_Stats)
# Basic violin plot
p <- ggplot(Central_Cluster_Stats_Combined, aes(x=dataset, y=entropy)) + 
  geom_violin()
p
p <- ggplot(Central_Cluster_Stats_Combined, aes(x=dataset, y=entropy)) + 
  geom_violin(trim=FALSE)

#Central_Cluster_Stats_Combined %>%
#  wilcox.test(entropy ~ dataset, alternative = "two.sided")
stat.test <- Central_Cluster_Stats_Combined %>% 
  wilcox_test(entropy ~ dataset) %>%
  add_significance()

stat.test <- stat.test %>% 
  mutate(y.position = 10.00)
# Plot
Central_Cluster_Stats_Combined %>%
  ggplot( aes(x=dataset, y=entropy)) +
  geom_violin(width=0.8, fill = "yellow") +
  geom_boxplot(width=0.1, color="black", alpha=0.2)+
  ylab("Entropy") +
  xlab("") +
  ylim(9.0, 10.2) +
  stat_pvalue_manual(stat.test) +
  theme_bw()
