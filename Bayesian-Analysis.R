setwd("~/DNA-methylation-paper")
library(tidyverse)
library(BayesFactor)
library(ggplot2)
library(bayestestR)

Entropy_data_control<-read.csv("control_pathway_entropy.csv",header=T,row.names=1,stringsAsFactors = T)
Entropy_data_kabuki<-read.csv("kabuki_pathway_entropy.csv",header=T,row.names=1,stringsAsFactors = T)
Entropy_data_control_galois<-read.csv("control_pathway_entropy_galois_cpg.csv",header=T,row.names=1,stringsAsFactors = T)
Entropy_data_kabuki_galois<-read.csv("kabuki_pathway_entropy_galois_cpg.csv",header=T,row.names=1,stringsAsFactors = T)

Entropy_data_control <- Entropy_data_control %>%
  mutate(Group = "Control")

Entropy_data_kabuki <- Entropy_data_kabuki %>%
  mutate(Group = "Kabuki")

Entropy_Data_Combined <- rbind(Entropy_data_control, Entropy_data_kabuki)

Entropy_data_control_galois <- Entropy_data_control_galois %>%
  mutate(Group = "Control")

Entropy_data_kabuki_galois <-  Entropy_data_kabuki_galois %>%
  mutate(Group = "Kabuki")

Entropy_Data_Combined_galois_cpgs <- rbind(Entropy_data_control_galois, Entropy_data_kabuki_galois)


#Human----
Entropy_subset_significant <- Entropy_Data_Combined %>%
  filter(Process == "Microtubule bundle formation" |
           Process == "Interleukin-17-mediated signaling pathway" |
           Process == "Sensory perception of light stimulus")


Entropy_Data_Combined
ggplot(Entropy_subset_significant,aes(x=Process,y=Entropy,col=Group))+geom_boxplot()+
  coord_flip()+
  theme_bw(base_size=18)+
  ylab(label="Entropy")+
  xlab(label="Gene Ontology Process")
  

ls<-list()
for(i in 1:nlevels(Entropy_Data_Combined$Process)){
  ls[[i]]<-Entropy_Data_Combined[Entropy_Data_Combined$Process %in% levels(Entropy_Data_Combined$Process)[i],]
  ls[[i]]$Process<-factor(ls[[i]]$Process)
}

res_combined<-list()
for(i in 1:nlevels(Entropy_Data_Combined$Process)){
  res_combined[[i]]<-data.frame("mu"=rep(NA,10000),"beta"=rep(NA,10000),
                       "sigma"=rep(NA,10000),
                       process=rep(levels(Entropy_Data_Combined$Process)[i],10000))
  bf<-ttestBF(formula=Entropy~Group,data=ls[[i]])
  chains = posterior(bf, iterations = 10000)
  res_combined[[i]]$mu<-as.numeric(chains[,1])
  res_combined[[i]]$beta<-as.numeric(chains[,2])
  res_combined[[i]]$sigma<-as.numeric(chains[,3])
}
#res_combined <- do.call(rbind, res_combined)

#res_combined %>%
#  ci(ci = 0.89, method = "HDI")

ls<-list()
for(i in 1:nlevels(Entropy_Data_Combined_galois_cpgs$Process)){
  ls[[i]]<-Entropy_Data_Combined_galois_cpgs[Entropy_Data_Combined_galois_cpgs$Process %in% levels(Entropy_Data_Combined_galois_cpgs$Process)[i],]
  ls[[i]]$Process<-factor(ls[[i]]$Process)
}

res_combined_galois_cpgs<-list()
for(i in 1:nlevels(Entropy_Data_Combined_galois_cpgs$Process)){
  res_combined_galois_cpgs[[i]]<-data.frame("mu"=rep(NA,10000),"beta"=rep(NA,10000),
                                "sigma"=rep(NA,10000),
                                process=rep(levels(Entropy_Data_Combined_galois_cpgs$Process)[i],10000))
  bf<-ttestBF(formula=Entropy~Group,data=ls[[i]])
  chains = posterior(bf, iterations = 10000)
  res_combined_galois_cpgs[[i]]$mu<-as.numeric(chains[,1])
  res_combined_galois_cpgs[[i]]$beta<-as.numeric(chains[,2])
  res_combined_galois_cpgs[[i]]$sigma<-as.numeric(chains[,3])
}

#res_combined_galois_cpgs <-  do.call(rbind, res_combined_galois_cpgs)

# 89th percent credibility intervals on beta values
ci_res_list_combined <-list()
for(i in 1:length(res_combined)){
  ci_calc <- data.frame(matrix(NA,nrow=1,ncol=4))
  ci_hdi <- ci(res_combined[[i]][["beta"]], ci = 0.89, method = "HDI")
  ci_calc$process <- res_combined[[i]][["process"]][1]
  ci_calc$CI <- ci_hdi$CI
  ci_calc$CI_low <- ci_hdi$CI_low
  ci_calc$CI_high <- ci_hdi$CI_high
  ci_res_list_combined[[length(ci_res_list_combined)+1]]<-ci_calc
}
ci_res_combined<-do.call(rbind,ci_res_list_combined)
ci_res_combined <- subset(ci_res_combined, select = -c(X1, X2, X3, X4))
res_combined <- do.call(rbind, res_combined)

ci_res_list_galois_cpgs <-list()
for(i in 1:length(res_combined_galois_cpgs)){
  ci_calc <- data.frame(matrix(NA,nrow=1,ncol=4))
  ci_hdi <- ci(res_combined_galois_cpgs[[i]][["beta"]], ci = 0.89, method = "HDI")
  ci_calc$process <- res_combined_galois_cpgs[[i]][["process"]][1]
  ci_calc$CI <- ci_hdi$CI
  ci_calc$CI_low <- ci_hdi$CI_low
  ci_calc$CI_high <- ci_hdi$CI_high
  ci_res_list_galois_cpgs[[length(ci_res_list_galois_cpgs)+1]]<-ci_calc
}

ci_res_galois_cpgs<-do.call(rbind,ci_res_list_galois_cpgs)
ci_res_galois_cpgs <- subset(ci_res_galois_cpgs, select = -c(X1, X2, X3, X4))
res_combined_galois_cpgs <- do.call(rbind, res_combined_galois_cpgs)

#Entropy_res<-do.call(rbind,res)

# 89th percent credibility intervals on beta values
#ci_res_list_combined <-list()
#for(i in 1:length(res)){
#  ci_calc <- data.frame(matrix(NA,nrow=1,ncol=4))
#  ci_hdi <- ci(res[[i]][["beta"]], ci = 0.89, method = "HDI")
#  ci_calc$process <- res[[i]][["process"]][1]
#  ci_calc$CI <- ci_hdi$CI
#  ci_calc$CI_low <- ci_hdi$CI_low
#  ci_calc$CI_high <- ci_hdi$CI_high
#  ci_res_list_combined[[length(ci_res_list_combined)+1]]<-ci_calc
#}



#ci_res <- do.call(rbind, ci_res_list)
#ci_res <- subset(ci_res, select = -c(X1, X2, X3, X4))

#Entropy_res_CI <- merge(Entropy_res, ci_res, by="process")
#Entropy_res %>%
#  group_by(process) %>%
#  summarise()
#Entropy_res %>%
  #group_by(process)%>%
  #summarise(
  #  ci_low = ci_res$CI_low[ci_res$process == process],
  #  ci_high = ci_res$CI_high[ci_res$process == process]) #%>%
  #ggplot()+ #aes(x=beta,y=reorder(process,beta,mean),group=process))+
  #geom_boxplot(aes(x=process, lower=ci_low, upper=ci_high))+ #x=beta)+
  #theme_bw(base_size=18)+
  #ylab(label="Gene Ontology Process")+
  #xlab(label="beta Entropy")




#head(Entropy_res)
#beta_ci_res <- merge(Entropy_res, ci_res, by="process")
#ci_res <- column_to_rownames(ci_res, "process")
#res[[1]][['beta']]

#Entropy_res_CI <- Entropy_res %>%
#  group_by(process) %>%
#  summarise(ci_low = ci_res$CI_low[ci_res$process == process],
#            ci_high = ci_res$CI_high[ci_res$process == process])


#Entropy_res_CI %>%
#  group_by(process) %>%
#  summarise(min = ci_res$CI_low[ci_res$process == process],
#            max = ci_res$CI_high[ci_res$process == process])


#res_combined %>%
#  group_by(process) %>%
#  summarise(mean_val = mean(beta),
#            ci_low = 

beta_combined <- ggplot(res_combined, aes(x=beta, y=process,group=process))+
  geom_violin()+
  geom_boxplot(width=0.15, xlower = ci_res_combined$CI_low, xupper = ci_res_combined$CI_high)+
  #stat_boxplot(geom = "boxplot", notchlower = ci_res$CI_low, notchupper = ci_res$CI_)+
  theme_bw(base_size=18)+
  ylab(label="Gene Ontology Process")+
  xlab(label="beta Entropy")# +
  #aes(y = reorder(process,mu,mean))

beta_combined

#reorder(res_combined$process, res_combined$beta, mean(res_combined$beta))



res_combined_galois_cpgs_reordered <- with(res_combined_galois_cpgs, reorder(process, beta, mean))

vector_reordered_galois_cpgs <- unlist(levels(res_combined_galois_cpgs_reordered))

ci_res_galois_cpgs <-ci_res_galois_cpgs %>%
  mutate(process =  factor(process, levels = vector_reordered_galois_cpgs)) %>%
  arrange(process)   



  


beta_galois_cpgs <- ggplot(res_combined_galois_cpgs, aes(x=beta,y=reorder(process,beta,mean),group=process))+
  geom_violin()+
  geom_boxplot(width=0.15, xlower = ci_res_galois_cpgs$CI_low, xupper = ci_res_galois_cpgs$CI_high)+
  #stat_boxplot(geom = "boxplot", notchlower = ci_res$CI_low, notchupper = ci_res$CI_)+
  theme_bw(base_size=18)+
  ylab(label="Gene Ontology Process")+
  xlab(label="beta Entropy")


beta_galois_cpgs

res_combined_reordered <- with(res_combined, reorder(process, beta, mean))

vector_reordered_combined <- unlist(levels(res_combined_reordered))

ci_res_combined <-ci_res_combined %>%
  mutate(process =  factor(process, levels = vector_reordered_combined)) %>%
  arrange(process)

beta_combined <- ggplot(res_combined, aes(x=beta,y=reorder(process,beta,mean),group=process))+
  geom_violin()+
  geom_boxplot(width=0.15, xlower = ci_res_combined$CI_low, xupper = ci_res_combined$CI_high)+
  #stat_boxplot(geom = "boxplot", notchlower = ci_res$CI_low, notchupper = ci_res$CI_)+
  theme_bw(base_size=18)+
  ylab(label="Gene Ontology Process")+
  xlab(label="beta Entropy")

beta_combined
res_combined_min_max <- res_combined %>%
  group_by(process) %>%
  summarise(min = min(beta),
            max = max(beta),
            median = median(beta))

res_combined_galois_cpgs_min_max <- res_combined_galois_cpgs %>%
  group_by(process) %>%
  summarise(min = min(beta),
            max = max(beta),
            median = median(beta))

CTRL_vs_KABUKI_TOP20 <- merge(res_combined_min_max, ci_res_combined, by="process")

CTRL_vs_KABUKI_GALOIS <- merge(res_combined_galois_cpgs_min_max, ci_res_galois_cpgs, by="process")

write.csv(CTRL_vs_KABUKI_TOP20, "CTRL_vs_KS_TOP20_GO_TERMS.csv", row.names = FALSE)
write.csv(CTRL_vs_KABUKI_GALOIS, "CTRL_vs_KS_GALOIS_GO_TERMS.csv", row.names = FALSE)
