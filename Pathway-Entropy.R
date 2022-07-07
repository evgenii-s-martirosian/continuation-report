#BiocManager::install("Biobase")
#BiocManager::install("GEOquery")
#BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
library(BiocManager)
library(Biobase)
library(GEOquery)
library(limma)
library(lumi)
library(gplots)
library(gtools)
library(mixOmics)
library(DMwR)
library(Boruta)
library(mlbench)
library(rattle)
library(BioQC) # package to calculate entropy
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(tidyverse)
#install.packages("ROCR")
#install.packages( "~/R/x86_64-pc-linux-gnu-library/4.1/DMwR_0.4.1.tar.gz", repos=NULL, type="source" )
#linux
# setwd("/media/mdehsfs4/PEC/Affiliated_people/Sara Cuvertino/DNAmethyl- blood/Butcher Cuvertino Sobreira")
#windows
setwd("~/DNA-methylation-paper/")

data("IlluminaHumanMethylation450kanno.ilmn12.hg19")
anno<-getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19) #HG19
anno<-as.data.frame(anno)
# setwd("C:/Users/mfbx9tg4/Desktop/Current Data/Kabuki Methylome Cuvertino Butcher Sobreira/")
comb<-read.table("Combined-Butcher-Cuvertino-and-Sobreira-Data.txt",sep="\t") # CpG data
status<-read.table("Combined-Butcher-Cuvertino-and-Sobreira-Gender-and-Disease-State.txt",sep="\t",header=T)

dat<-comb

probes<-read.table("Combined-Butcher-and-Cuvertino-Data_With-Normalisations_Chromosomes.txt",header=T,sep="\t") # CpG position and chromosome
probes<-probes[which(probes$CHR!="X"),] # Removes chromosome X
probes$CHR<-factor(probes$CHR,levels = c(seq(1,22,1))) # Creates a class factor for chromosomes 1-22

probe_anno<-anno[na.omit(match(probes$Variable.annotation.1,as.character(rownames(anno)))),]

dat<-dat[na.omit(match(probes$Variable.annotation.1,rownames(dat))),]
dat1 <- tibble::rownames_to_column(dat, "X")
status$Disease.State<-substr(status$Disease.State,2,nchar(as.character(status$Disease.State)))
go_terms <- read.csv("dna-met-paper-go-terms-combined-genes.txt", sep="\t", col.names = c("Gene_name", "Description", "Process"))
sample_size<-50
rn<-rownames(dat) # Extract CpG names from dat object
dat<-apply(dat,2,as.numeric) # Takes column from dat and converts it into type "numeric", creates matrix
rownames(dat)<-rn
Kabuki<-dat[,which(status$Disease.State=="Kabuki")] # Extracts Kabuki data
Control<-dat[,which(status$Disease.State=="Control")] # Extracts Control Data

# Creating Gene column with only one Gene Symbol
probe_anno <- probe_anno %>% 
  separate(UCSC_RefGene_Name, c("Gene", NA))

# Creating list of GO Term dataframes
go_terms <- split(go_terms, go_terms$Process)



set.seed(1)

res_list_control <- list()

# Control-------------------
for(i in 1:length(go_terms)){
  res <- data.frame(matrix(NA,nrow=2,ncol=2))
  colnames(res)<-c("Process", "Entropy")
  # Extracted CpGs within genes
  j <- probe_anno[probe_anno$Gene %in% go_terms[[i]][["Gene_name"]],]
  extracted_cpgs <- dat1 %>%
    filter(dat1$X %in% rownames(j))
  for(n in 1:100){
    rand<-sample(extracted_cpgs$X, size = sample_size,replace=F)
    cormat<-cor(t(Control[na.omit(match(rand,rownames(Control))),]),t(Control[-na.omit(match(rand,rownames(Control))),]))
    std<-sd(cormat)
    bin<-abs(cormat)
    bin[which(bin>=std)]<-1
    bin[which(bin<std)]<-0
    hyp<-bin%*%t(bin)
    
    hm<-heatmap.2(hyp,trace="none")
    dend<-as.hclust(hm$rowDendrogram)
    k<-2
    ct<- cutree(dend, k)
    Control_cc<-names(ct[which(ct==1)])
    res$Process[n]<- go_terms[[i]][["Process"]][1]
    res$Entropy[n]<-mean(apply(hyp[which(ct==1),which(ct==1)],2,entropy))
  }
  res_list_control[[length(res_list_control)+1]]<-res
}

control_pathway_entropy <- do.call(rbind, res_list_control)
write.csv(control_pathway_entropy, "control_pathway_entropy.csv")


res_list_kabuki <- list()

# Kabuki-------------------
for(i in 1:length(go_terms)){
  res <- data.frame(matrix(NA,nrow=2,ncol=2))
  colnames(res)<-c("Process", "Entropy")
  # Extracted CpGs within genes
  j <- probe_anno[probe_anno$Gene %in% go_terms[[i]][["Gene_name"]],]
  extracted_cpgs <- dat1 %>%
    filter(dat1$X %in% rownames(j))
  for(n in 1:100){
    rand<-sample(extracted_cpgs$X, size = sample_size,replace=F)
    cormat<-cor(t(Kabuki[na.omit(match(rand,rownames(Kabuki))),]),t(Kabuki[-na.omit(match(rand,rownames(Kabuki))),]))
    std<-sd(cormat)
    bin<-abs(cormat)
    bin[which(bin>=std)]<-1
    bin[which(bin<std)]<-0
    hyp<-bin%*%t(bin)
    
    hm<-heatmap.2(hyp,trace="none")
    dend<-as.hclust(hm$rowDendrogram)
    k<-2
    ct<- cutree(dend, k)
    Kabuki_cc<-names(ct[which(ct==1)])
    res$Process[n]<- go_terms[[i]][["Process"]][1]
    res$Entropy[n]<-mean(apply(hyp[which(ct==1),which(ct==1)],2,entropy))
  }
  res_list_kabuki[[length(res_list_kabuki)+1]]<-res
}

kabuki_pathway_entropy <- do.call(rbind, res_list_kabuki)
write.csv(kabuki_pathway_entropy, "kabuki_pathway_entropy.csv")
