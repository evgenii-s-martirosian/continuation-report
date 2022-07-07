# BiocManager:install("Biobase","GEOquery","limma","lumi","gplots","gtools","mixOmics","DMwR","Boruta","mlbench","rattle","BioQC","IlluminaHumanMethylation450kanno.ilmn12.hg19")
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
library(BioQC)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

setwd("/scratch/xxxxxem2/xxxxxem2/DNA-methylation-paper")

data("IlluminaHumanMethylation450kanno.ilmn12.hg19")
anno<-getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19) #HG19
anno<-as.data.frame(anno)
comb<-read.table("Combined-Butcher-Cuvertino-and-Sobreira-Data.txt",sep="\t")
status<-read.table("Combined-Butcher-Cuvertino-and-Sobreira-Gender-and-Disease-State.txt",sep="\t",header=T)

dat<-comb

probes<-read.table("./Combined-Butcher-and-Cuvertino-Data_With-Normalisations_Chromosomes.txt",header=T,sep="\t")
probes<-probes[which(probes$CHR!="X"),]
probes$CHR<-factor(probes$CHR,levels = c(seq(1,22,1)))

probe_anno<-anno[na.omit(match(probes$Variable.annotation.1,as.character(rownames(anno)))),]

dat<-dat[na.omit(match(probes$Variable.annotation.1,rownames(dat))),]
status$Disease.State<-substr(status$Disease.State,2,nchar(as.character(status$Disease.State)))

dmps<-read.csv("./Butcher-Cuvertino-Sobreira-DMP2.csv")
n<-length(which(dmps$adj.P.Val<0.0001))
rn<-rownames(dat)
dat<-apply(dat,2,as.numeric)
rownames(dat)<-rn
dat<-dat[na.omit(match(dmps$X,rownames(dat))),]
dmps<-dmps[na.omit(match(rownames(dat),dmps$X)),]

Kabuki<-dat[,which(status$Disease.State=="Kabuki")]
Control<-dat[,which(status$Disease.State=="Control")]

#Kabuki
cd<-cor(t(Kabuki[which(dmps$adj.P.Val<0.0001),]),t(Kabuki[-which(dmps$adj.P.Val<0.0001),]))
bin<-abs(cd)
bin[which(cd>sd(cd))]<-1
bin[which(bin!=1)]<-0
hyp<-bin%*%t(bin)

png("Kabuki_Hypernetwork_Heatmap.png")
hm<-heatmap.2(hyp,trace="none")
dev.off()
dend<-as.hclust(hm$rowDendrogram)
k<-2
ct<- cutree(dend, k)
Kabuki_cc<-names(ct[which(ct==1)])
Kabuki_gc<-cd[na.omit(match(Kabuki_cc,rownames(cd))),]
Kabuki_gc<-Kabuki_gc[,which(colSums(Kabuki_gc)>nrow(Kabuki_gc)*0.4)]
stats<-as.data.frame(matrix(NA,nrow=length(Kabuki_cc),ncol=1))
colnames(stats)<-c("entropy")
stats$entropy<-apply(hyp[na.omit(match(Kabuki_cc,rownames(hyp))),na.omit(match(Kabuki_cc,rownames(hyp)))],1,entropy)
write.csv(Kabuki_gc,"Kabuki_Galois_0,4.csv")
write.csv(stats,"Kabuki_Central_Cluster_Stats.csv")

#Control
cd<-cor(t(Control[which(dmps$adj.P.Val<0.0001),]),t(Control[-which(dmps$adj.P.Val<0.0001),]))
bin<-abs(cd)
bin[which(cd>sd(cd))]<-1
bin[which(bin!=1)]<-0
hyp<-bin%*%t(bin)

png("Control_Hypernetwork_Heatmap.png")
hm<-heatmap.2(hyp,trace="none")
dev.off()
dend<-as.hclust(hm$rowDendrogram)
k<-2
ct<- cutree(dend, k)
Control_cc<-names(ct[which(ct==1)])
Control_gc<-cd[na.omit(match(Control_cc,rownames(cd))),]
Control_gc<-Control_gc[,which(colSums(Control_gc)>nrow(Control_gc)*0.4)]
stats<-as.data.frame(matrix(NA,nrow=length(Control_cc),ncol=1))
colnames(stats)<-c("entropy")
stats$entropy<-apply(hyp[na.omit(match(Control_cc,rownames(hyp))),na.omit(match(Control_cc,rownames(hyp)))],1,entropy)
write.csv(Control_gc,"Control_Galois_0,4.csv")
write.csv(stats,"Control_Central_Cluster_Stats.csv")
