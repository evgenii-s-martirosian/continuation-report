library(tidyverse)
library(plyr)

chromatin_related_genes <- 
  read_tsv("C:/Users/44784/OneDrive - The University of Manchester/PhD_Genomics/Chromatin-related genes/GO-terms-new-22-04-22/chromatin_related_genes_ensembl_omim.tsv")

gnomad_lof_metrics <-
  read_tsv("C:/Users/44784/OneDrive - The University of Manchester/PhD_Genomics/Chromatin-related genes/GO-terms-new-22-04-22/gnomad.v2.1.1.lof_metrics.by_gene.txt")



gene_vector <- chromatin_related_genes$`Gene Name`

gnomad_lof_metrics <- gnomad_lof_metrics %>%
  select(gene, pLI, mis_z)

chromatin_genes_lof_metrics <- gnomad_lof_metrics %>%
  filter(gene %in% gene_vector) #%>%
  #unique()

chromatin_related_genes_and_pli <- 
  full_join(chromatin_related_genes, chromatin_genes_lof_metrics, by=c("Gene Name" = "gene"))
#
#chromatin_related_genes_and_pli <- unique(chromatin_related_genes_and_pli)


write_tsv(chromatin_related_genes_and_pli, "C:/Users/44784/OneDrive - The University of Manchester/PhD_Genomics/Chromatin-related genes/GO-terms-new-22-04-22/chromatin_related_genes_ensembl_omim_pli.tsv")

