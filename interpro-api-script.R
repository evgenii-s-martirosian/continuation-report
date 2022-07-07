library(tidyverse)
library("httr")
library("jsonlite")
library("tidyverse")
library("dplyr")

setwd("C:/Users/44784/OneDrive - The University of Manchester/PhD_Genomics/Chromatin-related genes/GO-terms-new-22-04-22/Chromatin_list_final")

chromatin_gene_list <- read_tsv("Final_Chromatin_Gene_List-with-UniProtKB-IDs.txt")
uniprotkb_ids_vec <- as_vector(na.omit(chromatin_gene_list$`UniProtKB ID`))


# API extract interpro entries
#https://www.ebi.ac.uk/interpro/api/protein/uniprot/Q9UI42/entry/integrated
base <- "https://www.ebi.ac.uk/interpro/api/protein/uniprot/"
endpoint <- "entry/integrated"

list_interpro_info_dfs <- list()
for(uniprot in uniprotkb_ids_vec){
  print(uniprot)
  call <- paste0(base, uniprot, "/", endpoint)
  get_interpro_info <- GET(call)
  get_interpro_info_text <- try(content(get_interpro_info, "text"))
  get_interpro_info_json <- try(fromJSON(get_interpro_info_text, flatten = TRUE))
  interpro_df <- try(select(get_interpro_info_json[["entry_subset"]], -entry_protein_locations))
  interpro_df$`UniProtKB ID` <- uniprot
  list_interpro_info_dfs[[length(list_interpro_info_dfs)+1]]<-interpro_df
}

list_interpro_info_dfs <- list_interpro_info_dfs[- 341]
interpro_info_df <- do.call(rbind, list_interpro_info_dfs)
Q86X51_df <- tibble(accession = NA,
                    protein_length = NA,
                    source_database = NA,
                    entry_type = NA,
                    entry_integrated = NA,
                    `UniProtKB ID` = "Q86X51")


interpro_info_df <- rbind(interpro_info_df, Q86X51_df)

write_tsv(interpro_info_df, "interpro_info_chromatin_genes.txt")

interpro_info_df_domain_only <- interpro_info_df%>%
  filter(entry_type == "domain")

write_tsv(interpro_info_df_domain_only, "interpro_info_chromatin_genes_domain_only.txt")

unique_interpro_ids <- unique(interpro_info_df_domain_only$entry_integrated)


#API extract interpro names
base <- "https://www.ebi.ac.uk/interpro/api/entry/interpro/"


list_of_domain_names <- list()
for (interpro_id in  unique_interpro_ids){
  print(interpro_id)
  call <- paste0(base, interpro_id, "/")
  get_interpro_info <- GET(call)
  get_interpro_info_text <- content(get_interpro_info, "text")
  get_interpro_info_json <- fromJSON(get_interpro_info_text, flatten = TRUE)
  interpro_name <- get_interpro_info_json[["metadata"]][["name"]]
  interpro_name$entry_integrated <- interpro_id
  list_of_domain_names[[length(list_of_domain_names)+1]]<-interpro_name
}

domain_names_df <- do.call(rbind, list_of_domain_names)
domain_names_df <- as_tibble(domain_names_df)
domain_names_df <- apply(domain_names_df,2,as.character)
domain_names_df <- as_tibble(domain_names_df)
write_tsv(domain_names_df, "chromatin-related-genes-unique-protein-domains.tsv")


# API extract backaground genes information
base <- "https://www.ebi.ac.uk/interpro/api/protein/reviewed/proteome/uniprot/"
proteome <- "UP000005640"

# Make a list for results
proteome_human_list <- list()

#First API request
call <- paste0(base, proteome, "/?page_size=200")
get_proteome_info <- GET(call)
get_proteome_info_text <- content(get_proteome_info, "text")
get_proteome_info_json <- fromJSON(get_proteome_info_text, flatten = TRUE)

call <- get_proteome_info_json[["next"]] # Extract call for subsequent API request

# Get results from the request 1

proteome_human_list[[length(proteome_human_list)+1]]<-get_proteome_info_json[["results"]] # Append results to list


while (call != "null"){
  get_proteome_info <- GET(call)
  get_proteome_info_text <- content(get_proteome_info, "text")
  get_proteome_info_json <- fromJSON(get_proteome_info_text, flatten = TRUE)
  
  call <- get_proteome_info_json[["next"]]
  proteome_human_list[[length(proteome_human_list)+1]]<-get_proteome_info_json[["results"]]
}
human_proteome <- do.call(rbind, proteome_human_list)

human_proteome <- human_proteome %>% 
  select(-proteomes)

write_tsv(human_proteome, "human_proteome.tsv")

# Get gene names for a human proteome
base <- "https://www.ebi.ac.uk/proteins/api/genecentric/"
proteome_uniprot_ids_vec <- as_vector(human_proteome$metadata.accession)

proteome_gene_name_list <- list()
for (id in proteome_uniprot_ids_vec){
  call <- paste0(base, id)
  get_gene_name <- GET(call)
  get_gene_name <- content(get_gene_name, "text")
  get_gene_name <- fromJSON(get_gene_name, flatten = TRUE)
  get_gene_name <- as_tibble(get_gene_name[["relatedGene"]])
  get_gene_name <- get_gene_name %>%
    try(filter(geneNameType == "Gene name"))
  proteome_gene_name_list[[length(proteome_gene_name_list)+1]]<-get_gene_name
}

proteome_gene_name_df <- do.call(rbind, proteome_gene_name_list)
proteome_gene_name_df <- proteome_gene_name_df %>%
  filter(geneNameType == "Gene name")

proteome_gene_name_df <- proteome_gene_name_df %>%
  rename(`UniProtKB ID` = accession)

# Download whole human proteome info
base <- "https://www.ebi.ac.uk/interpro/api/protein/uniprot/"
endpoint <- "entry/integrated"
uniprotkb_ids_vec <- as_vector(human_proteome$metadata.accession)


#call <- paste0(base, "A5LHX3", "/", endpoint)

list_interpro_info_dfs <- list()
for(uniprot in uniprotkb_ids_vec){
  print(uniprot)
  call <- paste0(base, uniprot, "/", endpoint)
  get_interpro_info <- GET(call)
  get_interpro_info_text <- try(content(get_interpro_info, "text"))
  get_interpro_info_json <- try(fromJSON(get_interpro_info_text, flatten = TRUE))
  interpro_df <- try(select(get_interpro_info_json[["entry_subset"]], -entry_protein_locations), silent = TRUE)
  if (length(interpro_df) > 2){
    interpro_df$`UniProtKB ID` <- uniprot
    list_interpro_info_dfs[[length(list_interpro_info_dfs)+1]]<-interpro_df
  } else { 
    print("subscript out of bounds")
  }
}



interpro_proteome_df <- do.call(rbind, list_interpro_info_dfs)
interpro_proteome_df <- interpro_proteome_df %>%
  filter(entry_type == "domain")


proteome_df <- left_join(interpro_proteome_df, proteome_gene_name_df, by = "UniProtKB ID")

# Pull protein domain names
interpro_ids_vec <- as_vector(proteome_df$entry_integrated)
domain_names_list <- list()
base <- "https://www.ebi.ac.uk/interpro/api/entry/interpro/"

for (id in interpro_ids_vec){
  call <- paste0(base, "ipr034248", "/")
  get_domain_name <- GET(call)
  get_domain_name <- content(get_domain_name, "text")
  get_domain_name <- fromJSON(get_domain_name, flatten = TRUE)
  get_domain_name <- get_domain_name[["metadata"]][["name"]]
  get_domain_name$entry_integrated <- id
  domain_names_list[[length(domain_names_list)+1]]<-get_domain_name
}
proteome_domain_names <- do.call(rbind, domain_names_list)
proteome_domain_names <- apply(proteome_domain_names,2,as.character)
proteome_domain_names <- as_tibble(proteome_domain_names)


proteome_df <- full_join(proteome_df, proteome_domain_names, by = "entry_integrated")
#call <- paste0(base, "ipr034248", "/")
#get_domain_name <- GET(call)
#get_domain_name <- content(get_domain_name, "text")
#get_domain_name <- fromJSON(get_domain_name, flatten = TRUE)
#get_domain_name[["metadata"]][["name"]]
#final_protein_domain_data <- full_join(domain_names_df, interpro_proteome_df, by = "entry_integrated")



#write_tsv(final_protein_domain_data, "final_protein_domain_data.tsv")
