#Require the package so you can use it
library("httr")
library("jsonlite")
library("tidyverse")
library("dplyr")

# Request structure
#base <- "http://api.geneontology.org/api/bioentity/function/"
#endpoint <- "genes"
#NCBITaxon <- "NCBITaxon:9606"
#GOTerm <- "GO:0044598"

# GET Request
#call <- paste(base, GOTerm, "/" ,endpoint, sep="")
#call

#get_goterm_info <- GET(call)
#get_goterm_info_text <- content(get_goterm_info, "text")
#get_goterm_info_text
#get_goterm_info_json <- fromJSON(get_goterm_info_text, flatten = TRUE)
#tibble_go_term <- as_tibble(get_goterm_info_json[["associations"]])
#tibble_go_term %>%
#  filter(subject.taxon.label == "Homo sapiens")

# Module 1. Analyse go terms synonymous to "histone" search word

go_terms_histone_mf_vector <-  c("GO:0072371", "GO:0072354", "GO:1990889", "GO:0042800", "GO:0042799",
                                 "GO:0042054", "GO:0017136", "GO:0120301", "GO:0120295", "GO:0120297",
                                 "GO:0160009", "GO:0160012", "GO:0004402", "GO:0004407", "GO:0018024", 
                                 "GO:0045129", "GO:0097372", "GO:0044020", "GO:0044022", "GO:0044025", 
                                 "GO:0044016", "GO:0044017", "GO:0044013", "GO:0046969", "GO:0046976", 
                                 "GO:0046974", "GO:0046975", "GO:0046972", "GO:0046970", "GO:0032931",
                                 "GO:0034647", "GO:0032454", "GO:0032452", "GO:0032453", "GO:0034739",
                                 "GO:0031151", "GO:0070612", "GO:0070611", "GO:0031078", "GO:0008469", 
                                 "GO:0043994", "GO:0043993", "GO:0043996", "GO:0043995", "GO:0043998",
                                 "GO:0043997", "GO:0043999", "GO:0043992", "GO:1990162", "GO:0000511",
                                 "GO:0140765", "GO:0140751", "GO:0140713", "GO:0106229", "GO:1990226",
                                 "GO:1990259", "GO:1990244", "GO:0106078", "GO:0032041", "GO:0071558",
                                 "GO:0046811", "GO:0032129", "GO:0035173", "GO:0035175", "GO:0035174",
                                 "GO:0035184", "GO:0035034", "GO:0035033", "GO:0010484", "GO:0010485",
                                 "GO:0062122", "GO:0036408", "GO:0140683", "GO:0140684", "GO:0140680",
                                 "GO:0140681", "GO:0140682", "GO:0140674", "GO:0140665", "GO:0140566",
                                 "GO:0051864", "GO:0035575", "GO:0035403", "GO:0035402", "GO:0035401",
                                 "GO:0035400", "GO:0140068", "GO:0140069", "GO:0033749", "GO:0033746",
                                 "GO:0035642", "GO:0097030", "GO:0140750", "GO:0140658")
                      
go_terms_histone_bf_vector <- c("GO:0070932", "GO:0070933", "GO:0072370", "GO:0072355", "GO:0044648",
                                "GO:1990596", "GO:1990678", "GO:1990679", "GO:1990619", "GO:0080182",
                                "GO:1990853", "GO:0016573", "GO:0016574", "GO:0016575", "GO:0016576",
                                "GO:0016577", "GO:0016578", "GO:0016570", "GO:0016571", "GO:0016572",
                                "GO:0043486", "GO:0098532", "GO:0070077", "GO:0070076", "GO:0070079",
                                "GO:0070078", "GO:0070544", "GO:0070535", "GO:0070537", "GO:0097692",
                                "GO:0097676", "GO:0097198", "GO:0097043", "GO:0097044", "GO:0044154",
                                "GO:0034770", "GO:0034773", "GO:0034772", "GO:0034771", "GO:0034721",
                                "GO:0034720", "GO:0034729", "GO:0070734", "GO:0033523", "GO:0033522",
                                "GO:0043972", "GO:0043971", "GO:0043974", "GO:0043973", "GO:0043977",
                                "GO:0043970", "GO:0043969", "GO:0043968", "GO:0043967", "GO:0043966",
                                "GO:0043990", "GO:0043991", "GO:0043983", "GO:0043982", "GO:0043985",
                                "GO:0043984", "GO:0043987", "GO:0043988", "GO:0043981", "GO:0043980",
                                "GO:0043979", "GO:1990164", "GO:0000412", "GO:0140372", "GO:0140373",
                                "GO:1990258", "GO:1990245", "GO:0106077", "GO:2000751", "GO:0071572",
                                "GO:0071557", "GO:0071110", "GO:0071894", "GO:0061647", "GO:0010452",
                                "GO:0010390", "GO:0036124", "GO:0036123", "GO:0036414", "GO:0036413",
                                "GO:0036351", "GO:0036352", "GO:0036353", "GO:0051567", "GO:0051568",
                                "GO:0035518", "GO:0035522", "GO:0035521", "GO:0035574", "GO:0035405",
                                "GO:0035404", "GO:0035409", "GO:0035408", "GO:0035407", "GO:0035406",
                                "GO:0034969", "GO:0034968", "GO:0034972", "GO:0034971", "GO:0034970",
                                "GO:0035978", "GO:0035616", "GO:0006334", "GO:0006337", "GO:0031508",
                                "GO:0031497", "GO:0034080")

tibble_go_terms_histone_mf_combined <- tibble('subject.label' = character(), 'provided_by' = list(), 'subject.taxon.label' = character(), 'object.label' = character())

for (go_term in go_terms_histone_mf_vector) {
  base <- "http://api.geneontology.org/api/bioentity/function/"
  endpoint <- "genes?rows=100000"
  call <- paste(base, go_term, "/" ,endpoint, sep="")
  get_goterm_info <- GET(call)
  get_goterm_info_text <- content(get_goterm_info, "text")
  get_goterm_info_json <- fromJSON(get_goterm_info_text, flatten = TRUE)
  tibble_go_term <- as_tibble(get_goterm_info_json[["associations"]])
  go_term_label <- gsub(" ",  "_", tibble_go_term$object.label[[1]])
  tibble_go_term <- tibble_go_term %>%
    filter(subject.taxon.label == "Homo sapiens") %>%
    select(subject.label, provided_by, subject.taxon.label, object.label)
  tibble_go_terms_histone_mf_combined <- rbind(tibble_go_terms_histone_mf_combined, tibble_go_term)
  #assign(paste0(go_term, "_gene_products"), tibble_go_term)
}

tibble_go_terms_histone_mf_combined <- tibble_go_terms_histone_mf_combined %>%
  distinct(subject.label, .keep_all = TRUE)

tibble_go_terms_histone_bf_combined <- tibble('subject.label' = character(), 'provided_by' = list(), 'subject.taxon.label' = character(), 'object.label' = character())

for (go_term in go_terms_histone_bf_vector) {
  base <- "http://api.geneontology.org/api/bioentity/function/"
  endpoint <- "genes?rows=100000"
  call <- paste(base, go_term, "/" ,endpoint, sep="")
  get_goterm_info <- GET(call)
  get_goterm_info_text <- content(get_goterm_info, "text")
  get_goterm_info_json <- fromJSON(get_goterm_info_text, flatten = TRUE)
  tibble_go_term <- as_tibble(get_goterm_info_json[["associations"]])
  go_term_label <- gsub(" ",  "_", tibble_go_term$object.label[[1]])
  tibble_go_term <- tibble_go_term %>%
    filter(subject.taxon.label == "Homo sapiens") %>%
    select(subject.label, provided_by, subject.taxon.label, object.label)
  tibble_go_terms_histone_bf_combined <- rbind(tibble_go_terms_histone_bf_combined, tibble_go_term)
  #assign(paste0(go_term, "_gene_products"), tibble_go_term)
}


tibble_go_terms_histone_bf_combined <- tibble_go_terms_histone_bf_combined %>%
  distinct(subject.label, .keep_all = TRUE)

tibble_go_terms_histone_mf_bf_combined <- rbind(tibble_go_terms_histone_bf_combined, tibble_go_terms_histone_mf_combined)

tibble_go_terms_histone_mf_bf_combined <- tibble_go_terms_histone_mf_bf_combined %>%
  distinct(subject.label, .keep_all = TRUE)

# Module 2. Analyse go terms synonymous to "nucleosome" search word
go_terms_nucleosome_mf_vector <- c("GO:0140750", "GO:0004402", "GO:0140713", "GO:0140674", "GO:0140658", "GO:0140463")

tibble_go_terms_nucleosome_mf_combined <- tibble('subject.label' = character(), 'provided_by' = list(), 'subject.taxon.label' = character(), 'object.label' = character())

for (go_term in go_terms_nucleosome_mf_vector) {
  base <- "http://api.geneontology.org/api/bioentity/function/"
  endpoint <- "genes?rows=100000"
  call <- paste(base, go_term, "/" ,endpoint, sep="")
  get_goterm_info <- GET(call)
  get_goterm_info_text <- content(get_goterm_info, "text")
  get_goterm_info_json <- fromJSON(get_goterm_info_text, flatten = TRUE)
  tibble_go_term <- as_tibble(get_goterm_info_json[["associations"]])
  go_term_label <- gsub(" ",  "_", tibble_go_term$object.label[[1]])
  tibble_go_term <- tibble_go_term %>%
    filter(subject.taxon.label == "Homo sapiens") %>%
    select(subject.label, provided_by, subject.taxon.label, object.label)
  tibble_go_terms_nucleosome_mf_combined <- rbind(tibble_go_terms_nucleosome_mf_combined, tibble_go_term)
  #assign(paste0(go_term, "_gene_products"), tibble_go_term)
}

tibble_go_terms_nucleosome_mf_combined <- tibble_go_terms_nucleosome_mf_combined %>%
  distinct(subject.label, .keep_all = TRUE)

go_terms_nucleosome_bf_vector <- c("GO:0006334", "GO:0006337", "GO:0006335", "GO:0006336", "GO:0034080")

tibble_go_terms_nucleosome_bf_combined <- tibble('subject.label' = character(), 'provided_by' = list(), 'subject.taxon.label' = character(), 'object.label' = character())

for (go_term in go_terms_nucleosome_bf_vector) {
  base <- "http://api.geneontology.org/api/bioentity/function/"
  endpoint <- "genes?rows=100000"
  call <- paste(base, go_term, "/" ,endpoint, sep="")
  get_goterm_info <- GET(call)
  get_goterm_info_text <- content(get_goterm_info, "text")
  get_goterm_info_json <- fromJSON(get_goterm_info_text, flatten = TRUE)
  tibble_go_term <- as_tibble(get_goterm_info_json[["associations"]])
  go_term_label <- gsub(" ",  "_", tibble_go_term$object.label[[1]])
  tibble_go_term <- tibble_go_term %>%
    filter(subject.taxon.label == "Homo sapiens") %>%
    select(subject.label, provided_by, subject.taxon.label, object.label)
  tibble_go_terms_nucleosome_bf_combined <- rbind(tibble_go_terms_nucleosome_bf_combined, tibble_go_term)
  #assign(paste0(go_term, "_gene_products"), tibble_go_term)
}


tibble_go_terms_nucleosome_bf_combined <- tibble_go_terms_nucleosome_bf_combined %>%
  distinct(subject.label, .keep_all = TRUE)

tibble_go_terms_nucleosome_mf_bf_combined <- rbind(tibble_go_terms_nucleosome_bf_combined, tibble_go_terms_nucleosome_mf_combined)

tibble_go_terms_nucleosome_mf_bf_combined <- tibble_go_terms_nucleosome_mf_bf_combined %>%
  distinct(subject.label, .keep_all = TRUE)

# Module 3. Analyse go terms synonymous to "chromatin" search word
go_terms_chromatin_mf_vector <- c("GO:0030527", "GO:0140658", "GO:0140587", "GO:0140463", "GO:0004407",
                                  "GO:0140750", "GO:0140751", "GO:0046811", "GO:0140585", "GO:0140586",
                                  "GO:0140566")

tibble_go_terms_chromatin_mf_combined <- tibble('subject.label' = character(), 'provided_by' = list(), 'subject.taxon.label' = character(), 'object.label' = character())

for (go_term in go_terms_chromatin_mf_vector) {
  base <- "http://api.geneontology.org/api/bioentity/function/"
  endpoint <- "genes?rows=100000"
  call <- paste(base, go_term, "/" ,endpoint, sep="")
  get_goterm_info <- GET(call)
  get_goterm_info_text <- content(get_goterm_info, "text")
  get_goterm_info_json <- fromJSON(get_goterm_info_text, flatten = TRUE)
  tibble_go_term <- as_tibble(get_goterm_info_json[["associations"]])
  go_term_label <- gsub(" ",  "_", tibble_go_term$object.label[[1]])
  tibble_go_term <- tibble_go_term %>%
    filter(subject.taxon.label == "Homo sapiens") %>%
    select(subject.label, provided_by, subject.taxon.label, object.label)
  tibble_go_terms_chromatin_mf_combined <- rbind(tibble_go_terms_chromatin_mf_combined, tibble_go_term)
  #assign(paste0(go_term, "_gene_products"), tibble_go_term)
}

tibble_go_terms_chromatin_mf_combined <- tibble_go_terms_chromatin_mf_combined %>%
  distinct(subject.label, .keep_all = TRUE)

go_terms_chromatin_bf_vector <- c("GO:0006335", "GO:0006338", "GO:0006336", "GO:1990700", "GO:0034723",
                                  "GO:0034724", "GO:0031497", "GO:0031055", "GO:0034080", "GO:0061641",
                                  "GO:0140673", "GO:0140588", "GO:0006334", "GO:0006337", "GO:0006346",
                                  "GO:0070828", "GO:0030466", "GO:0030261", "GO:0080188", "GO:1990893",
                                  "GO:0043486", "GO:0031507", "GO:0031508", "GO:0031509", "GO:0140698",
                                  "GO:0140464", "GO:0140462")

tibble_go_terms_chromatin_bf_combined <- tibble('subject.label' = character(), 'provided_by' = list(), 'subject.taxon.label' = character(), 'object.label' = character())

for (go_term in go_terms_chromatin_bf_vector) {
  base <- "http://api.geneontology.org/api/bioentity/function/"
  endpoint <- "genes?rows=100000"
  call <- paste(base, go_term, "/" ,endpoint, sep="")
  get_goterm_info <- GET(call)
  get_goterm_info_text <- content(get_goterm_info, "text")
  get_goterm_info_json <- fromJSON(get_goterm_info_text, flatten = TRUE)
  tibble_go_term <- as_tibble(get_goterm_info_json[["associations"]])
  go_term_label <- gsub(" ",  "_", tibble_go_term$object.label[[1]])
  tibble_go_term <- tibble_go_term %>%
    filter(subject.taxon.label == "Homo sapiens") %>%
    select(subject.label, provided_by, subject.taxon.label, object.label)
  tibble_go_terms_chromatin_bf_combined <- rbind(tibble_go_terms_chromatin_bf_combined, tibble_go_term)
  #assign(paste0(go_term, "_gene_products"), tibble_go_term)
}


tibble_go_terms_chromatin_bf_combined <- tibble_go_terms_chromatin_bf_combined %>%
  distinct(subject.label, .keep_all = TRUE)

tibble_go_terms_chromatin_mf_bf_combined <- rbind(tibble_go_terms_chromatin_bf_combined, tibble_go_terms_chromatin_mf_combined)

tibble_go_terms_chromatin_mf_bf_combined <- tibble_go_terms_chromatin_mf_bf_combined %>%
  distinct(subject.label, .keep_all = TRUE)

# Module 4. Combine all of the datasets
tibble_combined <- rbind(tibble_go_terms_histone_mf_bf_combined,
                         tibble_go_terms_chromatin_mf_bf_combined,
                         tibble_go_terms_nucleosome_mf_bf_combined)

tibble_combined <- tibble_combined %>%
  distinct(subject.label, .keep_all = TRUE) %>%
  filter(!grepl('_human', subject.label))

write_tsv(tibble_combined, file = "C:/Users/44784/OneDrive - The University of Manchester/PhD_Genomics/Chromatin-related genes/gene-list-previous-go-terms.tsv")

# Module 5. Using new GO Terms and fetching the data
new_go_terms <- c("GO:0030527", "GO:0035034", "GO:0035033", "GO:0140713", "GO:0140587",
                  "GO:0140463", "GO:0140707", "GO:0000118", "GO:0035514", "GO:0009008",
                  "GO:0046811", "GO:0000123", "GO:0035097", "GO:0035034", "GO:0035033",
                  "GO:0016514", "GO:0070822", "GO:0000786", "GO:0070603", "GO:0070209",
                  "GO:0061638", "GO:1990483", "GO:0035101", "GO:0070823", "GO:0033503",
                  "GO:0070210", "GO:0034967", "GO:0070211", "GO:1904173")


tibble_new_go_terms_combined <- tibble('subject.label' = character(), 'provided_by' = list(), 'subject.taxon.label' = character(), 'object.label' = character())

for (go_term in new_go_terms) {
  base <- "http://api.geneontology.org/api/bioentity/function/"
  endpoint <- "genes?rows=100000"
  call <- paste(base, go_term, "/" ,endpoint, sep="")
  get_goterm_info <- GET(call)
  get_goterm_info_text <- content(get_goterm_info, "text")
  get_goterm_info_json <- fromJSON(get_goterm_info_text, flatten = TRUE)
  tibble_go_term <- as_tibble(get_goterm_info_json[["associations"]])
  go_term_label <- gsub(" ",  "_", tibble_go_term$object.label[[1]])
  tibble_go_term <- tibble_go_term %>%
    filter(subject.taxon.label == "Homo sapiens") %>%
    select(subject.label, provided_by, subject.taxon.label, object.label)
  tibble_new_go_terms_combined <- rbind(tibble_new_go_terms_combined, tibble_go_term)
  #assign(paste0(go_term, "_gene_products"), tibble_go_term)
}

tibble_new_go_terms_combined <- tibble_new_go_terms_combined %>%
  distinct(subject.label, .keep_all = TRUE) %>%
  filter(!grepl('_human', subject.label))

write_tsv(tibble_new_go_terms_combined, file = "C:/Users/44784/OneDrive - The University of Manchester/PhD_Genomics/Chromatin-related genes/GO-terms-new-08-04-22/gene-list-new-go-terms.tsv")

# Module 6. Fetch basic info for the new GO 

tibble_new_go_terms_info <- tibble('goid' = character(), 'label' = character(), 'definition' = character())

for (go_term in new_go_terms) {
  base <- "http://api.geneontology.org/api/ontology/term/"
  call <- paste(base, go_term, sep="")
  get_goterm_info <- GET(call)
  get_goterm_info_text <- content(get_goterm_info, "text")
  get_goterm_info_json <- fromJSON(get_goterm_info_text, flatten = TRUE)
  goid <- get_goterm_info_json[["goid"]]
  label <- get_goterm_info_json[["label"]]
  definition <- get_goterm_info_json[["definition"]]
  get_goterm_info_tibble <- tibble(goid, label, definition)
  tibble_new_go_terms_info <- rbind(tibble_new_go_terms_info, get_goterm_info_tibble)
}
write_tsv(tibble_new_go_terms_info, file = "C:/Users/44784/OneDrive - The University of Manchester/PhD_Genomics/Chromatin-related genes/GO-terms-new-08-04-22/new-go-terms-info.tsv")


# Module 7. Using bigger GO-Terms based on the GO hierarchy visualisation of previosuly used go terms from 08-04-22
## Fetch basic info for the new Go-Terms

`22_04_22_go_terms` <- c("GO:0016570", "GO:0004407", "GO:0030527", "GO:0006325", "GO:0070603", "GO:0006338",
                         "GO:0006306", "GO:0044815", "GO:0031497", "GO:0035101", "GO:0016575", "GO:0004402",
                         "GO:0004407", "GO:0033503", "GO:0140713", "GO:0140707", "GO:0140587", "GO:0140463",
                         "GO:0000118", "GO:0004407", "GO:0070603", "GO:0061638", "GO:0031056", "GO:0046811",
                         "GO:0031063", "GO:0004407", "GO:0035033", "GO:0070209", "GO:0035514")
tibble_new_go_terms_info <- tibble('goid' = character(), 'label' = character(), 'definition' = character())

for (go_term in `22_04_22_go_terms`) {
  base <- "http://api.geneontology.org/api/ontology/term/"
  call <- paste(base, go_term, sep="")
  get_goterm_info <- GET(call)
  get_goterm_info_text <- content(get_goterm_info, "text")
  get_goterm_info_json <- fromJSON(get_goterm_info_text, flatten = TRUE)
  goid <- get_goterm_info_json[["goid"]]
  label <- get_goterm_info_json[["label"]]
  definition <- get_goterm_info_json[["definition"]]
  get_goterm_info_tibble <- tibble(goid, label, definition)
  tibble_new_go_terms_info <- rbind(tibble_new_go_terms_info, get_goterm_info_tibble)
}

tibble_new_go_terms_info <- tibble_new_go_terms_info %>%
  distinct(goid, .keep_all = TRUE)

write_tsv(tibble_new_go_terms_info, file = "C:/Users/44784/OneDrive - The University of Manchester/PhD_Genomics/Chromatin-related genes/GO-terms-new-22-04-22/new-go-terms-info-22-04-22.tsv")

# Using new GO Terms and fetching the data
tibble_new_go_terms_combined <- tibble('subject.label' = character(), 'provided_by' = list(), 'subject.taxon.label' = character(), 'object.label' = character())

for (go_term in `22_04_22_go_terms`) {
  base <- "http://api.geneontology.org/api/bioentity/function/"
  endpoint <- "genes?rows=100000"
  call <- paste(base, go_term, "/" ,endpoint, sep="")
  get_goterm_info <- GET(call)
  get_goterm_info_text <- content(get_goterm_info, "text")
  get_goterm_info_json <- fromJSON(get_goterm_info_text, flatten = TRUE)
  tibble_go_term <- as_tibble(get_goterm_info_json[["associations"]])
  go_term_label <- gsub(" ",  "_", tibble_go_term$object.label[[1]])
  tibble_go_term <- tibble_go_term %>%
    filter(subject.taxon.label == "Homo sapiens") %>%
    select(subject.label, provided_by, subject.taxon.label, object.label)
  tibble_new_go_terms_combined <- rbind(tibble_new_go_terms_combined, tibble_go_term)
  #assign(paste0(go_term, "_gene_products"), tibble_go_term)
}

tibble_new_go_terms_combined <- tibble_new_go_terms_combined %>%
  distinct(subject.label, .keep_all = TRUE) %>%
  filter(!grepl('_human', subject.label)) %>%
  filter(!provided_by == "RNAcentral")

unique_object_labels <- unique(tibble_new_go_terms_combined$object.label)

unique_object_labels_as_tibble <- as_tibble(unique_object_labels) # Go through all of this and remove genes based on object.labels
write_tsv(unique_object_labels_as_tibble, file = "C:/Users/44784/OneDrive - The University of Manchester/PhD_Genomics/Chromatin-related genes/GO-terms-new-22-04-22/unique-go-terms-22-04-22.tsv")

terms_to_exclude <- c("cohesin complex", "condensin complex", "CST complex", "DNA topoisomerase III-beta-TDRD3 complex", "dosage compensation", "dosage compensation by inactivation of X chromosome",
                      "epigenetic maintenance of chromatin in transcription-competent conformation", "genomic imprinting", "heterochromatin assembly by small RNA", "meiotic cohesin complex",
                      "mitotic cohesin complex", "Ndc80 complex", "negative regulation of gene expression, epigenetic", "NMS complex", "nuclear meiotic cohesin complex",
                      "nuclear mitotic cohesin complex", "nuclear telomere cap complex", "oxidative DNA demethylase activity", "programmed DNA elimination by chromosome breakage",
                      "random inactivation of X chromosome", "rDNA heterochromatin assembly", "regulation of gene expression by genomic imprinting", "regulation of gene expression, epigenetic",
                      "shelterin complex", "sperm DNA condensation", "spermatogenesis, exchange of chromosomal proteins")

excluded_terms_tibble_new <- tibble_new_go_terms_combined %>%
  filter(!object.label == "cohesin complex" & !object.label == "condensin complex" & !object.label == "CST complex" & !object.label == "DNA topoisomerase III-beta-TDRD3 complex" & !object.label == "dosage compensation" &
           !object.label == "dosage compensation by inactivation of X chromosome" & !object.label == "epigenetic maintenance of chromatin in transcription-competent conformation" & !object.label == "genomic imprinting" &
           !object.label == "heterochromatin assembly by small RNA" & !object.label == "meiotic cohesin complex" & !object.label == "mitotic cohesin complex" & !object.label == "Ndc80 complex" & !object.label == "negative regulation of gene expression, epigenetic" &
           !object.label == "NMS complex" & !object.label == "nuclear meiotic cohesin complex" & !object.label == "nuclear mitotic cohesin complex" & !object.label == "nuclear telomere cap complex" & !object.label == "oxidative DNA demethylase activity" &
           !object.label == "programmed DNA elimination by chromosome breakage" & !object.label == "random inactivation of X chromosome" & !object.label == "rDNA heterochromatin assembly" & !object.label == "regulation of gene expression by genomic imprinting" &
           !object.label == "regulation of gene expression, epigenetic" & !object.label == "shelterin complex" & !object.label == "sperm DNA condensation" & !object.label == "spermatogenesis, exchange of chromosomal proteins")

write_tsv(excluded_terms_tibble_new, file = "C:/Users/44784/OneDrive - The University of Manchester/PhD_Genomics/Chromatin-related genes/GO-terms-new-22-04-22/gene-list-22-04-22-filtered-terms.tsv")


# Workflow 22-04-22 (Also get super GO-Term chromatin and repeat this approach)
  # 1. Get GO-Terms from the printed out visualisations of 08-04-22 GO-Term hierarchies +
  # 2. Use API to fetch data on GO-TErms and associated human genes +
  # 3. Download list of Unique go-term labels from the fetched genes based on GO-terms from 22/04/22 +
  # 4. Manually go through the GO-term list from point 3 and choose what to keep
  # 5. Re-run point 2 with go terms from point 4

# Module 8. Use very large go terms and do analysis as in module 7.
`23_04_22_go_terms` <- c("GO:0005654", "GO:0000785", "GO:0006325", "GO:0006304")

tibble_new_go_terms_info <- tibble('goid' = character(), 'label' = character(), 'definition' = character())

for (go_term in `23_04_22_go_terms`) {
  base <- "http://api.geneontology.org/api/ontology/term/"
  call <- paste(base, go_term, sep="")
  get_goterm_info <- GET(call)
  get_goterm_info_text <- content(get_goterm_info, "text")
  get_goterm_info_json <- fromJSON(get_goterm_info_text, flatten = TRUE)
  goid <- get_goterm_info_json[["goid"]]
  label <- get_goterm_info_json[["label"]]
  definition <- get_goterm_info_json[["definition"]]
  get_goterm_info_tibble <- tibble(goid, label, definition)
  tibble_new_go_terms_info <- rbind(tibble_new_go_terms_info, get_goterm_info_tibble)
}

tibble_new_go_terms_info <- tibble_new_go_terms_info %>%
  distinct(goid, .keep_all = TRUE)

write_tsv(tibble_new_go_terms_info, file = "C:/Users/44784/OneDrive - The University of Manchester/PhD_Genomics/Chromatin-related genes/GO-terms-new-23-04-22/new-go-terms-info-23-04-22.tsv")

# Using new GO Terms and fetching the data
tibble_new_go_terms_combined <- tibble('subject.label' = character(), 'provided_by' = list(), 'subject.taxon.label' = character(), 'object.label' = character())

for (go_term in `23_04_22_go_terms`) {
  base <- "http://api.geneontology.org/api/bioentity/function/"
  endpoint <- "genes?rows=1000000000"
  call <- paste(base, go_term, "/" ,endpoint, sep="")
  get_goterm_info <- GET(call)
  get_goterm_info_text <- content(get_goterm_info, "text")
  get_goterm_info_json <- fromJSON(get_goterm_info_text, flatten = TRUE)
  tibble_go_term <- as_tibble(get_goterm_info_json[["associations"]])
  go_term_label <- gsub(" ",  "_", tibble_go_term$object.label[[1]])
  tibble_go_term <- tibble_go_term %>%
    filter(subject.taxon.label == "Homo sapiens") %>%
    select(subject.label, provided_by, subject.taxon.label, object.label)
  tibble_new_go_terms_combined <- rbind(tibble_new_go_terms_combined, tibble_go_term)
  #assign(paste0(go_term, "_gene_products"), tibble_go_term)
}

tibble_new_go_terms_combined <- tibble_new_go_terms_combined %>%
  distinct(subject.label, .keep_all = TRUE) %>%
  filter(!grepl('_human', subject.label)) %>%
  filter(!provided_by == "RNAcentral")

unique_object_labels <- unique(tibble_new_go_terms_combined$object.label)

unique_object_labels_as_tibble <- as_tibble(unique_object_labels) # Go through all of this and remove genes based on object.labels
write_tsv(unique_object_labels_as_tibble, file = "C:/Users/44784/OneDrive - The University of Manchester/PhD_Genomics/Chromatin-related genes/GO-terms-new-23-04-22/unique-go-terms-23-04-22.tsv")

