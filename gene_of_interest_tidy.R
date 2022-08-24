library(tidyverse) # V1.3.1

# This script
# 1) Takes the sheet were all the patients are intersected the genes of interest 
# 2) Splits them into data frames per gene of interest

# ===================== #
      ## Counts ##  
# ===================== #

setwd("/home/erin/Documents/Work/WGS_puberty_CNV/WGS_puberty_CNVs/results") 
intersected_genes_of_interest <- read.delim("all_samples_intersected_with_genes_of_interest.txt", sep = '\t' , header = FALSE)

## Create different dfs for displaying results 
# Count of all the CNVs per gene, of each SV type
count_cnv_per_gene <- intersected_genes_of_interest %>% 
  dplyr::select(V1, V2, V3, V4,V5, V7, V20) %>% # genomic location (gene), HGNC ID, HGNC symbol, sample ID, SV type
  distinct()  %>% #
  group_by(V5, V20) %>% #HGNC Symbol & SV type 
  mutate(count = n()) %>% 
  ungroup() %>% 
  group_by(V5, V20, count) %>% 
  summarise(note = paste0(V7, collapse = ","))

write.table(count_cnv_per_gene ,"all_samples_intersected_with_genes_of_interest_count.txt",sep ='\t', col.names = FALSE, row.names = FALSE, quote = FALSE)

# Make a list of all the unique genes to loop through
unique_genes <- unique(intersected_genes_of_interest$V5)
gene <- "ACTL6B"
# Create and save a txt file per unique gene
for (gene in unique_genes) {
  
  unqiue <- intersected_genes_of_interest %>% 
    dplyr::select(V8, V9, V10, V5, V20, V16, # Genomic location (SV), HGNC Symbol, SV type, no of callers supporting, genotype
                  V17, V22, V23, V24, V25, V26, V7) %>% # SV end chr, PE CI around end and start, starnd, SV callers, Sample ID
    distinct() %>%  
    filter(V5 == gene) 
  
  write.table(unqiue , (paste0("CNV_",gene,".txt")),sep ='\t', col.names = FALSE, row.names = FALSE, quote = FALSE)
}