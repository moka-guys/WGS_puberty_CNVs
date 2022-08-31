library(tidyverse) # V1.3.1
library(biomaRt) # V2.46.3
## This script 
# Needs to be run in sections, alongside gene_discovery.sh
# This script will say when to go and run other bits 
# These two scripts facilitate the annotation of genetic regions not found in DGV gold 
# and not found in the genes of interest given for this project 
# to facilitate gene discovery 

##==============================================
## STEP ONE ##

setwd("/home/erin/Documents/Work/WGS_puberty_CNV/WGS_puberty_CNVs/vcfs/")
sample_list <- readLines("sample_list.txt")

for (sample_ID in sample_list) {
  print(paste0("Working on sample", sample_ID))
  setwd("/home/erin/Documents/Work/WGS_puberty_CNV/WGS_puberty_CNVs/gene_discovery")
  # Load VCF 

  vcf  <- read.delim(paste0(sample_ID,".no_bed_or_dgv_overlap.txt"), sep = '\t' , header = FALSE) %>%  
    dplyr::select(V1, V2, V3, V4, V9, V13, V15, V16, V17, V18, V19) # drop column, not needed 

  names(vcf) <- c("CHROM",  "POS", "END", "ID", "GENOTYPE", "SVTYPE", "CHR2", 
                  "CIPO", "CIEND", "STRANDS", "CALLERS")
  
  # Filter samples to the by SV type and genotype/ 
  filtered_vcf_del_heterozygous <- vcf %>% 
    filter(SVTYPE == "DEL" & GENOTYPE == "0/1" ) %>% 
    dplyr::select("CHROM",  "POS", "END", "ID") %>% 
    distinct()
  
  filtered_vcf_del_homozygous <- vcf %>% 
    filter(SVTYPE == "DEL" & GENOTYPE == "1/1" ) %>% 
    dplyr::select("CHROM",  "POS", "END", "ID") %>% 
    distinct()
  
  filtered_vcf_inv_heterozygous <- vcf %>% 
    filter(SVTYPE == "INV" & GENOTYPE == "0/1" ) %>% 
    dplyr::select("CHROM",  "POS", "END", "ID") %>% 
    distinct()
  
  filtered_vcf_inv_homozygous <- vcf %>% 
    filter(SVTYPE == "INV" & GENOTYPE == "1/1" ) %>% 
    dplyr::select("CHROM",  "POS", "END", "ID") %>% 
    distinct()
  
  # Save!
  setwd("/home/erin/Documents/Work/WGS_puberty_CNV/WGS_puberty_CNVs/gene_discovery/to_merge")
  write.table(filtered_vcf_del_heterozygous , paste0(sample_ID, ".DELS.HETEROZYGOUS.txt"), sep ='\t', col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  write.table(filtered_vcf_del_homozygous , paste0(sample_ID, ".DELS.HOMOZYGOUS.txt"), sep ='\t', col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  write.table(filtered_vcf_inv_heterozygous , paste0(sample_ID, ".INVS.HETEROZYGOUS.txt"), sep ='\t', col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  write.table(filtered_vcf_inv_homozygous , paste0(sample_ID, ".INVS.HOMOZYGOUS.txt"), sep ='\t', col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  
}

##==============================================

### RUN step two of gene_discovery.sh #####

##==============================================
## STEP TWO ##

ensembl_genes <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = "uswest") # hsapiens_gene_ensembl = Human genes (GRCh38.p13) GRCh38.p13


setwd("/home/erin/Documents/Work/WGS_puberty_CNV/WGS_puberty_CNVs/vcfs/")
sample_list <- readLines("sample_list.txt")


for (sample_ID in sample_list) {
  print(paste0("Working on sample ", sample_ID))
  setwd("/home/erin/Documents/Work/WGS_puberty_CNV/WGS_puberty_CNVs/gene_discovery/to_merge/")
  
  # Load VCF
  patient_merged_dels <- read.delim(paste0(sample_ID,".dels.merged.txt"), sep = '\t' , header = FALSE) 
  
  patient_merged_invs  <- read.delim(paste0(sample_ID,".invs.merged.txt"), sep = '\t' , header = FALSE) 
   
  # Run once for inversions 
  ensembl_coordinates_tidy_invs <- patient_merged_invs %>% 
    mutate(V1 = str_remove(V1, "chr")) %>% 
    unite(genomic_location, c(V1,V2, V3), sep = ":") %>% # Join them into one for searching
    distinct() %>% 
    dplyr::select(genomic_location)
  
  # Annotate genes for each region
  genes_invs = getBM(attributes=c('chromosome_name','start_position','end_position','hgnc_id', 'hgnc_symbol'), 
                filters = c("chromosomal_region"),
                values = ensembl_coordinates_tidy_invs$genomic_location ,
                mart=ensembl_genes)
  
  
  # Run a second time for deletions 
  ensembl_coordinates_tidy_dels <- patient_merged_dels %>% 
    mutate(V1 = str_remove(V1, "chr")) %>% 
    unite(genomic_location, c(V1,V2, V3), sep = ":") %>% # Join them into one for searching
    distinct() %>% 
    dplyr::select(genomic_location)
  
  
  genes_dels = getBM(attributes=c('chromosome_name','start_position','end_position','hgnc_id', 'hgnc_symbol'), 
                     filters = c("chromosomal_region"),
                     values = ensembl_coordinates_tidy_invs$genomic_location ,
                     mart=ensembl_genes)
  
genes_dels_to_merge <- genes_dels %>% 
  mutate(sv_type = "DEL")  %>% 
  filter(str_detect(hgnc_id, "HGNC")) %>% # only keep results which are in coding regions 
  mutate(chromosome_name = str_c("chr", chromosome_name))
  
genes_invs_to_merge <- genes_invs %>% 
  mutate(sv_type = "INV")  %>% 
  filter(str_detect(hgnc_id, "HGNC")) %>% 
  mutate(chromosome_name = str_c("chr", chromosome_name))

  # Save to then intersect with unique CNV IDs
  setwd("/home/erin/Documents/Work/WGS_puberty_CNV/WGS_puberty_CNVs/gene_discovery/annotated")
  write.table(genes_invs_to_merge , paste0(sample_ID, ".INVS.gene.bed"), sep ='\t', col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  write.table(genes_dels_to_merge , paste0(sample_ID, ".DELS.gene.bed"), sep ='\t', col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  rm(vcf)
  rm(ensembl_coordinates_tidy)
  rm(genes_tidy)
 rm(genes)
  
} 

##==============================================

## Run step three gene_discovery.sh ###

##==============================================
## STEP THREE ##

setwd("/home/erin/Documents/Work/WGS_puberty_CNV/WGS_puberty_CNVs/vcfs/")
sample_list <- readLines("sample_list.txt")


for (sample_ID in sample_list) {
  
  setwd("/home/erin/Documents/Work/WGS_puberty_CNV/WGS_puberty_CNVs/gene_discovery/annotated")
  # Load gene annotated regions with no overlap
  annotated_invs  <- read.delim(paste0(sample_ID,".gene_bed_no_overlap_invs.txt"), sep = '\t' , header = FALSE) 
  annotated_dels  <- read.delim(paste0(sample_ID,".gene_bed_no_overlap_dels.txt"), sep = '\t' , header = FALSE) 
  
  annotated <- bind_rows(annotated_invs, annotated_dels )
  
  gene_contact_annotated <- annotated %>% 
    distinct() %>% 
    group_by(V4) %>%  # this is the unique ID for each CNV called 
    mutate(genes_concat = paste0(V9, collapse = " , ")) %>% 
    dplyr::select(V1, V2, V3, V4, V10, genes_concat) %>%  
    distinct()  
  
  write.table(gene_contact_annotated , paste0(sample_ID, ".gene_symbols_concat.txt"), sep ='\t', col.names = FALSE, row.names = FALSE, quote = FALSE)

}

# This step finds the genes (or collections of genes), which are found in at least 20 of the samples 

all_sample_genes = list()


for (sample_ID in sample_list) {
  
  setwd("/home/erin/Documents/Work/WGS_puberty_CNV/WGS_puberty_CNVs/gene_discovery/annotated")
  # Load gene annotated regions with no overlap
  annotated  <- read.delim(paste0(sample_ID,".gene_symbols_concat.txt"), sep = '\t' , header = FALSE) 
  
  annoted_tidy <- annotated %>% 
    dplyr::select(V6, V5) %>% 
    distinct() %>% 
    mutate(patient = sample_ID)

  all_sample_genes[[sample_ID]] <- annoted_tidy #
} 

df_bound_genes <- dplyr::bind_rows(all_sample_genes) %>% 
  group_by(V5, V6 ) %>% # include V5 here to split count by Dels or Inversions
  distinct() %>% 
  mutate(ascount = n()) %>% 
  group_by(V5, count, V6) %>% 
  summarise(
    note = paste0(patient, collapse = ",")
  )

names(df_bound_genes) <- c("GENE/S", "COUNT", "SV_TYPE", "IDs")

# Split by number to save on tabs in another script

first_df <- df_bound_genes %>% 
  filter(between(COUNT, 30,40)) %>% 
  dplyr::select("GENE/S", "COUNT", "SV_TYPE", "IDs")

second_df <- df_bound_genes %>% 
  filter(between(COUNT, 20, 39)) %>% 
  dplyr::select("GENE/S", "COUNT", "SV_TYPE", "IDs")

third_df <- df_bound_genes %>% 
  filter(between(COUNT, 10, 19)) %>% 
  dplyr::select("GENE/S", "COUNT", "SV_TYPE", "IDs")

fourth_df <- df_bound_genes %>% 
  filter(between(COUNT, 4, 9)) %>% 
  dplyr::select("GENE/S", "COUNT", "SV_TYPE", "IDs")


setwd("/home/erin/Documents/Work/WGS_puberty_CNV/WGS_puberty_CNVs/gene_discovery/results")

write.table(first_df , "gene_discovery_dgv_gold_30_40.txt", sep ='\t', col.names = TRUE, row.names = FALSE, quote = FALSE)

write.table(second_df , "gene_discovery_dgv_gold_20_29.txt", sep ='\t', col.names = TRUE, row.names = FALSE, quote = FALSE)

write.table(third_df , "gene_discovery_dgv_gold_10_19.txt", sep ='\t', col.names = TRUE, row.names = FALSE, quote = FALSE)

write.table(fourth_df , "gene_discovery_dgv_gold_4_9.txt", sep ='\t', col.names = TRUE, row.names = FALSE, quote = FALSE)

