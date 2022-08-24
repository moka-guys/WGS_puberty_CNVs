library(tidyverse) # V1.3.1

# This script
# 1) Pulls in all the various sample data  
# 2) Counts the rows in each file & saves as a csv

# ===================== #
   ## Counts ##  
# ===================== #

# Count the rows in each file of the filtering process 
# Save as a csv
setwd("/home/erin/Documents/Work/WGS_puberty_CNV/WGS_puberty_CNVs/vcfs/")
sample_list <- readLines("sample_list.txt")

df_list = list()

for (sample_ID in sample_list) {
  
  # Load data per patient
  setwd("/home/erin/Documents/Work/WGS_puberty_CNV/WGS_puberty_CNVs/vcfs/")
  vcf  <- read.delim(paste0(sample_ID,".unfiltered_vcf.txt"), sep = '\t' , header = FALSE)

  setwd("/home/erin/Documents/Work/WGS_puberty_CNV/WGS_puberty_CNVs/filtered_vcfs/")
  filtered_vcf  <- read.delim(paste0(sample_ID, ".filtered_vcf.txt"), sep = '\t' , header = FALSE) 
  dgv_overlap_regions <- read.delim(paste0(sample_ID,".dgv_gold_overlap.txt"), sep = '\t' , header = FALSE) 
  no_dgv_overlap_regions <- read.delim(paste0(sample_ID,".no_dgv_gold_overlap.txt"), sep = '\t' , header = FALSE)
  bed_filtered <- read.delim(paste0(sample_ID,".bed_filtered_regions_unsorted_dgv_gold.txt"), sep = '\t' , header = FALSE) 

  # Count data 
  
  vcf_count <- vcf %>%
    summarise(n = n()) %>% 
    mutate(file = "VCF",
           sample = sample_ID)
  
  filtered_vcf_count <-  filtered_vcf %>% 
    summarise(n = n()) %>% 
    mutate(file = "VCF_filtered", 
           sample = sample_ID)
  
  dgv_overlap_regions_counts <- dgv_overlap_regions %>% 
    summarise(n = n()) %>% 
    mutate(file = "dgv_overlap_regions",
           sample = sample_ID)  
 
  no_dgv_overlap_regions_counts <- no_dgv_overlap_regions %>% 
    summarise(n = n()) %>% 
    mutate(file = "no_dgv_overlap_regions",
           sample = sample_ID)
  
  bed_filtered_count <- bed_filtered  %>% 
    summarise(n = n()) %>% 
    mutate(file = "bed_file_filtered_regions",
           sample = sample_ID) 
  
  # Merge per sample
  merged_df <- rbind(vcf_count, filtered_vcf_count, dgv_overlap_regions_counts, no_dgv_overlap_regions_counts, bed_filtered_count)
  
  df_list[[sample_ID]] <- merged_df # Put into a list
  
}

# Combine the list together, pivot and sum
df_bound_counts <- dplyr::bind_rows(df_list) %>% 
  pivot_wider(id_cols = c(file), names_from = sample, values_from = n) %>% 
  mutate(sum = rowSums(across(where(is.numeric))))

setwd("/home/erin/Documents/Work/WGS_puberty_CNV/WGS_puberty_CNVs/results")
write.table(df_bound_counts , "result_counts.txt", sep ='\t', col.names = TRUE, row.names = FALSE, quote = FALSE) # Save df 