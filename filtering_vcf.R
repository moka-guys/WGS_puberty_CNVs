library(tidyverse) # V1.3.1

# This script
# 1) Tidies DGV Gold SV track 
# 2) Filters VCFs
# Filtering steps
# Remove CNVs which haven't passed
# Remove CNVs not called by at least two callers from Parliment
# Remove common variants using DGV gold  
## example VCF ##
# Starts at 31,782
# Pass takes it to 8,376
# removing calls only made by one caller 3,1782


# ===================== #
## Tidy DGV data  ##  
# ===================== #

# DGV gold standard variants is a curated set of variants from a number of studies from DGV
# Taken from http://dgv.tcag.ca/dgv/docs/DGV.GS.hg38.gff3
# Release date 2016-05-15
setwd("/home/erin/Documents/Work/SNP_array_liftover/dgv_gold_track_38") 
dgv_gold_df <- read.delim("DGV.GS.hg38.gff3", sep = '\t' , header = FALSE) 

col_names_two_less_dgv <- c( "Name", "variant_type",  "outer_start", "inner_start", "inner_end", "outer_end", "inner_rank", "num_variants",
                             "variants", "num_studies", "Studies", "num_platforms", "Platforms", "number_of_algorithms", "algorithms","num_samples", "samples", 
                             "Frequency","PopulationSummary","Number_of_unique_samples_tested")

col_names_dgv <- c("ID", "Name", "variant_type", "variant_sub_type", "outer_start", "inner_start", "inner_end", "outer_end", "inner_rank", "num_variants",
                   "variants", "num_studies", "Studies", "num_platforms", "Platforms", "number_of_algorithms", "algorithms","num_samples", "samples", 
                   "Frequency","PopulationSummary","Number_of_unique_samples_tested")


tidy_dgv_gold_df <- dgv_gold_df %>% 
  separate(V9, col_names_dgv, sep = ";") %>% # Separate the final column into multiple columns
  pivot_longer(all_of(col_names_two_less_dgv), names_to = "name", values_to = "value") %>% # Flip df to easily manipulate data
  mutate(value_alone = str_remove(value, ".*=")) %>% # Remove string up to "="
  select(-c(value)) %>% # Drop unnecessary column 
  pivot_wider(names_from = name, values_from = value_alone) %>% 
  mutate(variant_sub_type = str_remove(variant_sub_type, ".*=")) %>% 
  rename(chrom = V1, # rename 
         start = inner_start, # thick regions are the ones with the highest confidence 
         end = inner_end) %>% 
  mutate(num_variants = as.integer(num_variants),
         num_studies = as.integer(num_studies),
         num_samples = as.integer(num_samples),
         Number_of_unique_samples_tested = as.integer(Number_of_unique_samples_tested),
         frequency_round = round(((as.integer(num_samples)/as.integer(Number_of_unique_samples_tested))*100), digits = 3)) %>% 
  select(chrom, start, end, num_variants, num_studies, num_samples, Frequency, frequency_round, Number_of_unique_samples_tested, ID) %>%
  filter(frequency_round > 5.0) %>% # Aligns with ACMG guidance (2015) that variants a > 5% in control population are very likely benign 
  # take dataframe from 113,556 to 31,347
  select(chrom, start, end) %>% # Need to remove ID column if going to merg this after
  distinct()  # removes duplicated rows, generated from having thick and thin start/ends for each variant.
  # from 31,347 to 10,449

rm(dgv_gold_df)


# Save 
setwd("/home/erin/Documents/Work/WGS_puberty_CNV/WGS_puberty_CNVs/controls/")
write.table(tidy_dgv_gold_df , "DGV_GS_hg38_tidy.txt", sep ='\t', col.names = FALSE, row.names = FALSE, quote = FALSE) # Save df

#merged_dgv <- read.delim("DGV_GS_hg38_tidy_sorted_100_bp_regions_merged.txt", sep = '\t' , header = FALSE) 

# ===================== #
## # Filter patient VCF   ## 
# ===================== #

setwd("/home/erin/Documents/Work/WGS_puberty_CNV/WGS_puberty_CNVs/vcfs/")
sample_list <- readLines("sample_list.txt")
#sample_ID <- "300033_180909_I302_CL100084594_L1_HUMykkRAAAB-525_1_rmdup"

for (sample_ID in sample_list) {
  
  setwd("/home/erin/Documents/Work/WGS_puberty_CNV/WGS_puberty_CNVs/vcfs/")
  # Load VCF 
  vcf  <- read.delim(paste0(sample_ID,".combined.genotyped.vcf"), sep = '\t' , header = TRUE, skip = 3401) 

  names(vcf)[10] <- "sample"  # rename this column

  col_names_vcf <- c("SUPP", "SUPP_VEC" , "AVGLEN", "SVTYPE", "SVMETHOD", "CHR2", "END", 
               "CIPO", "CIEND", "STRANDS", "CALLERS")


filtered_vcf <- vcf %>% 
  filter(FILTER == "PASS") %>% # Remove CNVs which didn't pass filter 
  separate(INFO, col_names_vcf, sep = ";") %>% # pull the info column into many columns
  pivot_longer(all_of(col_names_vcf), names_to = "name", values_to = "value") %>% # Flip df to easily manipulate data
  mutate(value_alone = str_remove(value, ".*=")) %>% # Remove string up to "="
  select(-c(value)) %>% # Drop unnecessary column 
  pivot_wider(names_from = name, values_from = value_alone) %>% 
  filter(str_detect(CALLERS, ",")) %>% # , denotes that at least two CNV callers found this CNV
  select("X.CHROM",  "POS", "END", "ID", "REF", "ALT",
          "sample", "SUPP", "SUPP_VEC", "AVGLEN", "SVTYPE", "SVMETHOD", 
         "CHR2", "CIPO", "CIEND", "STRANDS","CALLERS") # re order columns to reflect BED file format

setwd("/home/erin/Documents/Work/WGS_puberty_CNV/WGS_puberty_CNVs/filtered_vcfs/")
  
write.table(filtered_vcf , paste0(sample_ID, ".filtered_vcf.txt"), sep ='\t', col.names = FALSE, row.names = FALSE, quote = FALSE) # Save df 

}
#### Run Bedtools (v2.26.0) in the command line via the bedtools_WGV_Puberty_CNV.sh script ###





