library(tidyverse) # V1.3.1

# This script
# 1) Filters VCFs produced by Parliament  
# Filtering steps
# Remove SV which haven't passed
# Remove SV not called by at least two callers from Parliament2
## example VCF ##
# Starts at 31,782
# Pass takes it to 8,376
# removing calls only made by one caller 3,1782

# ===================== #
## Filter VCFs ##  
# ===================== #


setwd("/home/erin/Documents/Work/WGS_puberty_CNV/WGS_puberty_CNVs/vcfs/")
sample_list <- readLines("sample_list.txt")

col_names_vcf <- c("X.CHROM",  "POS", "ID", "REF", "ALT", "QUAL", "FILTER")
col_names_vcf_pivot <- c("SUPP", "SUPP_VEC" , "AVGLEN", "SVTYPE", "SVMETHOD", "CHR2", "END", 
                         "CIPO", "CIEND", "STRANDS", "CALLERS")


for (sample_ID in sample_list) {
  
  setwd("/home/erin/Documents/Work/WGS_puberty_CNV/WGS_puberty_CNVs/vcfs/")
  # Load VCF 
  vcf  <- read.delim(paste0(sample_ID,".combined.genotyped.vcf"), sep = '\t' , header = FALSE, skip = 3402) 

  names(vcf)[1:7] <- col_names_vcf

# Create an unfiltered vcf to save as an easier to understand output
tidy_vcf <- vcf %>% 
  separate(V8, col_names_vcf_pivot, sep = ";") %>% # pull the info column into many columns
  pivot_longer(all_of(col_names_vcf_pivot), names_to = "name", values_to = "value") %>% # Flip df to easily manipulate data
  mutate(value_alone = str_remove(value, ".*=")) %>% # Remove string up to "="
  dplyr::select(-c(value)) %>% # Drop unnecessary column 
  pivot_wider(names_from = name, values_from = value_alone) %>% 
  separate(V10, c("GT", "drop"), sep = ":") %>%
  dplyr::select("X.CHROM",  "POS", "END", "ID", "REF", "ALT", "QUAL", "FILTER",
         "GT", "SUPP", "SUPP_VEC", "AVGLEN", "SVTYPE", "SVMETHOD", 
         "CHR2", "CIPO", "CIEND", "STRANDS","CALLERS")
  
filtered_vcf <- tidy_vcf %>% 
  filter(FILTER == "PASS") %>% # Remove CNVs which didn't pass filter 
  filter(SUPP >= 2) # .+ 2 callers agree on the SV at this location 


setwd("/home/erin/Documents/Work/WGS_puberty_CNV/WGS_puberty_CNVs/vcfs/")
write.table(tidy_vcf , paste0(sample_ID, ".unfiltered_vcf.txt"), sep ='\t', col.names = FALSE, row.names = FALSE, quote = FALSE) # Save df

setwd("/home/erin/Documents/Work/WGS_puberty_CNV/WGS_puberty_CNVs/filtered_vcfs/")
write.table(filtered_vcf , paste0(sample_ID, ".filtered_vcf.txt"), sep ='\t', col.names = FALSE, row.names = FALSE, quote = FALSE) # Save df 

}
#### Run Bedtools (v2.26.0) in the command line via the bedtools_WGV_Puberty_CNV.sh script ###





