library(tidyverse) # V1.3.1
library(biomaRt) # V2.46.3
library(readxl)  # V1.3.1

# This script
# 1) Annotates the transcript bed file with gene symbols 
# 2) Tidies DGV Gold SV track 

# ===================== #
## Annotate transcript bed file ##  
# ===================== #

## Load data ##

# adding gene symbols back on to the regions of interest 
setwd("/home/erin/Documents/Work/WGS_puberty_CNV/WGS_puberty_CNVs/reference") 
transcript_bed_file <- read.delim("genes_of_interest.bed", sep = '\t' , header = FALSE) 
hgnc_id_symbols <- read_excel("Candidates_DelayedPuberty_CNVanalysis.xlsx") 

hgnc_id_symbols$HGNC <- paste("HGNC:", hgnc_id_symbols$HGNC, sep="")

## Create mart for searching ##
ensembl_genes <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = "uswest") # hsapiens_gene_ensembl = Human genes (GRCh38.p13) GRCh38.p13

transcript_bed_file <- transcript_bed_file %>% 
  mutate(V1 = str_remove(V1, "chr")) %>% 
  unite(genomic_location, c(V1,V2, V3), sep = ":") %>% # Join them into one for searching
  distinct() %>% 
  mutate(strand = case_when(V6 == "+" ~ 1,
                             V6 == "-" ~ -1)) %>% # Change strand for searching in biomaRt
  dplyr::select(genomic_location, strand) 

# For each genomic location and strand return the attributes 
annotate_transcript_bed_file = getBM(attributes=c('chromosome_name','start_position','end_position','hgnc_id', 'hgnc_symbol', 'strand'), 
                   filters = c("chromosomal_region", "strand"),
                   values = list(transcript_bed_file$genomic_location, transcript_bed_file$strand)  ,
                   mart=ensembl_genes)

# The search above will have returned the chromosomal regions on both strands 
# So this joining with the original candidate HGNC IDs, removes unwanted data 
filtered_annotated_bed_file <- annotate_transcript_bed_file %>% 
  filter(str_detect(hgnc_id, "" )) %>% 
  mutate(strand = case_when(strand == "1" ~ "+",
                            strand == "-1" ~ "-")) %>% 
  semi_join(hgnc_id_symbols, by = c("hgnc_id" = "HGNC")) %>% 
  mutate(chromosome_name = str_c("chr", chromosome_name ))

# Save!

write.table(filtered_annotated_bed_file , "bed_file_with_gene.bed", sep ='\t',
            col.names = FALSE, row.names = FALSE, quote = FALSE) # Save df

# ===================== #
## Tidy DGV Gold Standard data  ##  
# ===================== #

# DGV gold standard variants is a curated set of variants from a number of studies from DGV
# Taken from http://dgv.tcag.ca/dgv/docs/DGV.GS.hg38.gff3
# Release date 2016-05-15
setwd("/home/erin/Documents/Work/WGS_puberty_CNV/WGS_puberty_CNVs/reference")
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
  dplyr::select(-c(value)) %>% # Drop unnecessary column 
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
  dplyr::select(chrom, start, end, num_variants, num_studies, num_samples, Frequency, frequency_round, Number_of_unique_samples_tested, ID) %>%
  filter(frequency_round > 5.0) %>% # Aligns with ACMG guidance (2015) that variants a > 5% in the population are very likely benign 
  # take dataframe from 113,556 to 31,347
  dplyr::select(chrom, start, end) %>% # Need to remove ID column if going to merg this after
  distinct()  # removes duplicated rows, generated from having thick and thin start/ends for each variant.
# from 31,347 to 10,449

# Save 
setwd("/home/erin/Documents/Work/WGS_puberty_CNV/WGS_puberty_CNVs/reference/")
write.table(tidy_dgv_gold_df , "DGV_GS_hg38_tidy.txt", sep ='\t', col.names = FALSE, row.names = FALSE, quote = FALSE) # Save df

# ===================== #
## Tidy DGV  data  ##  
# ===================== #
setwd("/home/erin/Documents/Work/WGS_puberty_CNV/WGS_puberty_CNVs/reference")

dgv_sv_standard <- read.delim("GRCh38_hg38_variants_2020-02-25.txt", sep = '\t' , header = TRUE) 

more_than_two_supporting_samples <- dgv_sv_standard %>% 
  filter(str_detect(chr, "" ),
         str_detect(supportingvariants, ",")) %>% # 827,036 to 202,513
  mutate(chr = str_c("chr", chr)) %>% 
  select(chr, start, end)

write.table(more_than_two_supporting_samples , "DGV_not_GS_hg38_tidy.txt", sep ='\t', col.names = FALSE, row.names = FALSE, quote = FALSE) # Save df

