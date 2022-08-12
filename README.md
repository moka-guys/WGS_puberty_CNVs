# Scripts for filtering WGS Puberty CNV data

## What do these scripts do? 

Filter WGS CNV Puberty delay VCFs in the following steps 

## Run filterign_vcf.R

#### Tidy DGV Gold standard variants 
  1) DGV gold standard variants (taken from http://dgv.tcag.ca/dgv/docs/DGV.GS.hg38.gff3, release date 2016-05-15)
  2) Takes the thick regions for the coordinates for each variant (which has the highest confidence)
  3) Removes variants not seen at > 5.0% of the population to align with ACMG guidance for frequency of likely benign variants in the population 
#### Filter VCFs
 1) Remove CNVs which didn't meet Parliament2 filtering steps 
 2) Remove CNVs which were only called by a single caller. CNVs called by ≥ callers are kept


## Run bedtools_WGS_Puberty_CNV.sh
#### Merge nearby DGV Gold regions 
 1) Merge regions of the DGV gold file produced by filtering_vcf.R so that regions within 100 bp are merged 
 
#### Merge nearby bed file regions 
  1) Transcripts for genes of interest in genes_of_interest.bed merged into regions which overlap 


#### Intersect VCF with various other files 
  1) Find results from filtered VCF which overlap with regions in the DGV Gold file. 50% overlap of DGV gold & variant from VCF must be met 
  2) Find results from filtered VCF which do not overlap with regions in the DGV Gold file. < 50% overlap of DGV gold & variant from VCF must be met
  3) Filter regions based on genes of interest bed file (genes_of_interest_sorted_merged.bed)
  
  
## Run merge_outputs.py

1) Pull in all the different data into dataframes 
2) Set them to output on a separate tab, of a single .xlxs sheet per patient 
3) Save!



## These scripts were developed by Viapath Genome Informatic