#!/bin/bash
# bedtools v2.26.0
# This script
# 1) Annotates the transcript bed file with gene symbols 
# 2) Tidies DGV Gold SV track 

# ===================== #
## Merge nearby DGV GS regions ## 
# ===================== #

# Sort into chromosome order before merging
sort -k1,1 -k2,2n DGV_GS_hg38_tidy.txt > DGV_GS_hg38_tidy_sorted.txt
# Merge sections together which overlap, merge regions which have 50 bp between them
# Take DGV gold file from 10,446 to 7,745
# -i inoput -d Maximum distance between features allowed for features to be merged.
bedtools merge -i DGV_GS_hg38_tidy_sorted.txt -d 50 > DGV_GS_hg38_tidy_sorted_50_bp_regions_merged.txt

# ===================== #
## Merge nearby BED file regions ## 
# ===================== #

sort -k1,1 -k2,2n genes_of_interest.bed > genes_of_interest_sorted.bed
# Merge sections together which overlap, merge regions which have 100 bp between them
# Take bedfile from 524 to 88
bedtools merge -i genes_of_interest_sorted.bed > genes_of_interest_sorted_merged.bed 

# ===================== #
## Intersect gene symbols bed file with transcript ## 
# ===================== #
# intersect and take the genomic locations from the transcript bed file, 
# adding on the HGNCID and gene symbol
bedtools intersect -a genes_of_interest_sorted_merged.bed -b bed_file_with_gene.bed -wa -wb \
  | cut -f1,2,3,7,8,9 > transcript_annotated_genes_of_interest.bed
