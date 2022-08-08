#!/usr/bin/bash -d

# This script
# Shows one off scripts run using for this project 
# Has a number of steps for filtering VCFs 
# Filtering steps 
# Annotate filtered VCFs for regions which have a 50% overlap with regions in DGV_Gold files
# Anotate fitlered VCFs for regions which DONT have a 50% overlap with regions in the DGV files
# Annotate the filtered VCF with regions which are present in the gene panel bed file  


# One off commands run 
#### Run Bedtools (v2.26.0) in the command line ###

# ===================== #
## Merge nearby DGV regions ## 
# ===================== #

# Sort into chromosome order before merg 
# 1) sort -k1,1 -k2,2n DGV_GS_hg38_tidy.txt > DGV_GS_hg38_tidy_sorted.txt
# Merge sections together which overlap, merge regions which have 100 bp between them
# Take DGV gold file from 10,446 to 7,735
# doing it at 50 bp makes it 7,745
# -i inoput -d Maximum distance between features allowed for features to be merged.
# 2) bedtools merge -i DGV_GS_hg38_tidy_sorted.txt -d 100 > DGV_GS_hg38_tidy_sorted_100_bp_regions_merged.txt

# ===================== #
## Merge BED file regions ## 
# ===================== #

# 1) sort -k1,1 -k2,2n genes_of_interest.bed > genes_of_interest_sorted.bed
# Merge sections together which overlap, merge regions which have 100 bp between them
# Take bedfile from 524 to 88
# 2) bedtools merge -i genes_of_interest_sorted.bed > genes_of_interest_sorted_merged.bed

 # ===================== #
##  bedtools intersects for VCFs ## 
 # ===================== #


list_of_samples=/home/erin/Documents/Work/WGS_puberty_CNV/WGS_puberty_CNVs/vcfs/sample_list.txt
path_to_filtered_vcfs=/home/erin/Documents/Work/WGS_puberty_CNV/WGS_puberty_CNVs/filtered_vcfs
path_to_dgv=/home/erin/Documents/Work/WGS_puberty_CNV/WGS_puberty_CNVs/dgv
dgv_tidy_file=DGV_GS_hg38_tidy_sorted_100_bp_regions_merged.txt
path_to_genes_bed=/home/erin/Documents/Work/WGS_puberty_CNV/WGS_puberty_CNVs/reference


while IFS="" read line || [[ -n "$line" ]];
do
        sample_ID=$(echo $line | cut -d " " -f1)
        echo $sample_ID 
        # Find regions which are in normal controls 
        # -f 	Minimum overlap required as a fraction of A. 
        # -r Require that the fraction of overlap be reciprocal for A and B. 
        # In other words, if -f is 0.50 and -r is used, this requires that B overlap at least 50% of A and that A also overlaps at least 50% of B.
        # - wa keeps original region for -a and -wb the same for -b
        bedtools intersect -a $path_to_filtered_vcfs/$sample_ID.filtered_vcf.txt -b $path_to_dgv/$dgv_tidy_file -f 0.5 -r -wa -wb >> $path_to_filtered_vcfs/dgv_overlap_$sample_ID.txt
        # Find regions which are not in normal controls 
        # -v finds all regions in -a which don't appear in -b and a 50% overlap, for both A & B -r
        bedtools intersect -a $path_to_filtered_vcfs/$sample_ID.filtered_vcf.txt -b $path_to_dgv/$dgv_tidy_file -v -f 0.5 -r >> $path_to_filtered_vcfs/no_dgv_overlap_$sample_ID.txt
        # Find filtered regions  which are in genes of interest
        # As there's multipule transcripts in the same bed file, there's lots of duplication in some of the file 
        bedtools intersect -a $path_to_genes_bed/genes_of_interest_sorted_merged.bed -b $path_to_filtered_vcfs/no_control_dgv_$sample_ID.txt  \
                              -wa -wb -f 0.1  >> $path_to_filtered_vcfs/control_filtered_genes_of_interest_little_f_10_$sample_ID.txt
      
done < $list_of_samples