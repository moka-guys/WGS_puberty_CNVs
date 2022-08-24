#!/bin/bash
# bedtools v2.26.0

# This script
# Has a number of steps for filtering VCFs 
# Filtering steps 
# Annotate filtered VCFs for regions which have a 50% overlap with regions in DGV_Gold files
# Anotate fitlered VCFs for regions which DONT have a 50% overlap with regions in the DGV files
# Annotate the filtered VCF with regions which are present in the gene panel bed file  
 # ===================== #
##  bedtools intersects for VCFs ## 
 # ===================== #
# Paths
list_of_samples=/home/erin/Documents/Work/WGS_puberty_CNV/WGS_puberty_CNVs/vcfs/sample_list.txt
path_to_filtered_vcfs=/home/erin/Documents/Work/WGS_puberty_CNV/WGS_puberty_CNVs/filtered_vcfs
path_to_reference=/home/erin/Documents/Work/WGS_puberty_CNV/WGS_puberty_CNVs/reference
dgv_tidy_file=DGV_GS_hg38_tidy_sorted_50_bp_regions_merged.txt

 while IFS="" read line || [[ -n "$line" ]];
do
        sample_ID=$(echo $line | cut -d " " -f1)
        echo $sample_ID
        # Find regions which are in DGV Gold standard variants
        # -f 	Minimum overlap required as a fraction of A.
        # -r Require that the fraction of overlap be reciprocal for A and B.
        # In other words, if -f is 0.50 and -r is used, this requires that B overlap at least 50% of A and that A also overlaps at least 50% of B.
        # - wa keeps original region for -a and -wb the same for -b
        bedtools intersect -a $path_to_filtered_vcfs/$sample_ID.filtered_vcf.txt -b $path_to_reference/$dgv_tidy_file -f 0.5 -r  > $path_to_filtered_vcfs/$sample_ID.dgv_gold_overlap.txt
        # Find regions which are not in DGV Gold standard variants
        # -v finds all regions in -a which don't appear in -b and a 50% overlap, for both A & B -r
        bedtools intersect -a $path_to_filtered_vcfs/$sample_ID.filtered_vcf.txt -b $path_to_reference/$dgv_tidy_file -v -f 0.5 -r > $path_to_filtered_vcfs/$sample_ID.no_dgv_gold_overlap.txt
        # Find filtered regions  which are in genes of interest
        bedtools intersect -a $path_to_reference/transcript_annotated_genes_of_interest.bed -b $path_to_filtered_vcfs/$sample_ID.filtered_vcf.txt  \
                              -wa -wb -f 0.1 > $path_to_filtered_vcfs/$sample_ID.bed_filtered_regions_unsorted_dgv_gold.txt


done < $list_of_samples

