#!/bin/bash
# bedtools v2.26.0

path_to_filtered_vcfs=/home/erin/Documents/Work/WGS_puberty_CNV/WGS_puberty_CNVs/filtered_vcfs
path_to_genes_bed=/home/erin/Documents/Work/WGS_puberty_CNV/WGS_puberty_CNVs/reference
path_to_annotated_sample_beds=/home/erin/Documents/Work/WGS_puberty_CNV/WGS_puberty_CNVs/gene_discovery/annotated
path_to_vcfs_gene_discovery=/home/erin/Documents/Work/WGS_puberty_CNV/WGS_puberty_CNVs/gene_discovery/
list_of_samples=/home/erin/Documents/Work/WGS_puberty_CNV/WGS_puberty_CNVs/vcfs/sample_list.txt
path_to_merge=/home/erin/Documents/Work/WGS_puberty_CNV/WGS_puberty_CNVs/gene_discovery/to_merge/


while IFS="" read line || [[ -n "$line" ]];
do
  sample_ID=$(echo $line | cut -d " " -f1)
  echo $sample_ID 
  
  ##==============================================
  ## STEP ONE ##
  
   # -v finds all regions in -a which don't appear in -b, with at least one bp overlap
  # bedtools intersect -a $path_to_filtered_vcfs/$sample_ID.no_dgv_gold_overlap.txt -b $path_to_genes_bed/transcript_annotated_genes_of_interest.bed \
 #                       -v > $path_to_vcfs_gene_discovery/$sample_ID.no_bed_or_dgv_overlap.txt
  
  ##==============================================
  ## STEP TWO ##
  
  # Put the gene name regions back on the no_dgv_overlap regions, should remove regions from the bed file?
  # Run this step before 
  # There's a lot of repeated CNV regions in the output file, where one CNV is repeated many times with a unique ID
  # This command merges the close regions together, to reduce this overlapping, but does remove additional information 
  # Regions have to overlap by at least one bp to be merged
  # Need to split the files per patient into inversion, deletions,  etc 
  # Second R rode stage
  #sort -k1,1 -k2,2n $path_to_merge/$sample_ID.DELS.HETEROZYGOUS.txt > $path_to_merge/$sample_ID.DELS.HETEROZYGOUS.sorted.txt
  #sort -k1,1 -k2,2n $path_to_merge/$sample_ID.DELS.HOMOZYGOUS.txt > $path_to_merge/$sample_ID.DELS.HOMOZYGOUS.sorted.txt
  
  #sort -k1,1 -k2,2n $path_to_merge/$sample_ID.INVS.HETEROZYGOUS.txt > $path_to_merge/$sample_ID.INVS.HETEROZYGOUS.sorted.txt
  #sort -k1,1 -k2,2n $path_to_merge/$sample_ID.INVS.HOMOZYGOUS.txt > $path_to_merge/$sample_ID.INVS.HOMOZYGOUS.sorted.txt
  
  
  #bedtools merge -i $path_to_merge/$sample_ID.DELS.HETEROZYGOUS.sorted.txt -c 4 -o collapse  > $path_to_merge/$sample_ID.DELS.HETEROZYGOUS.merged.txt
  #bedtools merge -i $path_to_merge/$sample_ID.DELS.HOMOZYGOUS.sorted.txt -c 4 -o collapse  > $path_to_merge/$sample_ID.DELS.HOMOZYGOUS.merged.txt
  
  #bedtools merge -i $path_to_merge/$sample_ID.INVS.HETEROZYGOUS.sorted.txt -c 4 -o collapse  > $path_to_merge/$sample_ID.INVS.HETEROZYGOUS.merged.txt
  #bedtools merge -i $path_to_merge/$sample_ID.INVS.HOMOZYGOUS.sorted.txt -c 4 -o collapse  > $path_to_merge/$sample_ID.INVS.HOMOZYGOUS.merged.txt
  
  # stick them back together as one file, and sort by genomic location
  
#  cat $path_to_merge/$sample_ID.DELS.HETEROZYGOUS.merged.txt $path_to_merge/$sample_ID.DELS.HOMOZYGOUS.merged.txt | \
#  sort -k1,1 -k2,2n > $path_to_merge/$sample_ID.dels.merged.txt
  
#  cat $path_to_merge/$sample_ID.INVS.HETEROZYGOUS.merged.txt $path_to_merge/$sample_ID.INVS.HOMOZYGOUS.merged.txt | \
#  sort -k1,1 -k2,2n > $path_to_merge/$sample_ID.invs.merged.txt
  
  
 ##==============================================
 ## STEP THREE ##
 
  # Add the genes back on to the merged outputs, to annotate the orginal data 
  
bedtools intersect -a $path_to_merge/$sample_ID.dels.merged.txt -b $path_to_annotated_sample_beds/$sample_ID.DELS.gene.bed \
    -wa -wb > $path_to_annotated_sample_beds/$sample_ID.gene_bed_no_overlap_dels.txt
    
 bedtools intersect -a $path_to_merge/$sample_ID.invs.merged.txt -b $path_to_annotated_sample_beds/$sample_ID.INVS.gene.bed \
    -wa -wb > $path_to_annotated_sample_beds/$sample_ID.gene_bed_no_overlap_invs.txt
  
##============================================== 

# Line below needs to alway be uncommented to make the loop run !
done < $list_of_samples

