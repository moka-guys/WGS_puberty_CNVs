# Scripts for filtering WGS Puberty CNV data

## What do these scripts do? 

Filter WGS CNV Puberty delay VCFs
Annotate CNVs which occurred in regions of interest 
Try to discover additional significant CNVs in genes which aren't in the genes of interest

## Order for running 

## Filtering steps
## Run tidy_reference_data.R

#### Annotate regions of interest with gene symbols 
1) Take genes_of_interest.bed (multiple transcripts per gene of interest)
2) Run biomaRt to add HGNC ID and HGNC Symbol for each transcript
Output: bed_file_with_gene.bed


#### Tidy DGV Gold standard variants 
  1) DGV gold standard variants (taken from http://dgv.tcag.ca/dgv/docs/DGV.GS.hg38.gff3, release date 2016-05-15)
  2) Takes the thick regions for the coordinates for each variant (which has the highest confidence)
  3) Removes variants not seen at > 5.0% of the population to align with ACMG guidance for frequency of likely benign variants in the population 
  
## Run tidy_reference_data_step_2.sh

#### Merge nearby DGV Gold standard variant regions 
 1) Use bedtools merge to merge regions which are < 50 bp distance from each other

#### Merge transcript genes of interest bed file
 1) Use bedtools merge to merge regions which overlap, giving the largest transcript per gene of interest 
 
#### Merge above file with gene symbols 
 1) Use bedtools intersect to add the gene symbols on to the largest transcript for each gene
 
## Run filtering_vcf.R

#### Filter Parliment2 VCF outputs
1) Remove SV which hasn't pass quality filtering 
2) Remove SV not call by at least callers 

## Run bedtools_WGS_Puberty_CNV.sh

#### Intersect VCF with various other files 
  1) Find results from filtered VCF which overlap with regions in the DGV Gold file. 50% overlap of DGV gold & variant from VCF must be met 
  2) Find results from filtered VCF which do not overlap with regions in the DGV Gold file. < 50% overlap of DGV gold & variant from VCF must be met
  3) Filter regions based on genes of interest bed file (transcript_annotated_genes_of_interest.bed)
  
## Run count_results.R
1) Count all the rows, per sample across the different file outputs
2) Save as a csv file

## Run merge_outputs.py
1) Pull in all the different data into dataframes 
2) Set them to output on a separate tab, of a single .xlxs sheet per sample 
3) Save!

## Run genes_of_interest_tidy.R
1) Take file of all samples intersected with transcript_annotated_genes_of_interest.bed
2) Count number of samples with CNVs in each region 
3) Create a df of samples with CNVs on a per gene basis 

## Run merge_gene_of_interest_results.py
1) Pull in all the different genes of interest data into dataframes
2) Set them to output on a separate tab, of a single .xlxs sheet  
3) Save!

## Gene discovery steps

## Run gene_discovery.sh step one 
1) Intersect regions which aren't in DGV gold & find regions which don't appear in genes of interest (transcript_annotated_genes_of_interest.bed)

## Run gene_discovery.R 
1) Prepares data to be merged, splitting on SV type and genotype 
   Lots of the regions called in the VCFS are overlapping, this reduces the number 

## Understanding outputs. 
### Common headers across outputs
QUAL: 
- LowQual,Description="Variant calls with this profile of supporting calls typically have a low overall precision">
- Unknown,Description="Insufficient quality evidence exists for calls of this type and support">
- Unconfirmed,Description="It was not possible to confirm this event by genotyping">
CIEND:
- Description="PE confidence interval around END""
CIPOS:
- Description="PE confidence interval around POS"
CHR2:
-Description="Chromosome for END coordinate in case of a translocation"
END:
- Description="End position of the structural variant"
AVGLEN:
- Description="Length of the SV"
SVMETHOD:
- Description="Method for generating this merged VCF file."
SVTYPE:
- Description="Type of the SV."
SUPP_VEC:
- Description="Vector of supporting samples."
SUPP:
-Description="Number of samples supporting the variant"
STRANDS:
- Description="Indicating the direction of the reads with respect to the type and breakpoint."
CALLERS:
- Description="Callers that support an ALT call at this position. To be included, the caller must have been confirmed by separate genotyping with SVTyper"
GT:
- Description="Genotype"

### Specific files
#### Sample_ID_output.xlsx is a per sample excel sheet made up of four tabs
Unfiltered VCF: This show all the variants in the VCF before filtering 
Filtered VCF: Variants remaining after filtering. QUAL must be PASS, >= two callers supporting SV. 
No DGV overlap regions: this is the filtered VCF from the above step, but variants remaining do not have >=50% overlap with known DGV Gold standard SV.
BED file overlap regions: Filtered VCF, but only showing regions which overlap with the transcript_annotated_genes_of_interest.bed file 

#### Genes_of_interest_results.xlsx is a 

## These scripts were developed by Viapath Genome Informatic