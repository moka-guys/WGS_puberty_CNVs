# Scripts for filtering WGS Puberty SV data

## What do these scripts do? 

### Filter WGS SV Puberty delay VCFs 
These scripts reduce the number of calls made by the Parliment SV caller by:
- Removing SV calls which did not pass QC 
- Removing SV calls which did not have >=2 callers agree on the SV
- Removing regions which appear in DGV Gold Standard variants 

### Annotate SVs which occurred in regions of interest
- Filtered VCFs are used to determine how often SVs are found in regions of interest (transcript_annotated_genes_of_interest.bed)

### Try to discover additional significant SVs in genes which aren't in the genes of interest (gene discovery)
- Take filtered VCFs & remove all regions which occur in regions of interest 
- Annotate the genes these SVs appear in 
- Create an out which shows all the gene annotated regions per patient SV occurred in 
- Creates a summary document of genes/ collections of genes which occurred in >= 20 of the samples in the study


The below information details how the scripts in the repository work

## Run tidy_reference_data.R

#### Annotate regions of interest with gene symbols 
1) Take genes_of_interest.bed (multiple transcripts per gene of interest)
2) Run biomaRt to add HGNC ID and HGNC Symbol for each transcript
Output: bed_file_with_gene.bed

#### Tidy DGV Gold standard variants 
DGV gold standard variants (taken from http://dgv.tcag.ca/dgv/docs/DGV.GS.hg38.gff3, release date 2016-05-15)
  1) Takes the thick regions for the coordinates for each variant (which has the highest confidence)
  2) Removes variants not seen at > 5.0% of the population to align with ACMG guidance for frequency of likely benign variants in the population 
  
## Run tidy_reference_data_step_2.sh
This script performs three separate tasks 

#### Merge nearby DGV Gold standard variant regions 
 1) Use bedtools merge to merge regions which are < 50 bp distance from each other

#### Merge transcript genes of interest bed file
 1) Use bedtools merge to merge regions which overlap, giving the largest transcript per gene of interest 
 
#### Merge above file with gene symbols 
 1) Use bedtools intersect to add the gene symbols on to the largest transcript for each gene
 
## Run filtering_vcf.R

#### Filter Parliment2 VCF outputs
1) Remove SVs which haven't passed quality filtering 
2) Remove SVs not called by at least two callers 

## Run bedtools_WGS_Puberty_CNV.sh

#### Intersect VCF with various other files 
  1) Find results from filtered VCF which overlap with regions in the DGV Gold file. 50% overlap of DGV gold & variant from VCF must be met 
  2) Find results from filtered VCF which do not overlap with regions in the DGV Gold file.  50% overlap of DGV gold & variant from VCF must be met
  3) Filter regions based on genes of interest bed file (transcript_annotated_genes_of_interest.bed)
  
## Run count_results.R
1) Count all the rows, per sample across the different file outputs
2) Save as a csv file

## Run merge_outputs.py
1) Pull in all the different data into dataframes 
2) Set them to output on a separate tab, of a single .xlxs sheet per sample 
3) Save!

## Run genes_of_interest_tidy.R
1) Takes a file of all samples intersected with transcript_annotated_genes_of_interest.bed
2) Count the number of samples with SVs in each region 
3) Create a dataframe of samples with SVs on a per gene basis 

## Run merge_gene_of_interest_results.py
1) Pull in all the different genes of interest data into dataframes
2) Set them to output on a separate tab of a single .xlxs sheet  
3) Save!

## Gene discovery steps

## Run gene_discovery.sh: Step one 
1) Intersect regions which aren't in DGV gold & find regions which don't appear in genes of interest (transcript_annotated_genes_of_interest.bed)

## Run gene_discovery.R: Step one
1) Prepares data to be merged, splitting on SV type (deletion or inversion) and genotype 

## Run gene_discovery.sh: Step two
There's a lot of repeated SV regions in the output file, where one SV is repeated many times with a unique ID
<br />
1) Sort the files saved in the previous step 
2) Merge regions which overlap by at least one base pair
3) Stick the different genotypes per SV type back together and sort again 

## Run gene_discovery.R: Step two
1) Search ensemble using the genomic coordinates for each merged SV, to find the gene/s associated with that region
2) Removes the regions which aren't associated with a gene 

## Run gene_discovery.sh: Step three
The step above (annotation with ensemble) removes all the information given from Parliment2, this section puts the genes from the previous step 
on to that date 
1) Use bedtools intersect to see which regions from the merged file (with the unique IDs) intersect with the genes found in the previous step

## Run gene_discovery.R: Step three
1) Merge the deletions and inversions back together 
2) Remove unnecessary columns 
3) Puts all the genes for the merged regions into one row per SV regions & the associated genes 
4) Puts all the genes (and combinations of genes per for larger SVs) into one dataframe, with the patients ID
5) Counts regions which overlap with each other 
6) Splits dataframes by count of samples the SV occurred in and saves 

## Run gene_discovery_merge_outputs.py
1) Pull in the output from the last script 
2) Save as tabs of the sample xlxs sheet


## Understanding outputs. 
### Common headers across outputs
<br />
QUAL: 
- LowQual,Description="Variant calls with this profile of supporting calls typically have a low overall precision">
- Unknown,Description="Insufficient quality evidence exists for calls of this type and support">
- Unconfirmed,Description="It was not possible to confirm this event by genotyping">
</br>
CIEND:
- Description="PE confidence interval around END"
<br />
CIPOS:
- Description="PE confidence interval around POS"
<br />
CHR2:
-Description="Chromosome for END coordinate in case of a translocation"
<br />
END:
- Description="End position of the structural variant"
<br />
AVGLEN:
- Description="Length of the SV"
<br />
SVMETHOD:
- Description="Method for generating this merged VCF file."
<br />
SVTYPE:
- Description="Type of the SV."
<br />
SUPP_VEC:
- Description="Vector of supporting samples."
<br />
SUPP:
-Description="Number of samples supporting the variant"
<br />
STRANDS:
- Description="Indicating the direction of the reads with respect to the type and breakpoint."
<br />
CALLERS:
- Description="Callers that support an ALT call at this position. To be included, the caller must have been confirmed by separate genotyping with SVTyper"
<br />
GT:
- Description="Genotype"

### Specific files

#### result_counts.txt
- Rows per different file type (denoting the filtering steps)

#### Sample_ID_output.xlsx is a per sample excel sheet made up of four tabs
Unfiltered VCF: This show all the variants in the VCF before filtering 
<br />
Filtered VCF: Variants remaining after filtering. QUAL must be PASS, >= two callers supporting SV. 
<br />
No DGV overlap regions: this is the filtered VCF from the above step, but variants remaining do not have >=50% overlap with known DGV Gold standard SV.
<br />
BED file overlap regions: Filtered VCF, but only showing regions which overlap with the transcript_annotated_genes_of_interest.bed file 

#### Genes_of_interest_results.xlsx 
- Count of SV per gene: counts the number of SVs detected per SV type, per gene 
- Gene tabs: Show SV genomic coordinate which has occurred in each gene, and other sample information 

#### sample_ID.gene.symbols.concat.txt 
- Genomic regions, with overlapping SVs (with the same genotype) merged and the gene/s on a per patient basis 

#### gene_discovery_dgv_gold.txt
- genes, count (in the sample population), SV type and IDs of the samples with this variant



## These scripts were developed by Viapath Genome Informatics