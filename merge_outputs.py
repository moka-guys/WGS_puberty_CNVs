import pandas as pd 


# This script #
# Pulls in all the various data created
# Merges it into one excel file 

column_names = ["CHROM.CNV",  "START.CNV", "END.CNV", "ID", "REF", "ALT",
        "SAMPLE", "SUPP", "SUPP_VEC", "AVGLEN", "SVTYPE", "SVMETHOD", 
         "CHR2", "CIPO", "CIEND", "STRANDS","CALLERS"]

column_names_bed_filtered = ["CHROM.BED", "START.BED","END.BED","CHROM.CNV",  
                            "START.CNV", "END.CNV", "ID", "REF", "ALT","sample",
                            "SUPP", "SUPP_VEC", "AVGLEN", "SVTYPE", "SVMETHOD", 
                            "CHR2", "CIPO", "CIEND", "STRANDS","CALLERS"]

path_to_directory = "/home/erin/Documents/Work/WGS_puberty_CNV/WGS_puberty_CNVs"


with open(path_to_directory + "/vcfs/sample_list.txt", "r") as sample_list:   
    for sample_ID in sample_list:
        # take the new line off the end of sample_ID
        sample_ID = sample_ID.strip()
        print("Working on sample: " + sample_ID)
        # Load relevant data
        vcf_unfiltered = pd.read_csv(path_to_directory + "/vcfs/" + sample_ID + 
                                    ".unfiltered_vcf.txt", 
                                    sep = "\t",
                                    names= column_names,
                                    )            

        vcf_filtered = pd.read_csv(path_to_directory + "/filtered_vcfs/" + sample_ID + 
                                    ".filtered_vcf.txt", 
                                    sep = "\t",
                                    names= column_names)

        vcf_no_dgv_overlap = pd.read_csv(path_to_directory + "/filtered_vcfs/" + sample_ID + 
                                        ".no_dgv_overlap.txt", 
                                        sep = "\t",
                                        names= column_names)

        vcf_bed_file_filtered = pd.read_csv(path_to_directory + "/filtered_vcfs/" + sample_ID +
                                             ".bed_filtered_regions.txt", 
                                            sep = "\t",
                                            names= column_names_bed_filtered) 
        writer = pd.ExcelWriter(path_to_directory + "/results/" + sample_ID + "_output.xlsx", engine="xlsxwriter")
        vcf_unfiltered.to_excel(writer, sheet_name='Unfiltered VCF', index=False)
        vcf_filtered.to_excel(writer, sheet_name='Filtered VCF', index=False)
        vcf_no_dgv_overlap.to_excel(writer, sheet_name='No DGV overlap regions', index=False)
        vcf_bed_file_filtered.to_excel(writer, sheet_name='BED file overlap regions', index=False)
        writer.save()
        print("Saved XLSX: " + sample_ID)


