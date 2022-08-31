import pandas as pd 

# This script
# 1) Pulls in all the various data created 
# 2) Merges it into one excel file 

column_names = ["GENE/S", "COUNT", "SV_TYPE", "IDs"]


path_to_directory = "/home/erin/Documents/Work/WGS_puberty_CNV/WGS_puberty_CNVs/gene_discovery/results"

df_one = pd.read_csv(path_to_directory + "/gene_discovery_dgv_gold_30_40.txt", 
                            sep = "\t")            

df_two = pd.read_csv(path_to_directory + "/gene_discovery_dgv_gold_20_29.txt", 
                            sep = "\t")     

df_three = pd.read_csv(path_to_directory + "/gene_discovery_dgv_gold_10_19.txt", 
                            sep = "\t")    


df_four = pd.read_csv(path_to_directory + "/gene_discovery_dgv_gold_4_9.txt", 
                            sep = "\t")     

writer = pd.ExcelWriter(path_to_directory + "/gene_discovery_output.xlsx", engine="xlsxwriter")
df_four.to_excel(writer, sheet_name='4 to 9 samples', index=False)
df_three.to_excel(writer, sheet_name='10 - 19 samples', index=False)
df_two.to_excel(writer, sheet_name='20 - 29 samples', index=False)
df_one.to_excel(writer, sheet_name='30 - 40 samples', index=False)
writer.save()
