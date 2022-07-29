import os
vcf_folder = "/home/aled/Documents/210224_sasha_howard/vcfs"
DGV_BED="/home/aled/Documents/210224_sasha_howard/DGV.GS.hg38.gff3"
candidate_gene_BED = "/home/aled/Documents/210224_sasha_howard/genes_of_interest.bed"
common_regions_gain = []
common_regions_loss = []
candidate_gene_list = []
filtered_vcfs_folder="/home/aled/Documents/210224_sasha_howard/filtered_vcfs"

#TODO check string comparisons for coordinate comparisons
def filter_DGV_gff():
	with open(DGV_BED,'r') as dgv_gff_path:
		for line in dgv_gff_path.readlines():
			chr,ignore1,ignore2,start,stop = line.split("\t")[0:5]
			info=line.split("\t")[8].split(";")
			type=None
			number_of_samples = False
			number_of_studies = False
			frequency = False
			for item in info:
				if item == "variant_sub_type=Gain":
					type="Gain"
				elif item == "variant_sub_type=Loss":
					type="Loss"
				if "Number_of_unique_samples_tested=" in  item:
					number = int(item.split("=")[1].rstrip())
					if number > 30:
						number_of_samples = True
				if "num_studies=" in  item:
					number = int(item.split("=")[1].rstrip())
					if number > 2:
						number_of_studies = True
				if "Frequency=" in  item:
					number = float(item.split("=")[1].split("%")[0])
					if number > 1.0:
						frequency = True
				
				
			if frequency and number_of_samples and number_of_studies:
				if type == "Gain":
					common_regions_gain.append((chr,start,stop))
				elif type == "Loss":
					common_regions_loss.append((chr,start,stop))
	print "losses", len(common_regions_loss)
	print "gains", len(common_regions_gain)



def create_candidate_gene_list():
	with open(candidate_gene_BED,'r') as candidate_gene_BED_path:
		for line in candidate_gene_BED_path:
			chr,start,stop=line.split("\t")[0:3]
			candidate_gene_list.append((chr,start,stop))

def check_for_overlap(list,tuple):
	test_chr,test_start,test_stop=tuple
	for line in list:
		list_chr,list_start,list_stop = line
		if list_chr == test_chr and list_stop<test_start and test_stop>list_start:
			return True
	return False

def check_for_complete_inclusion(list,tuple):
	test_chr,test_start,test_stop=tuple
	for line in list:
		list_chr,list_start,list_stop = line
		if list_chr == test_chr and test_start>list_start and test_stop<list_stop:
			#print "common variant filtered"
			return True
	return False




def parse_VCFs():
	for file in os.listdir(vcf_folder):
		with open(os.path.join(vcf_folder,file),'r') as vcf:
			vcf_list=vcf.readlines()
		filter_variant_count=0
		print "input vcf has %s lines" % (len)(vcf_list)
		for line in vcf_list:
			type=None
			if line.startswith("#"):
				pass
			else:
				filter_variant_count+=1
				CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,sample = line.split("\t")
				for element in INFO.split(";"):
					if element.startswith("END="):
						stop = element.split("=")[1]
				test_tuple = (CHROM, POS, stop)
				if "DUP" in ALT:
					type="Gain"
					common = check_for_complete_inclusion(common_regions_gain,test_tuple)
				if "DEL" in ALT:
					type="Loss"
					common = check_for_complete_inclusion(common_regions_loss,test_tuple)

				if FILTER != "PASS" or common:
					vcf_list.remove(line)
					filter_variant_count - 1
		print "qual + DGV filtered vcf has %s lines" % (len)(vcf_list)
		print "filter vcf count = %s" % (str(filter_variant_count))
		yield (vcf_list)


def candidate_overlap(vcf_list):
	for line in vcf_list:
		type=None
		if line.startswith("#"):
			pass
		else:
			CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,sample = line.split("\t")
			for element in INFO.split(";"):
				if element.startswith("END="):
					stop = element.split("=")[1]
			test_tuple = (CHROM, POS, stop)
			if check_for_overlap(candidate_gene_list,test_tuple):
				vcf_list.remove(line)
				#TODO add filter flag based on overlap with DGV
	print "candidate gene filtered vcf has %s lines" % (len)(vcf_list)
				

filter_DGV_gff()
create_candidate_gene_list()
for sample in parse_VCFs():

	candidate_overlap(sample)






