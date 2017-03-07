#!/usr/bin/env python
import sys,time,argparse

def main(args):
	print >>sys.stdout, "Start analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S")
	dic_iso_info = extract_iso_info(args.input)
	output_gpd(dic_iso_info,args.output,args.add_prefix)
	print >>sys.stdout, "Finish analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S")
	
def extract_iso_info(gtf_file): # extract the all information of the isoform
	dic_iso_info = {}
	flag = 0
	for line in gtf_file:
		if (line[0] != "#") and (line.strip().split("\t")[2] == "exon" or line.strip().split("\t")[2] == "CDS"):
			chr,source,feature,start,end,score,strand,frame,group = line.rstrip("\n").split("\t")
			gene_id = group.split("gene_id \"")[1].split("\";")[0]
			transcript_id = group.split("transcript_id \"")[1].split("\";")[0]
			chr_strand = chr + "&" + strand
			if chr_strand not in dic_iso_info.keys():
				dic_iso_info[chr_strand] = {}
				dic_iso_info[chr_strand][transcript_id] = {}
				dic_iso_info[chr_strand][transcript_id]["gene_id"] = gene_id
				dic_iso_info[chr_strand][transcript_id]["exon_start"] = []
				dic_iso_info[chr_strand][transcript_id]["exon_end"] = []
				dic_iso_info[chr_strand][transcript_id]["cds_start"] = []
				dic_iso_info[chr_strand][transcript_id]["cds_end"] = []
				if feature == "exon":
					dic_iso_info[chr_strand][transcript_id]["exon_start"].append(int(start)-1)
					dic_iso_info[chr_strand][transcript_id]["exon_end"].append(int(end))
				else:
					dic_iso_info[chr_strand][transcript_id]["cds_start"].append(int(start)-1)
					dic_iso_info[chr_strand][transcript_id]["cds_end"].append(int(end))
			else:
				if transcript_id not in dic_iso_info[chr_strand].keys():
					dic_iso_info[chr_strand][transcript_id] = {}
					dic_iso_info[chr_strand][transcript_id]["gene_id"] = gene_id
					dic_iso_info[chr_strand][transcript_id]["exon_start"] = []
					dic_iso_info[chr_strand][transcript_id]["exon_end"] = []
					dic_iso_info[chr_strand][transcript_id]["cds_start"] = []
					dic_iso_info[chr_strand][transcript_id]["cds_end"] = []
				else:
					flag = 1
				if feature == "exon":
					dic_iso_info[chr_strand][transcript_id]["exon_start"].append(int(start)-1)
					dic_iso_info[chr_strand][transcript_id]["exon_end"].append(int(end))
				else:
					dic_iso_info[chr_strand][transcript_id]["cds_start"].append(int(start)-1)
					dic_iso_info[chr_strand][transcript_id]["cds_end"].append(int(end))
	return dic_iso_info
	gtf_file.close()

def output_gpd(dic_iso_info,gpd_file,add_flag): # print each isoform by gpd format
	for chr_strand in dic_iso_info.keys():
		for iso in dic_iso_info[chr_strand].keys():
			exon_start_list = dic_iso_info[chr_strand][iso]["exon_start"]
			exon_end_list = dic_iso_info[chr_strand][iso]["exon_end"]
			cds_start_list = dic_iso_info[chr_strand][iso]["cds_start"]
			cds_end_list = dic_iso_info[chr_strand][iso]["cds_end"]
			exon_start_list.sort()
			exon_end_list.sort()
			exon_number = str(len(exon_start_list))
			exon_start = str(exon_start_list[0])
			exon_end = str(exon_end_list[-1])
			if dic_iso_info[chr_strand][iso]["cds_start"] != [] and dic_iso_info[chr_strand][iso]["cds_end"] != []:
				cds_start = str(min(cds_start_list))
				cds_end = str(max(cds_end_list))
			else:
				cds_start = "."
				cds_end = "."
			if add_flag == "yes":
				print >>gpd_file, "refgene_%s\trefiso_%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (dic_iso_info[chr_strand][iso]["gene_id"],iso,chr_strand.split("&")[0],chr_strand.split("&")[1],exon_start,exon_end,cds_start,cds_end,exon_number,(",".join(str(x) for x in exon_start_list) + ","),( ",".join(str(x) for x in exon_end_list) + ","))
			else:
				print >>gpd_file, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (dic_iso_info[chr_strand][iso]["gene_id"],iso,chr_strand.split("&")[0],chr_strand.split("&")[1],exon_start,exon_end,cds_start,cds_end,exon_number,(",".join(str(x) for x in exon_start_list) + ","),( ",".join(str(x) for x in exon_end_list) + ","))
	gpd_file.close()

def do_inputs():
	output_gpd_format = '''
1. gene id
2. isoform id
3. chromosome id
4. strand
5. TSS
6. TTS
7. start site of CDS
8. end site of CDS
9. exon count
10. exon start set
11. exon end set'''
	parser = argparse.ArgumentParser(description="Convert gtf to gpd. Note: suggest to use Ensembl or GENCODE because RefSeq contains some isoforms with same ID but different locations.",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-p','--add_prefix',choices=["yes","no"],default="yes",help="If add the prefix 'refgene_' and 'refiso_' to gene ID and isoform ID, respectively")
	parser.add_argument('-i','--input',type=argparse.FileType('r'),required=True,help="Input: gtf file")
	parser.add_argument('-o','--output',type=argparse.FileType('w'),required=True,help="Output: gpd file")
	args = parser.parse_args()
	return args

if __name__=="__main__":
	args = do_inputs()
	main(args)
