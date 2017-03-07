#!/usr/bin/env python
import sys,time,argparse

def main(args):
	print >>sys.stdout, "Start analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S")
	dic_chr_strand_en = generate_junction_set(args.input)
	output_gpd(dic_chr_strand_en,args.output)
	print >>sys.stdout, "Finish analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S")

def generate_junction_set(input_gpd):
	dic_chr_strand_en = {}
	for line in input_gpd:
		gene_id,isoform_id,chr,strand,tss,tts,cds_start,cds_end,exon_number,exon_start,exon_end = line.strip().split("\t")
		if int(exon_number) > 1:
			chr_strand_en = chr + "_" + strand + "_" + exon_number
			junction_info = "-".join(exon_start.split(",")[1:-1]) + "_" + "-".join(exon_end.split(",")[:-2])
			if chr_strand_en not in dic_chr_strand_en.keys():
				dic_chr_strand_en[chr_strand_en] = {}
				dic_chr_strand_en[chr_strand_en][junction_info] = {}
				dic_chr_strand_en[chr_strand_en][junction_info]["gene_id"] = set()
				dic_chr_strand_en[chr_strand_en][junction_info]["isoform_id"] = set()
				dic_chr_strand_en[chr_strand_en][junction_info]["tss"] = set()
				dic_chr_strand_en[chr_strand_en][junction_info]["tts"] = set()
				dic_chr_strand_en[chr_strand_en][junction_info]["gene_id"].add(gene_id)
				dic_chr_strand_en[chr_strand_en][junction_info]["isoform_id"].add(isoform_id)
				dic_chr_strand_en[chr_strand_en][junction_info]["tss"].add(int(tss))
				dic_chr_strand_en[chr_strand_en][junction_info]["tts"].add(int(tts))
			else:
				if junction_info not in dic_chr_strand_en[chr_strand_en].keys():
					dic_chr_strand_en[chr_strand_en][junction_info] = {}
					dic_chr_strand_en[chr_strand_en][junction_info]["gene_id"] = set()
					dic_chr_strand_en[chr_strand_en][junction_info]["isoform_id"] = set()
					dic_chr_strand_en[chr_strand_en][junction_info]["tss"] = set()
					dic_chr_strand_en[chr_strand_en][junction_info]["tts"] = set()
					dic_chr_strand_en[chr_strand_en][junction_info]["gene_id"].add(gene_id)
					dic_chr_strand_en[chr_strand_en][junction_info]["isoform_id"].add(isoform_id)
					dic_chr_strand_en[chr_strand_en][junction_info]["tss"].add(int(tss))
					dic_chr_strand_en[chr_strand_en][junction_info]["tts"].add(int(tts))
				else:
					dic_chr_strand_en[chr_strand_en][junction_info]["gene_id"].add(gene_id)
					dic_chr_strand_en[chr_strand_en][junction_info]["isoform_id"].add(isoform_id)
					dic_chr_strand_en[chr_strand_en][junction_info]["tss"].add(int(tss))
					dic_chr_strand_en[chr_strand_en][junction_info]["tts"].add(int(tts))
		else:
			pass
	return dic_chr_strand_en
	input_gpd.close()

def output_gpd(dic_chr_strand_en,concatenate_gpd):
	for chr_strand_en in dic_chr_strand_en.keys():
		chr,strand,en = chr_strand_en.split("_")
		for junction_info in dic_chr_strand_en[chr_strand_en].keys():
			gene_id = ",".join(dic_chr_strand_en[chr_strand_en][junction_info]["gene_id"])
			isoform_id = ",".join(dic_chr_strand_en[chr_strand_en][junction_info]["isoform_id"])
			tss = str(min(dic_chr_strand_en[chr_strand_en][junction_info]["tss"]))
			tts = str(max(dic_chr_strand_en[chr_strand_en][junction_info]["tts"]))
			gene_count = str(len(dic_chr_strand_en[chr_strand_en][junction_info]["gene_id"]))
			isoform_count = str(len(dic_chr_strand_en[chr_strand_en][junction_info]["isoform_id"]))
			exon_start = tss + "," + ",".join(junction_info.split("_")[0].split("-")) + ","
			exon_end = ",".join(junction_info.split("_")[1].split("-")) + "," + tts + ","
			print >>concatenate_gpd, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (gene_id,isoform_id,chr,strand,tss,tts,gene_count,isoform_count,en,exon_start,exon_end)
	concatenate_gpd.close()

def do_inputs():
	output_gpd_format = '''
1. gene id
2. isoform id
3. chromosome id
4. strand
5. TSS
6. TTS
7. gene count
8. isoform count
9. exon count
10. exon start set
11. exon end set'''
	parser = argparse.ArgumentParser(description="Function: concatenate multi-exon isoforms if they have the same junction conbination set.",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i','--input',type=argparse.FileType('r'),required=True,help="Input: gpd file")
	parser.add_argument('-o','--output',type=argparse.FileType('w'),required=True,help="Output: concatenated gpd file")
	args = parser.parse_args()
	return args

if __name__=="__main__":
	args = do_inputs()
	main(args)
