#!/usr/bin/env python
import sys,time,argparse

def main(args):
#	print >>sys.stdout, "Start analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S")
	dic_genome_size = extract_genome_size(args.genome_size)
	extendsion(dic_genome_size,args.length,args.input,args.output)
#	print >>sys.stdout, "Finish analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S")

def extract_genome_size(genome_size):
	dic_genome_size = {}
	for line in genome_size:
		chr,size = line.strip().split("\t")[:2]
		dic_genome_size[chr] = int(size)
	return dic_genome_size
	genome_size.close()

def replace_tss_and_tts_gpd(tss_fl,tts_fl,gpd):
	gene_id,isoform_id,chr,strand,tss,tts,cds_start,cds_end,exon_number,exon_start,exon_end = gpd.split("\t")
	tss_c = str(tss_fl)
	tts_c = str(tts_fl)
	cds_start_c = str(cds_start) 
	cds_end_c = str(cds_end)
	exon_start_c = exon_start.split(",")
	exon_start_c[0] = tss_fl
	exon_end_c = exon_end.split(",")
	exon_end_c[-2] = tts_fl
	gpd_c = "\t".join(gpd.split("\t")[:4]) + "\t" + tss_c + "\t" + tts_c + "\t" + cds_start_c + "\t" + cds_end_c + "\t" + str(exon_number) + "\t" + ",".join(exon_start_c) + "\t" + ",".join(exon_end_c)
	return gpd_c

def extendsion(dic_genome_size,extend_length,input_gpd,output_gpd):
	for line in input_gpd:
		gene_id,iso_id,chrom,strand,tss,tts,lr_c,lr_set,exon_number,exon_start,exon_end = line.rstrip("\n").split("\t")
		if (int(tss)-extend_length) >=0:
			tss_fl = str(int(tss)-extend_length)
		else:
			tss_fl = "0"
#		print line,
		if (int(tts) + extend_length) <= dic_genome_size[chrom]:
			tts_fl = str(int(tts)+extend_length)
		else:
			tts_fl = str(dic_genome_size[chrom])
#		print line,
		gpd = line.rstrip("\n")
		gpd_fl = replace_tss_and_tts_gpd(tss_fl,tts_fl,gpd)
		print >>output_gpd, gpd_fl
	input_gpd.close()
	output_gpd.close()

def do_inputs():
	parser = argparse.ArgumentParser(description="Function: extend the 5'end and 3'end of the isoform",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-l','--length',type=int,default=20,required=True,help="extend length")
	parser.add_argument('-g','--genome_size',type=argparse.FileType('r'),required=True,help="genome size tab-splitted file; first column is chromosome ID and second column is length of chromosome")
	parser.add_argument('-i','--input',type=argparse.FileType('r'),required=True,help="Input: gpd file")
	parser.add_argument('-o','--output',type=argparse.FileType('w'),required=True,help="Output: extended gpd file")
	args = parser.parse_args()
	return args

if __name__=="__main__":
	args = do_inputs()
	main(args)
