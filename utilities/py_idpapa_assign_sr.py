#!/usr/bin/env python
import sys,re,time,argparse
from multiprocessing import cpu_count,Pool

def main(args):
	print >>sys.stdout, "Start analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S")
	output_gpd = args.output
	iso_list = get_iso_info(args.isoform)
	p = Pool(processes=args.cpu)
	csize = 100
	results = p.imap(func=assignment,iterable=generate_tx(args.short_reads,iso_list),chunksize=csize)
	for res in results:
		if not res: continue
		output_gpd.write(res+"\n")
	output_gpd.close()
	print >>sys.stdout, "Finish analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S")

def generate_tx(input_sr,iso_list):
	z = 0
	for line in input_sr:
		z += 1
		yield (line,z,iso_list)

# align first mate without splice alignment
def align_first_mate_s(strand,iso_exon_start,iso_exon_end,sr_exon_start,sr_exon_end):
	indic = "mismatch"
	if strand == "+":
		if int(sr_exon_start.split(",")[0]) >= int(iso_exon_start.split(",")[-2]) and int(sr_exon_end.split(",")[0]) <= int(iso_exon_end.split(",")[-2]):
			indic = "match"
		else:
			indic = "mismatch"
	else:
		if int(sr_exon_start.split(",")[0]) >= int(iso_exon_start.split(",")[0]) and int(sr_exon_end.split(",")[0]) <= int(iso_exon_end.split(",")[0]):
			indic = "match"
		else:
			indic = "mismatch"
	return indic

# align first mate with splice alignment
def align_first_mate_m(strand,iso_exon_number,iso_exon_start,iso_exon_end,sr_exon_number,sr_exon_start,sr_exon_end):
	iso_junc_list = []
	sr_junc_list = []
	indic = "mismatch"
	for i in range(0,int(iso_exon_number)-1):
		iso_junc_list.append(iso_exon_end.split(",")[i])
		iso_junc_list.append(iso_exon_start.split(",")[i+1])
	iso_junc_set = "," + ",".join(iso_junc_list) + ","
	iso_whole_set = "," + iso_exon_start.split(",")[0] + iso_junc_set + iso_exon_end.split(",")[-2] + ","
	for i in range(0,int(sr_exon_number)-1):
		sr_junc_list.append(sr_exon_end.split(",")[i])
		sr_junc_list.append(sr_exon_start.split(",")[i+1])
	sr_junc_set = "," + ",".join(sr_junc_list) + ","
	if strand == "+":
		pattern = sr_junc_set + "$"
		if int(sr_exon_end.split(",")[-2]) <= int(iso_exon_end.split(",")[-2]) and re.search(pattern,iso_junc_set) and int(sr_exon_start.split(",")[0]) >= int(iso_whole_set.split(sr_junc_set)[0].split(",")[-1]):
			indic = "match"
		else:
			indic = "mismatch"
	else:
		pattern = "^" + sr_junc_set
		if int(sr_exon_start.split(",")[0]) >= int(iso_exon_start.split(",")[0]) and re.search(pattern,iso_junc_set) and int(sr_exon_end.split(",")[-2]) <= int(iso_whole_set.split(sr_junc_set)[1].split(",")[0]):
			indic = "match"
		else:
			indic = "mismatch"
	return indic

# align second mate without splice alignment
def align_second_mate_s(iso_exon_number,iso_exon_start,iso_exon_end,sr_exon_start,sr_exon_end):
	indic = "mismatch"
	if int(iso_exon_number) == 1:
		if int(sr_exon_start.split(",")[0]) >= int(iso_exon_start.split(",")[0]) and int(sr_exon_end.split(",")[0]) <= int(iso_exon_end.split(",")[0]):
			indic = "match"
		else:
			indic = "mismatch"
	else:
		for i in range(0,int(iso_exon_number)):
			if int(sr_exon_start.split(",")[0]) >= int(iso_exon_start.split(",")[i]) and int(sr_exon_end.split(",")[0]) <= int(iso_exon_end.split(",")[i]):
				indic = "match"
				break
			else:
				indic = "mismatch"
	return indic

# align second mate with splice alignment
def align_second_mate_m(iso_exon_number,iso_exon_start,iso_exon_end,sr_exon_number,sr_exon_start,sr_exon_end):
	iso_junc_list = []
	sr_junc_list = []
	indic = "mismatch"
	for i in range(0,int(iso_exon_number)-1):
		iso_junc_list.append(iso_exon_end.split(",")[i])
		iso_junc_list.append(iso_exon_start.split(",")[i+1])
	iso_junc_set = "," + ",".join(iso_junc_list) + ","
	iso_whole_set = "," + iso_exon_start.split(",")[0] + iso_junc_set + iso_exon_end.split(",")[-2] + ","
	for i in range(0,int(sr_exon_number)-1):
		sr_junc_list.append(sr_exon_end.split(",")[i])
		sr_junc_list.append(sr_exon_start.split(",")[i+1])
	sr_junc_set = "," + ",".join(sr_junc_list) + ","
	if re.search(sr_junc_set,iso_junc_set) and len(iso_whole_set.split(sr_junc_set)[0].split(","))%2 == 0 and int(sr_exon_start.split(",")[0]) >= int(iso_whole_set.split(sr_junc_set)[0].split(",")[-1]) and int(sr_exon_end.split(",")[-2]) <= int(iso_whole_set.split(sr_junc_set)[1].split(",")[0]):
		indic = "match"
	else:
		indic = "mismatch"
	return indic

# extract pseudo isoform information
def get_iso_info(iso_gpd):
	iso_list = []
	for line in iso_gpd:
		iso_list.append(line.strip())
	return iso_list
	iso_gpd.close()

def assignment(inputs):
	(line,z,iso_list) = inputs
	read_id,chr,strand,start,end,mapq_1,sf_1,exon_number_1,exon_start_1,exon_end_1,mapq_2,sf_2,exon_number_2,exon_start_2,exon_end_2 = line.rstrip("\n").split("\t")
	sr_info = line.rstrip("\n")
	sr_polya_iso = []
	for iso in iso_list:
		gene_id,isoform_id,iso_chr,iso_strand,tss,tts,cds_start,cds_end,exon_number,exon_start,exon_end = iso.split("\t")
		if iso_chr == chr and iso_strand == strand and int(tss) <= int(start) and int(tts) >= int(end) and int(exon_number) >= int(exon_number_1) and int(exon_number) >= int(exon_number_2):
			if int(exon_number_1) == 1 and int(exon_number_2) == 1:
				indic_1 = align_first_mate_s(strand,exon_start,exon_end,exon_start_1,exon_end_1)
				indic_2 = align_second_mate_s(exon_number,exon_start,exon_end,exon_start_2,exon_end_2)
				if indic_1 == "match" and indic_2 == "match":
					sr_polya_iso.append(isoform_id)
			elif int(exon_number_1) == 1 and int(exon_number_2) > 1:
				indic_1 = align_first_mate_s(strand,exon_start,exon_end,exon_start_1,exon_end_1)
				indic_2 = align_second_mate_m(exon_number,exon_start,exon_end,exon_number_2,exon_start_2,exon_end_2)
				if indic_1 == "match" and indic_2 == "match":
					sr_polya_iso.append(isoform_id)
			elif int(exon_number_1) > 1 and int(exon_number_2) == 1:
				indic_1 = align_first_mate_m(strand,exon_number,exon_start,exon_end,exon_number_1,exon_start_1,exon_end_1)
				indic_2 = align_second_mate_s(exon_number,exon_start,exon_end,exon_start_2,exon_end_2)
				if indic_1 == "match" and indic_2 == "match":
					sr_polya_iso.append(isoform_id)
			else:
				indic_1 = align_first_mate_m(strand,exon_number,exon_start,exon_end,exon_number_1,exon_start_1,exon_end_1)
				indic_2 = align_second_mate_m(exon_number,exon_start,exon_end,exon_number_2,exon_start_2,exon_end_2)
				if indic_1 == "match" and indic_2 == "match":
					sr_polya_iso.append(isoform_id)
	if sr_polya_iso != []:
		return line.rstrip("\n") + "\t" + ",".join(sr_polya_iso)
	else:
		return None

def do_inputs():
	output_gpd_format = '''
1. read id
2. chromosome
3. strand
4. start site of alignment of fragment
5. end site of alignment of fragment
6. MAPQ of read1 (mate1) 
7. Number of nucleotides that are softly-clipped by aligner (mate1)
8. exon number (mate1)
9. exon start set (mate1)
10. exon end set (mate1)
11. MAPQ of read1 (mate2) 
12. Number of nucleotides that are softly-clipped by aligner (mate2)
13. exon number (mate2)
14. exon start set (mate2)
15. exon end set (mate2)
16. isoform set containing this polyA site'''
	parser = argparse.ArgumentParser(description="Function: assign the polyA sites identified by short reads to specific isoforms",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-r','--short_reads',type=argparse.FileType('r'),required=True,help="Short reads gpd file")
	parser.add_argument('-i','--isoform',type=argparse.FileType('r'),required=True,help="Input: isoform gpd file")
	parser.add_argument('-o','--output',type=argparse.FileType('w'),required=True,help="Output: short reads with assigned isoforms")
	parser.add_argument('-p','--cpu',type=int,default=cpu_count(),help="Number of process")
	args = parser.parse_args()
	return args

if __name__=="__main__":
	args = do_inputs()
	main(args)
