#!/usr/bin/env python
import sys,time,argparse
from multiprocessing import cpu_count,Pool

def main(args):
	print >>sys.stdout, "Start analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S")
	output_gpd = args.output
	dic_lr_pa = extract_polya_from_lr(args.long_reads)
	p = Pool(processes=args.cpu)
	csize = 100
	if args.primer_csv:
		flnc,nflncpa = extract_long_read_primer_info(args.primer_csv)
		results = p.imap(func=assignment_with_csv,iterable=generate_tx_with_csv(args.input,flnc,nflncpa,dic_lr_pa),chunksize=csize)
		for res in results:
			if not res: continue
			output_gpd.write(res+"\n")
		output_gpd.close()
	else:
		results = p.imap(func=assignment_without_csv,iterable=generate_tx_without_csv(args.input,dic_lr_pa),chunksize=csize)
		for res in results:
			if not res: continue
			output_gpd.write(res+"\n")
		output_gpd.close()
	print >>sys.stdout, "Finish analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S")

def generate_tx_with_csv(input_gpd,flnc,nflncpa,dic_lr_pa):
	z = 0
	for line in input_gpd:
		z += 1
		yield (line,z,flnc,nflncpa,dic_lr_pa)

def generate_tx_without_csv(input_gpd,dic_lr_pa):
	z = 0
	for line in input_gpd:
		z += 1
		yield (line,z,dic_lr_pa)

def extract_long_read_primer_info(lr_primer_csv):
	flnc = []
	nflncpa = []
	head = 1
	for line in lr_primer_csv:
		if head:
			head -= 1
		else:
			if line.strip().endswith("0"):
				id,strand,fiveseen,polyAseen,threeseen,fiveend,polyAend,threeend,primer,chimera = line.rstrip("\n").split(",")
				if fiveseen == "1" and polyAseen == "1" and threeseen == "1":
					flnc.append(id)
				elif fiveseen == "1" and polyAseen == "1" and threeseen == "0":
					nflncpa.append(id)
				elif fiveseen == "0" and polyAseen == "1" and threeseen == "1":
					nflncpa.append(id)
				else:
					pass
	return flnc,nflncpa
	lr_primer_csv.close()

def extract_polya_from_lr(lr_gpd):
	dic_lr_pa = {}
	for line in lr_gpd:
		read_id,read_id2,chr,strand,tss,tts,cds_start,cds_end,exon_number,exon_start,exon_end = line.rstrip("\n").split("\t")
		if strand == "+":
			dic_lr_pa[read_id] = int(tts)
		else:
			dic_lr_pa[read_id] = int(tss)
	return dic_lr_pa
	lr_gpd.close()

def assignment_without_csv(inputs):
	(line,z,dic_lr_pa) = inputs
	gene_id,iso_id,chrom,strand,tss,tts,lr_set,lr_c,exon_number,exon_start,exon_end = line.rstrip("\n").split("\t")
	dic_ps_iso_pa = {}
	if int(lr_c) > 1:
		for lr_id in lr_set.split(","):
			if dic_lr_pa[lr_id] not in dic_ps_iso_pa.keys():
				dic_ps_iso_pa[dic_lr_pa[lr_id]] = {}
				dic_ps_iso_pa[dic_lr_pa[lr_id]]["L"] = 1
			else:
				dic_ps_iso_pa[dic_lr_pa[lr_id]]["L"] = 1
	else:
		lr_id = lr_set
		if dic_lr_pa[lr_id] not in dic_ps_iso_pa.keys():
			dic_ps_iso_pa[dic_lr_pa[lr_id]] = {}
			dic_ps_iso_pa[dic_lr_pa[lr_id]]["L"] = 1
		else:
			dic_ps_iso_pa[dic_lr_pa[lr_id]]["L"] = 1
	pa_list = dic_ps_iso_pa.keys()
	pa_list.sort()
	pa_list_set = []
	for pa in pa_list:
		pa_set = str(pa) + "_" + str(dic_ps_iso_pa[pa]["L"]) + "L"
		pa_list_set.append(pa_set)
	return line.strip() + "\t" + ",".join(pa_list_set)

def assignment_with_csv(inputs):
	(line,z,flnc,nflncpa,dic_lr_pa) = inputs
	gene_id,iso_id,chrom,strand,tss,tts,lr_set,lr_c,exon_number,exon_start,exon_end = line.rstrip("\n").split("\t")
	dic_ps_iso_pa = {}
	if int(lr_c) > 1:
		for lr_id in lr_set.split(","):
			if dic_lr_pa[lr_id] not in dic_ps_iso_pa.keys():
				dic_ps_iso_pa[dic_lr_pa[lr_id]] = {}
				dic_ps_iso_pa[dic_lr_pa[lr_id]]["F"] = 0
				dic_ps_iso_pa[dic_lr_pa[lr_id]]["N"] = 0
				dic_ps_iso_pa[dic_lr_pa[lr_id]]["P"] = 0				
				if lr_id in flnc:
					dic_ps_iso_pa[dic_lr_pa[lr_id]]["F"] += 1
				elif lr_id in nflncpa:
					dic_ps_iso_pa[dic_lr_pa[lr_id]]["N"] += 1
				else:
					dic_ps_iso_pa[dic_lr_pa[lr_id]]["P"] += 1
			else:
				if lr_id in flnc:
					dic_ps_iso_pa[dic_lr_pa[lr_id]]["F"] += 1
				elif lr_id in nflncpa:
					dic_ps_iso_pa[dic_lr_pa[lr_id]]["N"] += 1
				else:
					dic_ps_iso_pa[dic_lr_pa[lr_id]]["P"] += 1					                  
	else:
		lr_id = lr_set
		if dic_lr_pa[lr_id] not in dic_ps_iso_pa.keys():
			dic_ps_iso_pa[dic_lr_pa[lr_id]] = {}
			dic_ps_iso_pa[dic_lr_pa[lr_id]]["F"] = 0
			dic_ps_iso_pa[dic_lr_pa[lr_id]]["N"] = 0
			dic_ps_iso_pa[dic_lr_pa[lr_id]]["P"] = 0
			if lr_id in flnc:
				dic_ps_iso_pa[dic_lr_pa[lr_id]]["F"] += 1
			elif lr_id in nflncpa:
				dic_ps_iso_pa[dic_lr_pa[lr_id]]["N"] += 1
			else:
				dic_ps_iso_pa[dic_lr_pa[lr_id]]["P"] += 1
		else:
			if lr_id in flnc:
				dic_ps_iso_pa[dic_lr_pa[lr_id]]["F"] += 1
			elif lr_id in nflncpa:
				dic_ps_iso_pa[dic_lr_pa[lr_id]]["N"] += 1
			else:
				dic_ps_iso_pa[dic_lr_pa[lr_id]]["P"] += 1
	pa_list = dic_ps_iso_pa.keys()
	pa_list.sort()
	pa_list_set = []
	for pa in pa_list:
		pa_set = str(pa) + "_" + str(dic_ps_iso_pa[pa]["F"]) + "F" + str(dic_ps_iso_pa[pa]["N"]) + "N" + str(dic_ps_iso_pa[pa]["P"]) + "P"
		pa_list_set.append(pa_set)
	return line.strip() + "\t" + ",".join(pa_list_set)

def do_inputs():
	output_gpd_format = '''
1. gene id
2. isoform id
3. chromosome id
4. strand
5. TSS
6. TTS
7. long reads set
8. long reads count
9. exon count
10. exon start set
11. exon end set
12. polyA sites identified by long reads'''
	parser = argparse.ArgumentParser(description="Function: assign the polyA sites identified by long reads to specific isoforms",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-c','--primer_csv',type=argparse.FileType('r'),help="Optional: Primer csv file produced by Iso-Seq showing the statistics of long reads")
	parser.add_argument('-r','--long_reads',type=argparse.FileType('r'),required=True,help="Long reads gpd file")
	parser.add_argument('-i','--input',type=argparse.FileType('r'),required=True,help="Input: isoform gpd file")
	parser.add_argument('-o','--output',type=argparse.FileType('w'),required=True,help="Output: isoform with polyA gpd file")
	parser.add_argument('-p','--cpu',type=int,default=cpu_count(),help="Number of process")
	args = parser.parse_args()
	return args
if __name__=="__main__":
	args = do_inputs()
	main(args)
