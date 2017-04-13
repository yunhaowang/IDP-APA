#!/usr/bin/env python
import sys,time,argparse

def main(args):
#	print >>sys.stdout, "Start analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S")
	dic_iso_pa_set = extract_sr_pa(args.short_reads)
	merge_pa(dic_iso_pa_set,args.long_reads,args.output)
#	print >>sys.stdout, "Finish analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S")

def extract_sr_pa(sr_pa):
	dic_iso_pa = {}
	dic_iso_pa_set = {}
	for line in sr_pa: # decide if the polyA site is uniquely or multiply for specific isoform
		read_id,chr,strand,start,end,mapq_1,sf_1,exon_number_1,exon_start_1,exon_end_1,mapq_2,sf_2,exon_number_2,exon_start_2,exon_end_2,iso_set = line.rstrip("\n").split("\t")
		if iso_set != "":
			if strand == "+":
				pa = int(end)
			else:
				pa = int(start)
			if len(iso_set.split(",")) == 1:
				iso_id = iso_set
				if iso_id not in dic_iso_pa.keys():
					dic_iso_pa[iso_id] = {}
					if pa not in dic_iso_pa[iso_id].keys():
						dic_iso_pa[iso_id][pa] = {}
						dic_iso_pa[iso_id][pa]["S"] = 1
						dic_iso_pa[iso_id][pa]["M"] = 0
					else:
						dic_iso_pa[iso_id][pa]["S"] += 1
				else:
					if pa not in dic_iso_pa[iso_id].keys():
						dic_iso_pa[iso_id][pa] = {}
						dic_iso_pa[iso_id][pa]["S"] = 1
						dic_iso_pa[iso_id][pa]["M"] = 0
					else:
						dic_iso_pa[iso_id][pa]["S"] += 1
			else:
				for iso_id in iso_set.split(","):
					if iso_id not in dic_iso_pa.keys():
						dic_iso_pa[iso_id] = {}
						if pa not in dic_iso_pa[iso_id].keys():
							dic_iso_pa[iso_id][pa] = {}
							dic_iso_pa[iso_id][pa]["M"] = 1
							dic_iso_pa[iso_id][pa]["S"] = 0
						else:
							dic_iso_pa[iso_id][pa]["M"] += 1
					else:
						if pa not in dic_iso_pa[iso_id].keys():
							dic_iso_pa[iso_id][pa] = {}
							dic_iso_pa[iso_id][pa]["M"] = 1
							dic_iso_pa[iso_id][pa]["S"] = 0
						else:
							dic_iso_pa[iso_id][pa]["M"] += 1
	for iso_id in dic_iso_pa.keys(): # sort the polyA position
		pa_list = dic_iso_pa[iso_id].keys()
		pa_list.sort()
		pa_set_list = []
		for pa in pa_list:
			pa_set = str(pa) + "_" + str(dic_iso_pa[iso_id][pa]["S"]) + "S" + str(dic_iso_pa[iso_id][pa]["M"]) + "M"
			pa_set_list.append(pa_set)
		dic_iso_pa_set[iso_id] = ",".join(pa_set_list)
	return dic_iso_pa_set
	sr_pa.close()

def merge_pa(dic_iso_pa_set,lr_pa,output_gpd):
	for line in lr_pa:
		iso_id = line.rstrip("\n").split("\t")[1].split(",")[0] # for merge isoform id
		if iso_id in dic_iso_pa_set.keys():
			print >>output_gpd, line.strip() + "\t" + dic_iso_pa_set[iso_id]
		else:
			print >>output_gpd, line.strip() + "\t" + ""
	lr_pa.close()
	output_gpd.close()

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
12. polyA sites identified by long reads
13. polyA sites identified by short reads'''
	parser = argparse.ArgumentParser(description="Function: merge the polyA sites identified by short and long reads for each isoforms",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-s','--short_reads',type=argparse.FileType('r'),required=True,help="Short reads with assigned isoform")
	parser.add_argument('-l','--long_reads',type=argparse.FileType('r'),required=True,help="Isoform with polyA identified by long reads")
	parser.add_argument('-o','--output',type=argparse.FileType('w'),required=True,help="Output: isoform with polyA identified by short reads and long reads")
	args = parser.parse_args()
	return args

if __name__=="__main__":
	args = do_inputs()
	main(args)
