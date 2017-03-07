#!/usr/bin/env python
import sys,re,time,argparse

def main(args):
	print >>sys.stdout, "Start analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S")
	if args.primer_csv == "yes":
		optimization_with_csv(args.input,args.output,args.sr_count,args.lr_count,args.distance)
	else:
		optimization_without_csv(args.input,args.output,args.sr_count,args.lr_count,args.distance)
	print >>sys.stdout, "Finish analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S")

def optimization_with_csv(input_gpd,output_gpd,sr_c,lr_c,distance):
	for line in input_gpd:
		gene_id,isoform_id,chr,strand,tss,tts,cds_start,cds_end,exon_number,exon_start,exon_end,lr_pa_set,sr_pa_set = line.rstrip("\n").split("\t")
		if sr_pa_set == "":
			print >>output_gpd, line.rstrip("\n") + "\t\tNA"
	
		else:
			if sr_pa_set.split(",") == 1:
				sr_pa = sr_pa_set
				sr_pa_pos = int(re.findall(r"^(\d+)_",sr_pa)[0])
				sr_pa_pos_s = int(re.findall(r"(\d+)S",sr_pa)[0])
				sr_pa_pos_m = int(re.findall(r"(\d+)M",sr_pa)[0])
				lr_pa_pos_c = {}
				lr_pa_pos_c["F"] = 0
				lr_pa_pos_c["N"] = 0
				for lr_pa in lr_pa_set.split(","):
					lr_pa_pos = int(re.findall(r"^(\d+)_",lr_pa)[0])
					if lr_pa_pos >= (sr_pa_pos - distance) and lr_pa_pos <= (sr_pa_pos + distance):
						lr_pa_pos_c["F"] += int(re.findall(r"(\d+)F",lr_pa)[0])
						lr_pa_pos_c["N"] += int(re.findall(r"(\d+)N",lr_pa)[0])
				if (sr_pa_pos_s + sr_pa_pos_m) >= sr_c and (lr_pa_pos_c["F"]+lr_pa_pos_c["N"]) >= lr_c:
					print >>output_gpd, line.strip() + "\t" + sr_pa + "_" + str(lr_pa_pos_c["F"]) + "F" + str(lr_pa_pos_c["N"]) + "N;" + "\tPA"
				else:
					print >>output_gpd, line.strip() + "\t" + ";" + sr_pa + "_" + str(lr_pa_pos_c["F"]) + "F" + str(lr_pa_pos_c["N"]) + "N" + "\tNA"
			else:
				dic_sr_lr_positive = {}
				dic_sr_lr_negative = {}
				sr_lr_positive_list = []
				sr_lr_negative_list = []
				sr_lr_positive_list_set = []
				sr_lr_negative_list_set = []
				for sr_pa in sr_pa_set.split(","):
					sr_pa_pos = int(re.findall(r"^(\d+)_",sr_pa)[0])
					sr_pa_pos_s = int(re.findall(r"(\d+)S",sr_pa)[0])
					sr_pa_pos_m = int(re.findall(r"(\d+)M",sr_pa)[0])
					lr_pa_pos_c = {}
					lr_pa_pos_c["F"] = 0
					lr_pa_pos_c["N"] = 0
					for lr_pa in lr_pa_set.split(","):
						lr_pa_pos = int(re.findall(r"^(\d+)_",lr_pa)[0])
						if lr_pa_pos >= (sr_pa_pos - distance) and lr_pa_pos <= (sr_pa_pos + distance):
							lr_pa_pos_c["F"] += int(re.findall(r"(\d+)F",lr_pa)[0])
							lr_pa_pos_c["N"] += int(re.findall(r"(\d+)N",lr_pa)[0])
					if (sr_pa_pos_s + sr_pa_pos_m) >= sr_c and (lr_pa_pos_c["F"]+lr_pa_pos_c["N"]) >= lr_c:
						dic_sr_lr_positive[sr_pa_pos] = sr_pa+"_"+str(lr_pa_pos_c["F"]) + "F" + str(lr_pa_pos_c["N"]) + "N"
						sr_lr_positive_list.append(sr_pa_pos)
					else:
						dic_sr_lr_negative[sr_pa_pos] = sr_pa+"_"+str(lr_pa_pos_c["F"]) + "F" + str(lr_pa_pos_c["N"]) + "N"
						sr_lr_negative_list.append(sr_pa_pos)
				for sr_lr_positive in sr_lr_positive_list:
					sr_lr_positive_list_set.append(dic_sr_lr_positive[sr_lr_positive])
				for sr_lr_negative in sr_lr_negative_list:
					sr_lr_negative_list_set.append(dic_sr_lr_negative[sr_lr_negative])
				if len(sr_lr_positive_list_set) > 1:
					print >>output_gpd, line.strip() + "\t" + ",".join(sr_lr_positive_list_set) + ";" + ",".join(sr_lr_negative_list_set) + "\tAPA"
				elif len(sr_lr_positive_list_set) == 1:
					print >>output_gpd, line.strip() + "\t" + ",".join(sr_lr_positive_list_set) + ";" + ",".join(sr_lr_negative_list_set) + "\tPA"
				else:
					print >>output_gpd, line.strip() + "\t" + ",".join(sr_lr_positive_list_set) + ";" + ",".join(sr_lr_negative_list_set) + "\tNA"
	input_gpd.close()
	output_gpd.close()

def optimization_without_csv(input_gpd,output_gpd,sr_c,lr_c,distance):
	for line in input_gpd:
		gene_id,isoform_id,chr,strand,tss,tts,cds_start,cds_end,exon_number,exon_start,exon_end,lr_pa_set,sr_pa_set = line.rstrip("\n").split("\t")
		if sr_pa_set == "":
			print >>output_gpd, line.rstrip("\n") + "\t\tNA"
		else:
			if sr_pa_set.split(",") == 1:
				sr_pa = sr_pa_set
				sr_pa_pos = int(re.findall(r"^(\d+)_",sr_pa)[0])
				sr_pa_pos_l = int(re.findall(r"(\d+)L",sr_pa)[0])
				lr_pa_pos_c = {}
				lr_pa_pos_c["L"] = 0
				for lr_pa in lr_pa_set.split(","):
					lr_pa_pos = int(re.findall(r"^(\d+)_",lr_pa)[0])
					if lr_pa_pos >= (sr_pa_pos - distance) and lr_pa_pos <= (sr_pa_pos + distance):
						lr_pa_pos_c["L"] += int(re.findall(r"(\d+)L",lr_pa)[0])
				if (sr_pa_pos_s + sr_pa_pos_m) >= sr_c and (lr_pa_pos_c["L"]) >= lr_c:
					print >>output_gpd, line.strip() + "\t" + sr_pa + "_" + str(lr_pa_pos_c["L"]) + "L;" + "\tPA"
				else:
					print >>output_gpd, line.strip() + "\t" + ";" + sr_pa + "_" + str(lr_pa_pos_c["L"]) + "L" + "\tNA"
			else:
				dic_sr_lr_positive = {}
				dic_sr_lr_negative = {}
				sr_lr_positive_list = []
				sr_lr_negative_list = []
				sr_lr_positive_list_set = []
				sr_lr_negative_list_set = []
				for sr_pa in sr_pa_set.split(","):
					sr_pa_pos = int(re.findall(r"^(\d+)_",sr_pa)[0])
					sr_pa_pos_s = int(re.findall(r"(\d+)S",sr_pa)[0])
					sr_pa_pos_m = int(re.findall(r"(\d+)M",sr_pa)[0])
					lr_pa_pos_c = {}
					lr_pa_pos_c["L"] = 0
					for lr_pa in lr_pa_set.split(","):
						lr_pa_pos = int(re.findall(r"^(\d+)_",lr_pa)[0])
						if lr_pa_pos >= (sr_pa_pos - distance) and lr_pa_pos <= (sr_pa_pos + distance):
							lr_pa_pos_c["L"] += int(re.findall(r"(\d+)L",lr_pa)[0])
					if (sr_pa_pos_s + sr_pa_pos_m) >= sr_c and (lr_pa_pos_c["L"]) >= lr_c:
						dic_sr_lr_positive[sr_pa_pos] = sr_pa+"_"+str(lr_pa_pos_c["L"]) + "L"
						sr_lr_positive_list.append(sr_pa_pos)
					else:
						dic_sr_lr_negative[sr_pa_pos] = sr_pa+"_"+str(lr_pa_pos_c["L"]) + "L"
						sr_lr_negative_list.append(sr_pa_pos)
				for sr_lr_positive in sr_lr_positive_list:
					sr_lr_positive_list_set.append(dic_sr_lr_positive[sr_lr_positive])
				for sr_lr_negative in sr_lr_negative_list:
					sr_lr_negative_list_set.append(dic_sr_lr_negative[sr_lr_negative])
				if len(sr_lr_positive_list_set) > 1:
					print >>output_gpd, line.strip() + "\t" + ",".join(sr_lr_positive_list_set) + ";" + ",".join(sr_lr_negative_list_set) + "\tAPA"
				elif len(sr_lr_positive_list_set) == 1:
					print >>output_gpd, line.strip() + "\t" + ",".join(sr_lr_positive_list_set) + ";" + ",".join(sr_lr_negative_list_set) + "\tPA"
				else:
					print >>output_gpd, line.strip() + "\t" + ",".join(sr_lr_positive_list_set) + ";" + ",".join(sr_lr_negative_list_set) + "\tNA"
	input_gpd.close()
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
13. polyA sites identified by short reads
14. polyA sites merged by long reads and short reads
15. polyA type: (1) NA, no avaiable polyA site; (2) PA, one polyA site; (3) APA, multiple polyA sites'''
	parser = argparse.ArgumentParser(description="Function: optimize the polyA sites for each isoforms",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-c','--primer_csv',choices=["yes","no"],required=True,help="If the primer csv is used for long reads polyA identification, polyA sites tagged with 'F/N/P'; or not, polyA sites tagged with 'L'.")
	parser.add_argument('-s','--sr_count',type=int,default=1,help="Number of supported short reads")
	parser.add_argument('-l','--lr_count',type=int,default=1,help="Number of supported long reads")
	parser.add_argument('-d','--distance',type=int,default=10,help="Distance between short reads-defined polyA site and long reads-defined one")
	parser.add_argument('-i','--input',type=argparse.FileType('r'),required=True,help="Input: isoform with merged polyA sites")
	parser.add_argument('-o','--output',type=argparse.FileType('w'),required=True,help="Output: isoform with optimzed polyA sites")
	args = parser.parse_args()
	return args

if __name__=="__main__":
	args = do_inputs()
	main(args)	
