#!/usr/bin/env python
import sys,time,argparse
from multiprocessing import cpu_count,Pool

def main(args):
	print >>sys.stdout, "Start analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S")
	st = args.st
	s5 = args.s5
	s3 = args.s3
	lr_file = args.input
	polished_lr_file = args.output
	dic_junction = extract_junction_from_annotation(args.anno)
	p = Pool(processes=args.cpu)
	csize = 100
	results = p.imap(func=polish,iterable=generate_tx(lr_file,dic_junction,st,s5,s3),chunksize=csize)
	for res in results:
		if not res: continue
		polished_lr_file.write(res+"\n")
	polished_lr_file.close()
	print >>sys.stdout, "Finish analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S")

def generate_tx(inf,dic_junction,spiceSite_tolerance,softClip_5,softClip_3):
	z = 0
	for line in inf:
		z += 1
		yield (line,z,dic_junction,spiceSite_tolerance,softClip_5,softClip_3)

#=== extract junction site from known gene annotation library ===
def extract_junction_from_annotation(anno_file_list):
	dic_chr_type_junction = {}
	tts_list = []
	for anno_file in anno_file_list:
		for line in anno_file:
			gene_id,isoform_id,chrom,strand,tss,tts,cds_start,cds_end,exon_number,exon_start,exon_end = line.strip().split("\t")
			tts_list.append(int(tts))
			if strand == "+" and int(exon_number) != 1:
				if chrom not in dic_chr_type_junction.keys():
					dic_chr_type_junction[chrom] = {}
					for type in ["sorted_plus5","sorted_plus3","sorted_minus5","sorted_minus3"]:
						dic_chr_type_junction[chrom][type] = set()
#				continue
				for p5 in exon_end.split(",")[:-2]:
					dic_chr_type_junction[chrom]["sorted_plus5"].add(int(p5))
				for p3 in exon_start.split(",")[1:-1]:
					dic_chr_type_junction[chrom]["sorted_plus3"].add(int(p3))
			elif strand == "-" and int(exon_number) != 1:
				if chrom not in dic_chr_type_junction.keys():
					dic_chr_type_junction[chrom] = {}
					for type in ["sorted_plus5","sorted_plus3","sorted_minus5","sorted_minus3"]:
						dic_chr_type_junction[chrom][type] = set()
#				continue
				for p5 in exon_end.split(",")[:-2]:
					dic_chr_type_junction[chrom]["sorted_minus5"].add(int(p5))
				for p3 in exon_start.split(",")[1:-1]:
					dic_chr_type_junction[chrom]["sorted_minus3"].add(int(p3))
			else:
				pass
#			anno_file.close()
	for chr in dic_chr_type_junction.keys():
		for tp in dic_chr_type_junction[chr].keys():
			dic_chr_type_junction[chr][tp] = list(dic_chr_type_junction[chr][tp])
			dic_chr_type_junction[chr][tp].append(max(tts_list))
			dic_chr_type_junction[chr][tp].sort()
	return dic_chr_type_junction

#==== get the nearest known junction site for each junction of long read ===
def getNearest(x,list,tol):
	if list != [] and len(list) >=2:
		ab = list[-1]
		y = x
		for i in list:
			if abs(x-i) < ab:
				ab = abs(x-i)
				y = i
			else:
				break
		if abs(y-x) <= int(tol):
			return y
		else:
			return x
	else:
		return x

def polish(inputs):
	(line,z,dic_junction,spiceSite_tolerance,softClip_5,softClip_3) = inputs
	gene_id,isoform_id,chrom,strand,tss,tts,mapq,sf,exon_number,exon_start,exon_end = line.rstrip("\n").split("\t")
	polished_line = ""
	if (chrom in dic_junction.keys()) and (int(exon_number)!=1): # chromosome in annotation library and exon number > 1.
		gpd_start = []
		gpd_end = []
		gpd_start.append(tss)
		if strand == "+":
			if int(sf.split("_")[0]) <= softClip_5 and int(sf.split("_")[1]) <= softClip_3:
				for i in range(1,int(exon_number)):
					p5 = int(exon_end.split(",")[i-1])
					p3 = int(exon_start.split(",")[i])
					p5_c = getNearest(p5,dic_junction[chrom]["sorted_plus5"],spiceSite_tolerance)
					p3_c = getNearest(p3,dic_junction[chrom]["sorted_plus3"],spiceSite_tolerance)
					if p3_c > p5_c:
						gpd_start.append(str(p3_c))
						gpd_end.append(str(p5_c))
					else:
						print >>sys.stderr, line,
						print >>sys.stderr, "error: p3 < p5 !!! please check the splice site tolerance or exon/intron size; suggest to use smaller tolerance"
						sys.exit()
				gpd_end.append(tts)
				gpd_start.append("")
				gpd_end.append("")
				polished_line = "\t".join(line.rstrip("\n").split("\t")[:-2]) + "\t" + ",".join(gpd_start) + "\t" + ",".join(gpd_end)
		elif strand == "-":
			if int(sf.split("_")[0]) <= softClip_3 and int(sf.split("_")[1]) <= softClip_5:
				for i in range(1,int(exon_number)):
					p5 = int(exon_end.split(",")[i-1])
					p3 = int(exon_start.split(",")[i])
					p5_c = getNearest(p5,dic_junction[chrom]["sorted_minus5"],spiceSite_tolerance)
					p3_c = getNearest(p3,dic_junction[chrom]["sorted_minus3"],spiceSite_tolerance)
					if p3_c > p5_c:
						gpd_start.append(str(p3_c))
						gpd_end.append(str(p5_c))
					else:
						print >>sys.stderr, line,
						print >>sys.stderr, "error: p3 < p5 !!! please check the splice site tolerance or exon/intron size; suggest use smaller tolerance"
						sys.exit()
				gpd_end.append(tts)
				gpd_start.append("")
				gpd_end.append("")
				polished_line = "\t".join(line.rstrip("\n").split("\t")[:-2]) + "\t" + ",".join(gpd_start) + "\t" + ",".join(gpd_end)
		else:
			pass
	else:
		if strand == "+":
			if int(sf.split("_")[0]) <= softClip_5 and int(sf.split("_")[1]) <= softClip_3:
				polished_line = line.rstrip("\n")
		elif strand == "-":
			if int(sf.split("_")[0]) <= softClip_3 and int(sf.split("_")[1]) <= softClip_5:
				polished_line = line.rstrip("\n")
		else:
			pass
	if polished_line != "":
		return polished_line
	else:
		return None

def do_inputs():
	parser = argparse.ArgumentParser(description="Polish LRs by annotation",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-t','--st',type=int,default=5,help="splice site tolerance (bp)")
	parser.add_argument('-5','--s5',type=int,default=100,help="length of soft-clipped 5'end (bp)")
	parser.add_argument('-3','--s3',type=int,default=20,help="length of soft-clipped 3'end (bp)")
	parser.add_argument('-a','--anno',type=argparse.FileType('r'),required=True,nargs="+",help="Input: annotation gpd file")
	parser.add_argument('-i','--input',type=argparse.FileType('r'),required=True,help="Input: LRs gpd file")
	parser.add_argument('-o','--output',type=argparse.FileType('w'),required=True,help="Output: polished LRs gpd file")
	parser.add_argument('-p','--cpu',type=int,default=cpu_count(),help="Number of threads")
	args = parser.parse_args()
	return args

if __name__=="__main__":
	args = do_inputs()
	main(args)
