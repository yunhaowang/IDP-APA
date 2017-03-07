#!/usr/bin/env python
import sys,re,time,argparse

def main(args):
	print >>sys.stdout, "Start analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S")
	output_gpd(args.input,args.output,args.length,args.ratio)
	print >>sys.stdout, "Finish analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S")

def extract_exon_length_from_cigar(cigar):
	cigar_m = ["0"] + re.findall(r"(\d+)M",cigar)
	cigar_d = ["0"] + re.findall(r"(\d+)D",cigar)
	cigar_m_s,cigar_d_s = [0,0]
	for m in cigar_m:
		cigar_m_s += int(m)
	for d in cigar_d:
		cigar_d_s += int(d)
	exon_length = cigar_m_s+cigar_d_s
	return exon_length

def extract_soft_clip_from_cigar(cigar):
	cigar_5 = ["0"] + re.findall(r"^(\d+)S",cigar)
	cigar_3 = ["0"] + re.findall(r"(\d+)S$",cigar)
	cigar_5_s,cigar_3_s = [0,0]
	for s5 in cigar_5:
		cigar_5_s += int(s5)
	for s3 in cigar_3:
		cigar_3_s += int(s3)
	return cigar_5_s,cigar_3_s

def sam2gpd(sam):
	qname,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq = sam.split("\t")[:10]
	tag = "\t".join(sam.strip().split("\t")[11:])
	s5,s3 = extract_soft_clip_from_cigar(cigar)
	sf = str(s5)+"_"+str(s3)
	strand = (re.search(r"XS:A:(\S)",tag)).group(1)
	cigar_n_l = 0
	exon_length = 0
	exon_start = int(pos)-1
	exon_end = 0
	exon_start_list = []
	exon_end_list = []
	if "N" in cigar:
		for exon in cigar.split("N"):
			exon = exon + "N"
			exon_start = exon_start + exon_length + cigar_n_l
			exon_length = extract_exon_length_from_cigar(exon)
			exon_end = exon_start + exon_length
			if re.search(r"(\d+)N",exon):
				cigar_n_l = int((re.search(r"(\d+)N",exon)).group(1))
			exon_start_list.append(str(exon_start))
			exon_end_list.append(str(exon_end))
	else:
		exon_start = exon_start
		exon_length = extract_exon_length_from_cigar(cigar)
		exon_end = exon_start + exon_length
		exon_start_list.append(str(exon_start))
		exon_end_list.append(str(exon_end))
	exon_start_list.append("")
	exon_end_list.append("")
	gpd = [qname,qname,rname,strand,str(int(pos)-1),str(exon_end),mapq,sf,str(len(exon_start_list)-1),",".join(exon_start_list),",".join(exon_end_list)]
	return gpd

read1_flag = "no_polya"
def output_gpd(sam_file,gpd_file,polya_len,a_ratio):
	global read1_flag
	for line in sam_file:
		if line[0] != "@":
			qname,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq = line.strip().split("\t")[:10]
			if rname != "*":
				if flag == "83":
					if re.search(r"(\d+)S$",cigar):
						sf_n = int((re.search(r"(\d+)S$",cigar)).group(1))
						sf_seq = seq[-sf_n:]
						if sf_n >= polya_len and float(sf_seq.count("A"))/float(sf_n) >= a_ratio:
							read1_info = sam2gpd(line.strip())
							read1_flag = "83"
						else:
							read1_flag = "no_polya"
					else:
						read1_flag = "no_polya"
				elif flag == "99":
					if re.search(r"^(\d+)S",cigar):
						sf_n = int((re.search(r"^(\d+)S",cigar)).group(1))
						sf_seq = seq[:sf_n]
						if sf_n >= polya_len and float(sf_seq.count("T"))/float(sf_n) >= a_ratio:
							read1_info = sam2gpd(line.strip())
							read1_flag = "99"
						else:
							read1_flag = "no_polya"
					else:
						read1_flag = "no_polya"
				else:
					if read1_flag == "83":
						read2_info = sam2gpd(line.strip())
						tts_both = read1_info[5]
						info_first = read1_info[2:]
						tss_both = read2_info[4]
						info_first[2] = tss_both
						info_second = read2_info[-5:]
						info_both = info_first + info_second
						print >>gpd_file,qname + "\t" + "\t".join(info_both)
					elif read1_flag == "99":
						read2_info = sam2gpd(line.strip())
						tss_both = read1_info[4]
						info_first = read1_info[2:]
						tts_both = read2_info[5]
						info_first[3] = tts_both
						info_second = read2_info[-5:]
						info_both = info_first + info_second
						print >>gpd_file,qname + "\t" + "\t".join(info_both)
					else:
						pass
	sam_file.close()
	gpd_file.close()

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
15. exon end set (mate2)'''
	parser = argparse.ArgumentParser(description="Function: convert sam to gpd.",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-l','--length',type=int,default=5,help="PolyA tail length")
	parser.add_argument('-r','--ratio',type=float,default=0.8,help="Ratio of A in polyA tail")
	parser.add_argument('-i','--input',type=argparse.FileType('r'),required=True,help="Input: sam file")
	parser.add_argument('-o','--output',type=argparse.FileType('w'),required=True,help="Output: gpd file")
	args = parser.parse_args()
	return args

if __name__=="__main__":
	args = do_inputs()
	main(args)

