#!/usr/bin/env python
import sys,re,time,argparse

def main(args):
	print >>sys.stdout, "Start analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S")
	convert(args.input,args.output)
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

def convert(sam_file,gpd_file):
	for line in sam_file:
		if line[0] != "@":
			qname,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq = line.strip().split("\t")[:10]
			tag = "\t".join(line.strip().split("\t")[11:])
			if rname != "*" and re.search(r"XS:A:(\S)",tag):
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
				print >>gpd_file, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (qname,qname,rname,strand,str(int(pos)-1),str(exon_end),mapq,sf,str(len(exon_start_list)-1),",".join(exon_start_list),",".join(exon_end_list))
	sam_file.close()
	gpd_file.close()

def do_inputs():
	output_gpd_format = '''
1. read id
2. read id
3. chromosome id
4. strand
5. start site of alignment
6. end site of alignment
7. MAPQ 
8. Number of nucleotides that are softly-clipped by aligner; left_right
9. exon count
10. exon start set
11. exon end set'''
	parser = argparse.ArgumentParser(description="Function: convert sam to gpd.",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i','--input',type=argparse.FileType('r'),required=True,help="Input: sam file")
	parser.add_argument('-o','--output',type=argparse.FileType('w'),required=True,help="Output: gpd file")
	args = parser.parse_args()
	return args

if __name__=="__main__":
	args = do_inputs()
	main(args)
