#!/usr/bin/env python
import sys,time,argparse

def main(args):
#	print >>sys.stdout, "Start analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S")
	separate(args.type_gpd,args.input,args.multi_exon,args.singlton_plus,args.singleton_minus)
#	print >>sys.stdout, "Finish analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S")

def separate(type,input_gpd,mlt_gpd,sgt_p_gpd,sgt_m_gpd):
	if type == "annotation":
		for line in input_gpd:
			gene_id,iso_id,chr,strand,exon_start,exon_end,cds_start,cds_end,exon_number,exon_start_set,exon_end_set = line.rstrip("\n").split("\t")
			gene_id = "refgene_" + gene_id
			iso_id = "refiso_" + iso_id
			if int(exon_number) == 1:
				if strand == "+":
					print >>sgt_p_gpd, gene_id + "\t" + iso_id + "\t" + "\t".join(line.rstrip("\n").split("\t")[2:])
				else:
					print >>sgt_m_gpd, gene_id + "\t" + iso_id + "\t" + "\t".join(line.rstrip("\n").split("\t")[2:])
			else:
				print >>mlt_gpd, gene_id + "\t" + iso_id + "\t" + "\t".join(line.rstrip("\n").split("\t")[2:])
	else:
		for line in input_gpd:
			gene_id,iso_id,chr,strand,exon_start,exon_end,cds_start,cds_end,exon_number,exon_start_set,exon_end_set = line.rstrip("\n").split("\t")
			if int(exon_number) == 1:
				if strand == "+":
					print >>sgt_p_gpd, gene_id + "\t" + iso_id + "\t" + "\t".join(line.rstrip("\n").split("\t")[2:])
				else:
					print >>sgt_m_gpd, gene_id + "\t" + iso_id + "\t" + "\t".join(line.rstrip("\n").split("\t")[2:])
			else:
				print >>mlt_gpd, gene_id + "\t" + iso_id + "\t" + "\t".join(line.rstrip("\n").split("\t")[2:])
	input_gpd.close()
	mlt_gpd.close()
	sgt_p_gpd.close()
	sgt_m_gpd.close()

def do_inputs():
	parser = argparse.ArgumentParser(description="Function: classify gpd file into three groups:(1)isoform with multiple exon;(2) isoform with single exon and plus strand; and (3) isoform with single exon and minus strand.",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-t','--type_gpd',type=str,choices=["annotation","alignment"],required=True,help="Type of input file, if 'annotation', add the prefix 'refgene_' and 'refiso_' to gene ID and isoform ID, respectively.")
	parser.add_argument('-i','--input',type=argparse.FileType('r'),required=True,help="Input: gpd file")
	parser.add_argument('-m','--multi_exon',type=argparse.FileType('w'),required=True,help="Output: gpd file containing multi-exon isoforms")
	parser.add_argument('--singlton_plus',type=argparse.FileType('w'),required=True,help="Output: gpd file containing plus-strand singleton isoforms")
	parser.add_argument('--singleton_minus',type=argparse.FileType('w'),required=True,help="Output: gpd file containing minus-strand singleton isoforms")
	args = parser.parse_args()
	return args

if __name__=="__main__":
	args = do_inputs()
	main(args)
