#!/usr/bin/env python
import sys,re,time,argparse

def main(args):
	print >>sys.stdout, "Start analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S")
	stat_list,known_sgt_na,known_sgt_pa,known_sgt_apa,known_mlt_na,known_mlt_pa,known_mlt_apa,novel_sgt_na,novel_sgt_pa,novel_sgt_apa,novel_mlt_na,novel_mlt_pa ,novel_mlt_apa = stat_results(args.input)
#	print stat_list
	print >>args.output, "Sample\APA type\tNA_no_polyA\tPA_with_one_polyA\tAPA_with_multiple_polyA\tTotal isoform number"
	print >>args.output, "%s\t%s\t%s\t%s\t%s" % ("Known singleton isoform",str(known_sgt_na),str(known_sgt_pa),str(known_sgt_apa),str(known_sgt_na+known_sgt_pa+known_sgt_apa))
	print >>args.output, "%s\t%s\t%s\t%s\t%s" % ("Known multi-exon isoform", str(known_mlt_na),str(known_mlt_pa),str(known_mlt_apa),str(known_mlt_na+known_mlt_pa+known_mlt_apa))
	print >>args.output, "%s\t%s\t%s\t%s\t%s" % ("Novel singleton isoform", str(novel_sgt_na),str(novel_sgt_pa),str(novel_sgt_apa),str(novel_sgt_na+novel_sgt_pa+novel_sgt_apa))
	print >>args.output, "%s\t%s\t%s\t%s\t%s" % ("Novel multi-exon isoform",str(novel_mlt_na),str(novel_mlt_pa),str(novel_mlt_apa),str(novel_mlt_na+novel_mlt_pa+novel_mlt_apa))
	print >>args.output, "%s\t%s\t%s\t%s\t%s" % ("Total isoform number",str(known_sgt_na+known_mlt_na+novel_sgt_na+novel_mlt_na),str(known_sgt_pa+known_mlt_pa+novel_sgt_pa+novel_mlt_pa),str(known_sgt_apa+known_mlt_apa+novel_sgt_apa+novel_mlt_apa),str(sum(stat_list)))
	print >>sys.stdout, "Finish analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S")

def stat_results(input_gpd):
	head = 1
	known_sgt_na = 0
	known_sgt_pa = 0
	known_sgt_apa = 0
	known_mlt_na = 0
	known_mlt_pa = 0
	known_mlt_apa = 0
	novel_sgt_na = 0
	novel_sgt_pa = 0
	novel_sgt_apa = 0
	novel_mlt_na = 0
	novel_mlt_pa = 0
	novel_mlt_apa = 0
	for line in input_gpd:
		if head:
			head -= 1
			continue
		gene_id,iso_id,chr,strand,tss,tts,cds_s,cds_e,exon_number,exon_start,exon_end,lr_pa_set,sr_pa_set,pa_set,pa_type = line.rstrip("\n").split("\t")
		if int(exon_number) == 1:
			if "novel_sgt_iso_" in iso_id:
				if pa_type == "NA":
					novel_sgt_na += 1
				elif pa_type == "PA":
					novel_sgt_pa += 1
				else:
					novel_sgt_apa += 1
			else:
				if pa_type == "NA":
					known_sgt_na += 1
				elif pa_type == "PA":
					known_sgt_pa += 1
				else:
					known_sgt_apa += 1
		else:
			if "novel_mlt_iso_" in iso_id:
				if pa_type == "NA":
					novel_mlt_na += 1
				elif pa_type == "PA":
					novel_mlt_pa += 1
				else:
					novel_mlt_apa += 1
			else:
				if pa_type == "NA":
					known_mlt_na += 1
				elif pa_type == "PA":
					known_mlt_pa += 1
				else:
					known_mlt_apa += 1
	stat_list = [known_sgt_na,known_sgt_pa,known_sgt_apa,known_mlt_na,known_mlt_pa,known_mlt_apa,novel_sgt_na,novel_sgt_pa,novel_sgt_apa,novel_mlt_na,novel_mlt_pa,novel_mlt_apa]
	return stat_list,known_sgt_na,known_sgt_pa,known_sgt_apa,known_mlt_na,known_mlt_pa,known_mlt_apa,novel_sgt_na,novel_sgt_pa,novel_sgt_apa,novel_mlt_na,novel_mlt_pa,novel_mlt_apa
	input_gpd.close()

def do_inputs():
	parser = argparse.ArgumentParser(description="Function: stat constructed isoforms and identified polyA sites",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i','--input',type=argparse.FileType('r'),required=True,help="Input: optimized isoforms with polyA site, gpd file")
	parser.add_argument('-o','--output',type=argparse.FileType('w'),required=True,help="Output: final output file")
	args = parser.parse_args()
	return args

if __name__=="__main__":
	args = do_inputs()
	main(args)	
