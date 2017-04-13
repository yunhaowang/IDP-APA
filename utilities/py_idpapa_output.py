#!/usr/bin/env python
import sys,re,time,argparse

def main(args):
#	print >>sys.stdout, "Start analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S")
	dic_iso_cds_start,dic_iso_cds_end = extract_cds(args.anno)
	output_modified_gpd(args.input,args.output,dic_iso_cds_start,dic_iso_cds_end)
#	print >>sys.stdout, "Finish analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S")

def extract_cds(input_anno):
	dic_iso_cds_start = {}
	dic_iso_cds_end = {}
	for line in input_anno:
		gene,iso,chr,strand,start,end,cds_start,cds_end = line.rstrip("\n").split("\t")[:8]
		chr_str_iso = chr + "&" + strand + "&" + iso
		if cds_start != "." and cds_end != ".":
			dic_iso_cds_start[chr_str_iso] = int(cds_start)
			dic_iso_cds_end[chr_str_iso] = int(cds_end)
		else:
			pass
	return dic_iso_cds_start,dic_iso_cds_end
	input_anno.close()

def output_modified_gpd(input_gpd,output_gpd,dic_iso_cds_start,dic_iso_cds_end):
	print >>output_gpd, "gene_id\tisoform_id\tchromosome\tstrand\tTSS\tTES\tCDS_start\tCDS_end\texon_number\texon_start_positions\texon_end_positions\tlong_reads_polyA_set\tshort_reads_polyA_set\tpolyA_site_set\tAPA_type"
	for line in input_gpd:
		gene_id,iso_id_set,chr,strand,tss,tts,lr_set,lr_count,exon_number,exon_start,exon_end,lr_pa_set,sr_pa_set,pa_set,pa_type = line.rstrip("\n").split("\t")
		boundary_set = set()
		for pa in lr_pa_set.split(","):
			boundary_set.add(int(pa.split("_")[0]))
		if sr_pa_set != "":
			for pa in sr_pa_set.split(","):
				boundary_set.add(int(pa.split("_")[0]))
		else:
			flag = 1
		if strand == "+":
			exon_start_c = exon_start.split(",")
			exon_start_c[0] = str(int(tss)+20)
			exon_end_c = exon_end.split(",")
			exon_end_c[-2] = str(max(boundary_set))
			tss = str(int(tss)+20)
			tts = str(max(boundary_set))
		else:
			exon_start_c = exon_start.split(",")
			exon_start_c[0] = str(min(boundary_set))
			exon_end_c = exon_end.split(",")
			exon_end_c[-2] = str(int(tts)-20)
			tss = str(min(boundary_set))
			tts = str(int(tts)-20)

		if "refiso_" in iso_id_set:
#			if int(lr_count) < lr_known_count:
#				continue
			cds_start_set = set()
			cds_end_set = set()
			for iso_id in iso_id_set.split(","):
				if (chr + "&" + strand + "&" + iso_id) in dic_iso_cds_start.keys():
					cds_start_set.add(dic_iso_cds_start[chr + "&" + strand + "&" + iso_id])
					cds_end_set.add(dic_iso_cds_end[chr + "&" + strand + "&" + iso_id])
				else:
					pass
			if cds_start_set != set():
				print >>output_gpd, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (re.sub("refgene_","",gene_id),re.sub("refiso_","",iso_id_set),chr,strand,tss,tts,str(min(cds_start_set)),str(max(cds_end_set)),exon_number,",".join(exon_start_c),",".join(exon_end_c),lr_pa_set,sr_pa_set,pa_set,pa_type)
			else:
				print >>output_gpd, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (re.sub("refgene_","",gene_id),re.sub("refiso_","",iso_id_set),chr,strand,tss,tts,".",".",exon_number,",".join(exon_start_c),",".join(exon_end_c),lr_pa_set,sr_pa_set,pa_set,pa_type)
		else:
#			if int(lr_count) < lr_novel_count:
#				continue
			print >>output_gpd, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (re.sub("refgene_","",gene_id),re.sub("refiso_","",iso_id_set),chr,strand,tss,tts,".",".",exon_number,",".join(exon_start_c),",".join(exon_end_c),lr_pa_set,sr_pa_set,pa_set,pa_type)
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
7. CDS start
8. CDS end
9. exon count
10. exon start set
11. exon end set
12. polyA sites identified by long reads
13. polyA sites identified by short reads
14. polyA sites merged by long reads and short reads
15. polyA type: (1) NA, no avaiable polyA site; (2) PA, one polyA site; (3) APA, multiple polyA sites'''
	parser = argparse.ArgumentParser(description="Function: output final file",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-a','--anno',type=argparse.FileType('r'),required=True,help="Annotation file, gpd file")
	parser.add_argument('-i','--input',type=argparse.FileType('r'),required=True,help="Input: optimized isoforms with polyA site, gpd file")
	parser.add_argument('-o','--output',type=argparse.FileType('w'),required=True,help="Output: final output, modified gpd file")
	args = parser.parse_args()
	return args

if __name__=="__main__":
	args = do_inputs()
	main(args)	
