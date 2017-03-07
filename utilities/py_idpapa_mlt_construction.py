#!/usr/bin/env python
import sys,time,argparse

def main(args):
	print >>sys.stdout, "Start analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S")
	dic_junction_set = extract_junction_set_from_annotation(args.anno)
	if args.short_reads:
		dic_junction_sr = extract_junction_from_short_reads(args.short_reads)
		dic_junction_anno,dic_chr_strand_iso = extract_junction_from_annotation(args.anno)
		dic_junction = {}
		chr_set = set(dic_junction_sr.keys()+dic_junction_anno.keys())
		for chr in chr_set:
			dic_junction[chr] = {}
			if chr in dic_junction_sr.keys() and chr in dic_junction_anno.keys():
				for type in dic_junction_sr[chr].keys():
					dic_junction[chr][type] = set(list(dic_junction_sr[chr][type]) + list(dic_junction_anno[chr][type]))
			elif chr in dic_junction_sr.keys() and chr not in dic_junction_anno.keys():
				for type in dic_junction_sr[chr].keys():
					dic_junction[chr][type] = dic_junction_sr[chr][type]
			elif chr not in dic_junction_sr.keys() and chr in dic_junction_anno.keys():
				for type in dic_junction_anno[chr].keys():
					dic_junction[chr][type] = dic_junction_anno[chr][type]
			else:
				pass
	else:
		dic_junction,dic_chr_strand_iso = extract_junction_from_annotation(args.anno)
#	print >>sys.stdout, dic_junction_set
#	print >>sys.stdout, dic_junction
	construction(dic_junction,dic_chr_strand_iso,dic_junction_set,args.input,args.output,args.lr_known,args.lr_novel)
	print >>sys.stdout, "Finish analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S")

#=== extract junction site from known gene annotation library ===
def extract_junction_from_annotation(anno_gpd):
	anno_file = open(anno_gpd,"r")
	dic_chr_type_junction = {}
	dic_chr_strand_iso = {}
	flag = 0
	for line in anno_file:
		gene_id,isoform_id,chrom,strand,tss,tts,cds_start,cds_end,exon_number,exon_start,exon_end = line.rstrip("\n").split("\t")
		if int(exon_number) > 1:
			chr_strand = chrom + "_" + strand
			if chr_strand not in dic_chr_strand_iso.keys():
				dic_chr_strand_iso[chr_strand] = {}
#				continue
			else:
				flag = 1
			if strand == "+":
				if chrom not in dic_chr_type_junction.keys():
					dic_chr_type_junction[chrom] = {}
					for type in ["sorted_plus5","sorted_plus3","sorted_minus5","sorted_minus3"]:
						dic_chr_type_junction[chrom][type] = set()
#					continue
				for p5 in exon_end.split(",")[:-2]:
					dic_chr_type_junction[chrom]["sorted_plus5"].add(int(p5))
				for p3 in exon_start.split(",")[1:-1]:
					dic_chr_type_junction[chrom]["sorted_plus3"].add(int(p3))
			else:
				if chrom not in dic_chr_type_junction.keys():
					dic_chr_type_junction[chrom] = {}
					for type in ["sorted_plus5","sorted_plus3","sorted_minus5","sorted_minus3"]:
						dic_chr_type_junction[chrom][type] = set()
#					continue
				for p5 in exon_end.split(",")[:-2]:
					dic_chr_type_junction[chrom]["sorted_minus5"].add(int(p5))
				for p3 in exon_start.split(",")[1:-1]:
					dic_chr_type_junction[chrom]["sorted_minus3"].add(int(p3))
			dic_chr_strand_iso[chr_strand][isoform_id] = {}
			dic_chr_strand_iso[chr_strand][isoform_id]["gene_id"] = gene_id
			dic_chr_strand_iso[chr_strand][isoform_id]["tss"] = int(tss)
			dic_chr_strand_iso[chr_strand][isoform_id]["tts"] = int(tts)
	return dic_chr_type_junction,dic_chr_strand_iso
	anno_file.close()

#=== extract junction site from short reads ===
def extract_junction_from_short_reads(sr_gpd):
	dic_chr_type_junction = {}
	for line in sr_gpd:
		gene_id,isoform_id,chrom,strand,tss,tts,cds_start,cds_end,exon_number,exon_start,exon_end = line.rstrip("\n").split("\t")
		if strand == "+" and int(exon_number) != 1:
			if chrom not in dic_chr_type_junction.keys():
				dic_chr_type_junction[chrom] = {}
				for type in ["sorted_plus5","sorted_plus3","sorted_minus5","sorted_minus3"]:
					dic_chr_type_junction[chrom][type] = set()
			for p5 in exon_end.split(",")[:-2]:
				dic_chr_type_junction[chrom]["sorted_plus5"].add(int(p5))
			for p3 in exon_start.split(",")[1:-1]:
				dic_chr_type_junction[chrom]["sorted_plus3"].add(int(p3))
		elif strand == "-" and int(exon_number) != 1:
			if chrom not in dic_chr_type_junction.keys():
				dic_chr_type_junction[chrom] = {}
				for type in ["sorted_plus5","sorted_plus3","sorted_minus5","sorted_minus3"]:
					dic_chr_type_junction[chrom][type] = set()
			for p5 in exon_end.split(",")[:-2]:
				dic_chr_type_junction[chrom]["sorted_minus5"].add(int(p5))
			for p3 in exon_start.split(",")[1:-1]:
				dic_chr_type_junction[chrom]["sorted_minus3"].add(int(p3))
		else:
			pass
	return dic_chr_type_junction
	sr_gpd.close()

#=== extract junction set from known gene annotation library ===
def extract_junction_set_from_annotation(anno_gpd):
	anno_file = open(anno_gpd,"r")
	dic_junction_set = {}
	for line in anno_file:
		gene_id,isoform_id,chrom,strand,tss,tts,cds_start,cds_end,exon_number,exon_start,exon_end = line.rstrip("\n").split("\t")
		if int(exon_number) > 1:
			chr_strand = chrom + "_" + strand
			if chr_strand not in dic_junction_set.keys():
				dic_junction_set[chr_strand] = set()
				dic_junction_set[chr_strand].add("-".join(exon_start.split(",")[1:-1])+"_"+"-".join(exon_end.split(",")[:-2]))
			else:
				dic_junction_set[chr_strand].add("-".join(exon_start.split(",")[1:-1])+"_"+"-".join(exon_end.split(",")[:-2]))
	return dic_junction_set
	anno_file.close()

#=== test if all splice sites of spliced aligned long read are in annotation library ===
def test_splice_sites(lr_mult_gpd,dic_junction):
	lr_id,lr_name,chrom,strand,tss,tts,mapq,sf,exon_number,exon_start,exon_end = lr_mult_gpd.split("\t")
	flag_p5 = "no"
	flag_p3 = "no"
	if chrom in dic_junction.keys() and strand == "+":
		for p5 in exon_end.split(",")[:-2]:
			if int(p5) in dic_junction[chrom]["sorted_plus5"]:
				flag_p5 = "yes"
			else:
				flag_p5 = "no"
				break
		for p3 in exon_start.split(",")[1:-1]:
			if int(p3) in dic_junction[chrom]["sorted_plus3"]:
				flag_p3 = "yes"
			else:
				flag_p3 = "no"
				break
	elif chrom in dic_junction.keys() and strand == "-":
		for p5 in exon_end.split(",")[:-2]:
			if int(p5) in dic_junction[chrom]["sorted_minus5"]:
				flag_p5 = "yes"
			else:
				flag_p5 = "no"
				break
		for p3 in exon_start.split(",")[1:-1]:
			if int(p3) in dic_junction[chrom]["sorted_minus3"]:
				flag_p3 = "yes"
			else:
				flag_p3 = "no"
				break
	else:
		flag_p5 = "no"
		flag_p3 = "no"
	return flag_p5,flag_p3	

def assign_to_gene(mlt_gpd,dic_chr_strand_iso):
	lr_id,lr_name,chrom,strand,tss,tts,mapq,sf,exon_number,exon_start,exon_end = mlt_gpd.split("\t")
	chr_strand = chrom + "_" + strand
	gene_id = ""
	dic_read_ref_identity = {}
	if chr_strand in dic_chr_strand_iso.keys():
		for iso in dic_chr_strand_iso[chr_strand].keys():
			if int(tss) < dic_chr_strand_iso[chr_strand][iso]["tts"] and int(tts) > dic_chr_strand_iso[chr_strand][iso]["tss"]:
				read_ref_identity = float(min(int(tts),dic_chr_strand_iso[chr_strand][iso]["tts"])-max(int(tss),dic_chr_strand_iso[chr_strand][iso]["tss"]))/float(max(int(tts),dic_chr_strand_iso[chr_strand][iso]["tts"])-min(int(tss),dic_chr_strand_iso[chr_strand][iso]["tss"]))
				dic_read_ref_identity[read_ref_identity] = dic_chr_strand_iso[chr_strand][iso]["gene_id"]
			else:
				pass
		if dic_read_ref_identity != {}:
			gene_id = dic_read_ref_identity[max(dic_read_ref_identity.keys())]
		else:
			gene_id = ""
	else:
		gene_id = ""
	return gene_id

def construction(dic_junction,dic_chr_strand_iso,dic_junction_set,concat_gpd,constructed_gpd,lr_count_known,lr_count_novel): # construct isoform
	mlt_novel_iso_index = 0
	mlt_novel_loci_index = 0
	for line in concat_gpd:
		gene_set,iso_set,chrom,strand,tss,tts,gene_number,iso_number,exon_number,exon_start,exon_end = line.rstrip("\n").split("\t")
		if "refgene_" in gene_set: # known isoform
			refgene_set = set()
			refiso_set = set()
			lr_set = set()
			for gene_id in gene_set.split(","):
				if gene_id.startswith("refgene_"):
					refgene_set.add(gene_id)
				else:
					lr_set.add(gene_id)
			if len(lr_set) >= lr_count_known: # number of supported long reads
				for iso_id in iso_set.split(","):
					if iso_id.startswith("refiso_"):
						refiso_set.add(iso_id)
				refgene_set = list(refgene_set)
				refiso_set = list(refiso_set)
				refgene_set.sort()
				refiso_set.sort()
				print >>constructed_gpd, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (",".join(refgene_set),",".join(refiso_set),chrom,strand,tss,tts,",".join(lr_set),str(len(lr_set)),exon_number,exon_start,exon_end)
		else: # novel isoform
			if int(gene_number) < lr_count_novel: # number of supported long reads 
				continue
			lr_gpd = line.rstrip("\n")
			flag_p5,flag_p3 = test_splice_sites(lr_gpd,dic_junction)
			if flag_p5 == "yes" and flag_p3 == "yes": # if all splice sites are annotated by annotated library
				flag = "non-overlap"
				p5 = "-".join(exon_start.split(",")[1:-1])
				p3 = "-".join(exon_end.split(",")[:-2])
				chr_strand = chrom+"_"+strand
				if chr_strand in dic_junction_set.keys():
					for p5_p3 in dic_junction_set[chr_strand]:
						if p5 in p5_p3.split("_")[0] and p3 in p5_p3.split("_")[1]: # if the junction combination set is a subset of any known isoform
							flag = "overlap"
							break
						else:
							flag = "non-overlap"
				else:
					flag = "non-overlap"
				if flag == "non-overlap":
					mlt_novel_iso_index += 1
					mlt_novel_iso_name = "novel_mlt_iso_" + str(mlt_novel_iso_index)
					gene_id = assign_to_gene(lr_gpd,dic_chr_strand_iso)
					if gene_id == "": # novel loci
						mlt_novel_loci_index += 1
						mlt_novel_loci_name = "novel_mlt_loci_" + str(mlt_novel_loci_index)
						print >>constructed_gpd, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (mlt_novel_loci_name,mlt_novel_iso_name,chrom,strand,tss,tts,gene_set,gene_number,exon_number,exon_start,exon_end)
					else: # known loci and novel isoform
						print >>constructed_gpd, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (gene_id,mlt_novel_iso_name,chrom,strand,tss,tts,gene_set,gene_number,exon_number,exon_start,exon_end)
				else:
					pass
	
	concat_gpd.close()
	constructed_gpd.close()

def do_inputs():
	output_gpd_format = '''
1. gene id
2. isoform id
3. chromosome id
4. strand
5. TSS
6. TTS
7. long read set
8. long read count
9. exon count
10. exon start set
11. exon end set'''
	parser = argparse.ArgumentParser(description="Function: construct multi-exon isoform",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-a','--anno',type=str,required=True,help="Input: annotation libray (gpd file)")
	parser.add_argument('-s','--short_reads',type=argparse.FileType('r'),help="Optional input: short reads for confirming splice sites (gpd file)")
	parser.add_argument('--lr_known',type=int,default=1,help="number of long reads for supporting a known isoform")
	parser.add_argument('--lr_novel',type=int,default=1,help="number of long reads for supporting a novel isoform")
	parser.add_argument('-i','--input',type=argparse.FileType('r'),required=True,help="Input: concatenated singleton isoform (gpd file)")
	parser.add_argument('-o','--output',type=argparse.FileType('w'),required=True,help="Output: constructed isoform (gpd file)")
	args = parser.parse_args()
	return args

if __name__=="__main__":
	args = do_inputs()
	main(args)
