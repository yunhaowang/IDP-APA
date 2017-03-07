#!/usr/bin/env python
import sys,time,argparse

def main(args):
	print >>sys.stdout, "Start analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S")
	construction(args.anno,args.input,args.output,args.lr_known,args.lr_novel)
	print >>sys.stdout, "Finish analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S")

#=== extract gene region from known gene annotation library ===
def extract_gene_region_from_annotation(anno_gpd):
	dic_chr_gene_region = {}
	for line in anno_gpd:
		gene_id,isoform_id,chrom,strand,tss,tts,cds_start,cds_end,exon_number,exon_start,exon_end = line.strip().split("\t")
		chr_strand = chrom
		if chr_strand not in dic_chr_gene_region.keys():
			dic_chr_gene_region[chr_strand] = set()
			dic_chr_gene_region[chr_strand].add(tss+"_"+tts)
		else:
			dic_chr_gene_region[chr_strand].add(tss+"_"+tts)
	anno_gpd.close()
	return dic_chr_gene_region

def construction(anno_gpd,concat_gpd,constructed_gpd,lr_count_known,lr_count_novel): # construct isoform
	dic_chr_gene_region = extract_gene_region_from_annotation(anno_gpd)
	sgt_novel_index = 0
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
			flag = "non-overlap"
			chr_strand = chrom
			if chr_strand in dic_chr_gene_region.keys():
				for tss_tts in dic_chr_gene_region[chr_strand]:
					if int(tss) < int(tss_tts.split("_")[1]) and int(tts) > int(tss_tts.split("_")[0]): # overlap with known gene region
						flag = "overlap"
						break
					else:
						flag = "non-overlap"
			else:
				flag = "non-overlap"
			if flag == "non-overlap": # output novel isoform
				sgt_novel_index += 1
				sgt_novel_loci = "novel_loci_sgt_" + str(sgt_novel_index)
				sgt_novel_name = "novel_sgt_iso_" + str(sgt_novel_index)
				print >>constructed_gpd, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (sgt_novel_loci,sgt_novel_name,chrom,strand,tss,tts,gene_set,gene_number,exon_number,exon_start,exon_end)
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
	parser = argparse.ArgumentParser(description="Function: construct singleton isoform",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-a','--anno',type=argparse.FileType('r'),required=True,help="Input: annotation libray (gpd file)")
	parser.add_argument('--lr_known',type=int,default=1,help="number of long reads for supporting a known isoform")
	parser.add_argument('--lr_novel',type=int,default=1,help="number of long reads for supporting a novel isoform")
	parser.add_argument('-i','--input',type=argparse.FileType('r'),required=True,help="Input: concatenated singleton isoform (gpd file)")
	parser.add_argument('-o','--output',type=argparse.FileType('w'),required=True,help="Output: constructed isoform (gpd file)")
	args = parser.parse_args()
	return args

if __name__=="__main__":
	args = do_inputs()
	main(args)
