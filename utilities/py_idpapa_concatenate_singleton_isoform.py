#!/usr/bin/env python
import sys,time,argparse

def main(args):
	print >>sys.stdout, "Start analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S")
	concatenate(args.input,args.output)
	print >>sys.stdout, "Finish analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S")

def concatenate(anno_gpd,concatenate_anno_gpd):
	dic_sgt = {}
	chr_set = set()
	for line in anno_gpd:
		gene_id,isoform_id,chr,strand,tss,tts,cds_start,cds_end,exon_number,exon_start,exon_end = line.rstrip("\n").split("\t")
		if chr not in chr_set:
			if chr_set != set():
				print >>concatenate_anno_gpd, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s,\t%s," % (",".join(dic_sgt["gene_id"]),",".join(dic_sgt["isoform_id"]),dic_sgt["chr"],dic_sgt["strand"],dic_sgt["tss"],dic_sgt["tts"],str(len(dic_sgt["gene_id"])),str(len(dic_sgt["isoform_id"])),"1",dic_sgt["tss"],dic_sgt["tts"])
			chr_set.add(chr)
			dic_sgt["gene_id"] = set()
			dic_sgt["isoform_id"] = set()
			dic_sgt["chr"] = chr
			dic_sgt["strand"] = strand
			dic_sgt["tss"] = tss
			dic_sgt["tts"] = tts
			dic_sgt["gene_id"].add(gene_id)
			dic_sgt["isoform_id"].add(isoform_id)
		else:
			if int(tss) < int(dic_sgt["tts"]):
				dic_sgt["tts"] = str(max(int(tts),int(dic_sgt["tts"])))
				dic_sgt["gene_id"].add(gene_id)
				dic_sgt["isoform_id"].add(isoform_id)
			else:
				print >>concatenate_anno_gpd, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s,\t%s," % (",".join(dic_sgt["gene_id"]),",".join(dic_sgt["isoform_id"]),dic_sgt["chr"],dic_sgt["strand"],dic_sgt["tss"],dic_sgt["tts"],str(len(dic_sgt["gene_id"])),str(len(dic_sgt["isoform_id"])),"1",dic_sgt["tss"],dic_sgt["tts"])
				dic_sgt["gene_id"] = set()
				dic_sgt["isoform_id"] = set()
				dic_sgt["chr"] = chr
				dic_sgt["strand"] = strand
				dic_sgt["tss"] = tss
				dic_sgt["tts"] = tts
				dic_sgt["gene_id"].add(gene_id)
				dic_sgt["isoform_id"].add(isoform_id)

	print >>concatenate_anno_gpd, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s,\t%s," % (",".join(dic_sgt["gene_id"]),",".join(dic_sgt["isoform_id"]),dic_sgt["chr"],dic_sgt["strand"],dic_sgt["tss"],dic_sgt["tts"],str(len(dic_sgt["gene_id"])),str(len(dic_sgt["isoform_id"])),"1",dic_sgt["tss"],dic_sgt["tts"])
	anno_gpd.close()
	concatenate_anno_gpd.close()

def do_inputs():
	output_gpd_format = '''
1. gene id
2. isoform id
3. chromosome id
4. strand
5. TSS
6. TTS
7. gene count
8. isoform count
9. exon count
10. exon start set
11. exon end set'''
	parser = argparse.ArgumentParser(description="Function: concatenate singleton isoforms if they have the same strand and overlap. Note: (1) gpd file only includes one strand (plus or minus) information (2) gpd file must be sorted by chromosome,start,end using the command 'sort -k3,3 -k5,5n -k6,6n'",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i','--input',type=argparse.FileType('r'),required=True,help="Input: gpd file")
	parser.add_argument('-o','--output',type=argparse.FileType('w'),required=True,help="Output: concatenated gpd file")
	args = parser.parse_args()
	return args

if __name__=="__main__":
	args = do_inputs()
	main(args)
