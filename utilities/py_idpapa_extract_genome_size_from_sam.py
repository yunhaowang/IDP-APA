#!/usr/bin/env python
import sys,time,argparse

def main(args):
#	print >>sys.stdout, "Start analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S")
	extraction(args.input,args.output)
#	print >>sys.stdout, "Finish analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S")

def extraction(input_sam,output_txt):
	for line in input_sam:
		if line.startswith("@"):
			if line.startswith("@SQ"):
				sq,sn,ln = line.strip().split("\t")
				print >>output_txt, sn.split(":")[1] + "\t" + ln.split(":")[1]
			else:
				pass
		else:
			break
	input_sam.close()
	output_txt.close()

def do_inputs():
	output_format = '''
1. chromosome ID
2. length of chromosome'''
	parser = argparse.ArgumentParser(description="Function: extract genome size from sam.",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i','--input',type=argparse.FileType('r'),required=True,help="Input: sam file")
	parser.add_argument('-o','--output',type=argparse.FileType('w'),required=True,help="Output: genome size file")
	args = parser.parse_args()
	return args

if __name__=="__main__":
	args = do_inputs()
	main(args)
