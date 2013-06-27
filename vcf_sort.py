#!/usr/bin/env python
#Usage of running this script by itself is <input_vcf> <output_vcf>
import os, sys, re

def tryint(s):
	try:
		return int(s)
	except:
		return s

def alphanum_key(s):
	return [ tryint(c) for c in re.split('([0-9]+)', s) ]

def sort_nicely(l):
	l.sort(key=alphanum_key)

if __name__ == "__main__":
	if len(sys.argv) < 3:
		print "Usage: <input_path> <output_path>\n"
		print "Takes input_path and sorts it numerically to output_path. Input is a VCF file."
		sys.exit(1)
	inputPath = sys.argv[1]
	outputPath = sys.argv[2]
	
	headerLines = []
	regularLines = []
	for line in open(inputPath):
		if line.startswith("#"):
			headerLines.append(line.rstrip("\n"))
		else:
			regularLines.append(line.rstrip("\n"))
	
	outputFile = open(outputPath, "w")
	
	sort_nicely(regularLines)
	
	if len(headerLines) > 0:
		outputFile.write("\n".join(headerLines) + "\n")
	
	if len(regularLines) > 0:
		outputFile.write("\n".join(regularLines) + "\n")
	
	outputFile.close()
