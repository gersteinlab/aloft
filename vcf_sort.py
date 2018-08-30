#!/usr/bin/env python
#Usage of running this script by itself is <input_vcf> <output_vcf>
import os, sys, re

def sortVCFLines(lines):
	compiledRE = re.compile('([0-9]+)')
	def compareFunc(line):
		#find second tab index
		tabIndex = line.find("\t")
		tabIndex += line[tabIndex+1:].find("\t") + 1
		return [int(c) if c.isdigit() else c for c in compiledRE.split(line[:tabIndex])]
	lines.sort(key=compareFunc)

def sortVCF(inputPath, outputPath):
	regularLines = []
	outputFile = open(outputPath, "w")
	inputFile = open(inputPath)
	for line in inputFile:
		if line.startswith("#"):
			outputFile.write(line)
		else:
			regularLines.append(line.rstrip())
	inputFile.close()

	sortVCFLines(regularLines)
	
	for line in regularLines:
		outputFile.write(line + "\n")
	
	outputFile.close()

if __name__ == "__main__":
	if len(sys.argv) < 3:
		print("Usage: <input_path> <output_path>\n")
		print("Takes input_path and sorts it numerically to output_path. Input is a VCF file.")
		sys.exit(1)

	sortVCF(sys.argv[1], sys.argv[2])
