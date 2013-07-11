#!/usr/bin/env python
#You can just run VAT. Usage is <input_vcf> <vat_output> <annotation_interval_input> <annotation_sequence_input> This script will take care of sorting the input_vcf file numerically.

import os, sys
from subprocess import Popen, PIPE
import re
from vcf_sort import *

VAT_BIN_PATH = "vat/vat-2.0.1-install/bin/"

def run_vat(arguments):
	#example as input path: '/net/gerstein/sb238/ftw/finnish/Finns.nogeno.vcf.gz'
	# or a .vcf is fine too
	print('Parsing VCF file...')
	inputPath = arguments[1]
	vatOutputPath = arguments[2]
	annotationIntervalPath = arguments[3]
	annotationSequencePath = arguments[4]
	
	try:
		open(inputPath)
	except:
		print("Failed to open " + inputPath)
		sys.exit(1)
	
	if inputPath.endswith(".gz"):
		zcatPipe = Popen(['zcat', inputPath], stdout=PIPE)
		inputFile = zcatPipe.stdout
	else:
		inputFile = open(inputPath)
	
	TEMP_SNP_PATH = os.path.join(os.path.split(vatOutputPath)[0],"snpinput_temp")
	TEMP_INDEL_PATH = os.path.join(os.path.split(vatOutputPath)[0],"indelinput_temp")
	
	snpInputFile = open(TEMP_SNP_PATH, "w")
	indelInputFile = open(TEMP_INDEL_PATH, "w")
	
	foundHeader = False
	foundID = True
	numberOfMissingComponents = 0
	normalHeaderComponents = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
	for lineBytes in inputFile:
		line = lineBytes.decode()
		lineComponents = line.rstrip("\n").rstrip("\t").split("\t")

		if line.startswith("#"):
			if line.startswith("#CHROM"):
				if lineComponents[2] != 'ID':
					foundID = False
				if len(lineComponents) < 8:
					numberOfMissingComponents = 8 - len(lineComponents)
					if not foundID:
						numberOfMissingComponents -= 1

				if not foundID:
					lineComponents = lineComponents[0:2] + ['ID'] + lineComponents[2:]

				if numberOfMissingComponents > 0:
					lineComponents += normalHeaderComponents[-numberOfMissingComponents:]

				snpInputFile.write('\t'.join(lineComponents) + "\n")
				indelInputFile.write('\t'.join(lineComponents) + "\n")

				foundHeader = True
			else:
				snpInputFile.write(line)
				indelInputFile.write(line)
		elif not foundHeader:
			snpInputFile.write("#" + "\t".join(normalHeaderComponents) + "\n")
			indelInputFile.write("#" + "\t".join(normalHeaderComponents) + "\n")
			foundHeader = True
		if not line.startswith("#"):
			if not foundID:
				lineComponents = lineComponents[0:2] + ['NA'] + lineComponents[2:]
			if numberOfMissingComponents > 0:
				lineComponents += ['NA'] * numberOfMissingComponents
			refComponents = lineComponents[3].split(",")
			altComponents = lineComponents[4].split(",")

			foundSnp = False
			foundIndel = False

			for index in range(len(refComponents)):
				refComponent = refComponents[index]
				altComponent = altComponents[index]
				if len(refComponent) == 1 and len(altComponent) == 1 and not foundSnp:
					snpInputFile.write("\t".join(lineComponents) + "\n")
					foundSnp = True
					if foundIndel:
						break

				if (len(refComponent) > 1 or len(altComponent) > 1) and not foundIndel:
					indelInputFile.write("\t".join(lineComponents) + "\n")
					foundIndel = True
					if foundSnp:
						break
	
	snpInputFile.close()
	indelInputFile.close()
	
	snpInputFile = open(TEMP_SNP_PATH, "r")
	indelInputFile = open(TEMP_INDEL_PATH, "r")
	
	print("Running snpMapper...")
	try:
		snpMapperPipe = Popen([os.path.join(VAT_BIN_PATH, 'snpMapper'), annotationIntervalPath, annotationSequencePath], stdout=PIPE, stdin=snpInputFile)
	except:
		print("ERROR: Failed to open snpMapper")
		sys.exit(1)
	
	numSnp = 0
	numIndel = 0
	
	sortedLines = []
	
	vcfOutputFile = open(vatOutputPath, "w")
	for lineBytes in snpMapperPipe.stdout:
		line = lineBytes.decode()
		if line.startswith("#"):
			vcfOutputFile.write(line)
		else:
			numSnp += 1
			sortedLines.append(line.rstrip("\n"))
	
	snpInputFile.close()
	
	os.remove(TEMP_SNP_PATH)
	
	print("Running indelMapper...")
	try:
		indelMapperPipe = Popen([os.path.join(VAT_BIN_PATH, 'indelMapper'), annotationIntervalPath, annotationSequencePath], stdout=PIPE, stdin=indelInputFile)
	except:
		print("ERROR: Failed to open indelMapper.")
		sys.exit(1)
	
	for lineBytes in indelMapperPipe.stdout:
		line = lineBytes.decode()
		if not line.startswith("#"):
			sortedLines.append(line.rstrip("\n"))
			numIndel += 1
	
	indelInputFile.close()
	
	os.remove(TEMP_INDEL_PATH)
	
	print("Writing out VAT file...")
	#Sort and remove duplicate entries
	sort_nicely(sortedLines)
	lastLine = None
	for line in sortedLines:
		if line != lastLine:
			vcfOutputFile.write(line + "\n")
		lastLine = line
	
	vcfOutputFile.close()
	
	print("There were %d snp lines and %d indel lines" % (numSnp, numIndel))

if __name__ == "__main__":
	run_vat(sys.argv)
