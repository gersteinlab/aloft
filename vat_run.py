#!/usr/bin/env python
#You can just run VAT without running aloft if you want.

import os, sys
from subprocess import Popen, PIPE
import re
from vcf_sort import *
import gzip
from common import printError
import platform

def run_vat(arguments, forceVerbose=False):
	VAT_BIN_PATH = os.path.join('vat-bin', platform.system() + "_" + platform.machine())
	
	snpMapperPath = os.path.join(VAT_BIN_PATH, 'snpMapper')
	indelMapperPath = os.path.join(VAT_BIN_PATH, 'indelMapper')

	if not os.path.exists(snpMapperPath) or not os.path.exists(indelMapperPath):
		printError("VAT is not installed correctly - please see INSTALL")

	#example as input path: '/net/gerstein/sb238/ftw/finnish/Finns.nogeno.vcf.gz'
	# or a .vcf is fine too
	try:
		inputPath = arguments[1]
		vatOutputPath = arguments[2]
		annotationIntervalPath = arguments[3]
		annotationSequencePath = arguments[4]
	except:
		printError("Failed to parse arguments\nUsage is <input_vcf> <vat_output> <annotation_interval_input> <annotation_sequence_input> <verbosity_level>\nThis program will take care of sorting the input_vcf file numerically.\nFor verbosity_level you must pass in 0 (indicating no verbosity) or 1 (indicating verbosity)")
	verbose = False
	if forceVerbose or (5 < len(arguments) and int(arguments[5]) > 0):
		verbose = True

	if verbose: print('Parsing VCF file...')
	
	try:
		if inputPath.endswith(".gz"):
			inputFile = gzip.open(inputPath, 'rb')
		else:
			inputFile = open(inputPath, 'rb')
	except:
		printError("Failed to open %s" % (inputPath))
	
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
			if line.startswith("#CHR"):
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
	
	if verbose: print("Running snpMapper...")
	try:
		snpMapperPipe = Popen([snpMapperPath, annotationIntervalPath, annotationSequencePath], stdout=PIPE, stdin=snpInputFile)
	except:
		printError("Failed to open snpMapper")
	
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
	
	if verbose: print("Running indelMapper...")
	try:
		indelMapperPipe = Popen([indelMapperPath, annotationIntervalPath, annotationSequencePath], stdout=PIPE, stdin=indelInputFile)
	except:
		printError("Failed to open indelMapper")
	
	for lineBytes in indelMapperPipe.stdout:
		line = lineBytes.decode()
		if not line.startswith("#"):
			sortedLines.append(line.rstrip("\n"))
			numIndel += 1
	
	indelInputFile.close()
	
	os.remove(TEMP_INDEL_PATH)
	
	if verbose: print("Writing out VAT file...")
	#Sort and remove duplicate entries
	sort_nicely(sortedLines)
	lastLine = None
	for line in sortedLines:
		if line != lastLine:
			vcfOutputFile.write(line + "\n")
		lastLine = line
	
	vcfOutputFile.close()
	
	if verbose: print("Finishing VAT. There were %d snp lines and %d indel lines in the output" % (numSnp, numIndel))

if __name__ == "__main__":
	if len(sys.argv) < 6:
		print("Usage: %s <input_vcf> <output_vat> <interval_file> <sequences_file> <verbosity_level>" % sys.argv[0])
		print("%s will automatically take care of sorting the input_vcf file numerically" % sys.argv[0])
		print("For verbosity_level you must pass in 0 (indicating no verbosity) or 1 (indicating verbosity)")
		sys.exit(1)
	run_vat(sys.argv)
