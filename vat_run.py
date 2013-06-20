#You can just run VAT. Usage is <input_vcf> <vat_output> <annotation_interval_input> <annotation_sequence_input> This script will take care of sorting the input_vcf file numerically.

import os, sys
from subprocess import Popen, PIPE
import re

def tryint(s):
	try:
		return int(s)
	except:
		return s

def alphanum_key(s):
	return [ tryint(c) for c in re.split('([0-9]+)', s) ]

def sort_nicely(l):
	l.sort(key=alphanum_key)

def run_vat(arguments):
	#inputPath = '/net/gerstein/sb238/ftw/finnish/Finns.nogeno.vcf.gz'
	# or a .vcf
	print 'Parsing VCF file...'
	inputPath = arguments[1]
	vatOutputPath = arguments[2]
	annotationIntervalPath = arguments[3]
	annotationSequencePath = arguments[4]
	
	try:
		open(inputPath)
	except:
		print "Failed to open " + inputPath
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
	for line in inputFile:
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
	
	print "Running snpMapper..."
	#snpMapperPipe = Popen(['/net/gerstein/sb238/vat/bin/snpMapper', '/net/gerstein/sb238/MNP/gencode12.CCDSinclusive.interval', '/net/gerstein/sb238/MNP/gencode12.CCDSinclusive.fa'], stdout=PIPE, stdin=snpInputFile)
	try:
		snpMapperPipe = Popen(['/net/gerstein/sb238/vat/bin/snpMapper', annotationIntervalPath, annotationSequencePath], stdout=PIPE, stdin=snpInputFile)
	except:
		print "ERROR: Failed to open snpMapper. Is snpMapper in your PATH?"
		sys.exit(1)
	
	numSnp = 0
	numIndel = 0
	
	sortedLines = []
	
	vcfOutputFile = open(vatOutputPath, "w")
	for line in snpMapperPipe.stdout:
		if line.startswith("#"):
			vcfOutputFile.write(line)
		else:
			numSnp += 1
			sortedLines.append(line.rstrip("\n"))
	
	snpInputFile.close()
	
	os.remove(TEMP_SNP_PATH)
	
	print "Running indelMapper..."
	#indelMapperPipe = Popen(['/net/gerstein/sb238/vat/bin/indelMapper', '/net/gerstein/sb238/MNP/gencode12.CCDSinclusive.interval', '/net/gerstein/sb238/MNP/gencode12.CCDSinclusive.fa'], stdout=PIPE, stdin=indelInputFile)
	try:
		indelMapperPipe = Popen(['/net/gerstein/sb238/vat/bin/indelMapper', annotationIntervalPath, annotationSequencePath], stdout=PIPE, stdin=indelInputFile)
	except:
		print "ERROR: Failed to open indelMapper. Is indelMapper in your PATH?"
		sys.exit(1)
	
	for line in indelMapperPipe.stdout:
		if not line.startswith("#"):
			sortedLines.append(line.rstrip("\n"))
			numIndel += 1
	
	indelInputFile.close()
	
	os.remove(TEMP_INDEL_PATH)
	
	print "Writing out VAT file..."
	#Sort and remove duplicate entries
	sort_nicely(sortedLines)
	lastLine = None
	for line in sortedLines:
		if line != lastLine:
			vcfOutputFile.write(line + "\n")
		lastLine = line
	
	vcfOutputFile.close()
	
	print "There were %d snp lines and %d indel lines" % (numSnp, numIndel)

if __name__ == "__main__":
	run_vat(sys.argv)
