import gerprate
import os
import re

#disopred data is in form disorderedResidues:residueCount:disorderedResiduesAfterPrematureStop:residueCountAfterPrematureStop:disorderedResiduePercentage:disorderedResiduesAfterPrematureStopPercentage:stopPosition
#can (preferably) pass in variant info in first argument instead of entire line
def getDisopredDataFromLine(line, disopredSequencesPath, variantType):
	newData = "."
	if variantType in line:
		prematureStopIndex = line.index(variantType)
		lineComponents = line[prematureStopIndex:].split(":")
		transcriptID = lineComponents[3]
		stopPosition = int(lineComponents[4].split("_")[2])
		try:
			disoFilePath = os.path.join(disopredSequencesPath, "%s.diso" % (transcriptID))
			disoFile = open(disoFilePath)

			#this file is in a terrible format, skip first 5 lines
			for _ in range(5):
				disoFile.readline()

			disorderedResidues = 0
			residueCount = 0

			disorderedResiduesAfterPrematureStop = 0
			residueCountAfterPrematureStop = 0

			for disoLine in disoFile:
				if disoLine.strip():
					disoLineComponents = [component for component in re.split(r'[\t ]', disoLine.strip()) if component]
					residueNumber = int(disoLineComponents[0])
					if disoLineComponents[2] == '*':
						disorderedResidues += 1
						if residueNumber >= stopPosition:
							disorderedResiduesAfterPrematureStop += 1

					residueCount += 1
					if residueNumber >= stopPosition:
						residueCountAfterPrematureStop += 1

			newData = "/".join([str(disorderedResidues), str(residueCount), str(disorderedResiduesAfterPrematureStop), str(residueCountAfterPrematureStop), "%.2f" % (100.0 * disorderedResidues / residueCount), "%.2f" % (100.0 * disorderedResiduesAfterPrematureStop / residueCountAfterPrematureStop), str(stopPosition)])
		except IOError:
			#print "Skipping transcript %s" % (transcriptID)
			pass
	return newData

def getGERPelements(elementFile):
	return [(int(eline.split('\t')[0]),int(eline.split('\t')[1]), float(eline.split('\t')[3])) for eline in elementFile]

def getCodingExonIntervals(annotationIntervalPath):
	codingExonIntervals = {}
	
	for line in open(annotationIntervalPath):
		lineComponents = line.split("\t")
		transcript = lineComponents[0].split("|")[1]
		chromosome = lineComponents[1].split("chr")[-1]
		intervalBeginComponents = lineComponents[6].split(",")
		intervalEndComponents = lineComponents[7].split(",")

		if not chromosome in codingExonIntervals:
			codingExonIntervals[chromosome] = {}
			
		#convert intervals from 0-based starting to 1-based starting
		intervals = [(int(intervalBeginComponents[intervalIndex])+1, int(intervalEndComponents[intervalIndex])) for intervalIndex in range(len(intervalBeginComponents))]
		codingExonIntervals[chromosome][transcript] = intervals

	return codingExonIntervals

def buildGerpRates(GERPratepath, GERPratecachepath, chromosome):
	ratefilepath = os.path.join(GERPratepath,'chr' +chromosome+'.maf.rates')
	cachepath = os.path.join(GERPratecachepath, chromosome+'.maf.rates_cached')
	try:
		#build cache file if it does not already exist using gerprate.so module
		gerprate.buildList(ratefilepath, cachepath)
	except:
		print ratefilepath + " could not be opened."
		print "Exiting program."
		sys.exit(1)

def getGerpScore(start, length):
	return sum([gerprate.floatInList(start+rateindex) for rateindex in range(length)]) / length

#binary search to find GERP element
def findGERPelementIndex(elements, start, end):
	low = 0
	high = len(elements) - 1
	while low<=high:
		mid = (low+high)/2
		if start > elements[mid][1]:
			low = mid+1
		elif end < elements[mid][0]:
			high = mid-1
		else:
			break
	if low <= high:
		return mid
	return -1
