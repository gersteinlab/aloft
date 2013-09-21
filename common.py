import os
import re
import struct
import subprocess
import sys

def printError(error, exit=True):
	sys.stderr.write("Error: %s\n" % error)
	if exit:
		sys.stderr.write("Exiting..\n")
		sys.exit(1)

def getTruncatedExons(exons, start, direction):
	truncatedExons = None
	stopExonIndex = 0
	for block in exons:
		if block[0] <= start and block[1] >= start:
			if direction == '+':
				truncatedExons = exons[stopExonIndex:]
			elif direction == '-':
				truncatedExons = exons[0:stopExonIndex+1]
			break
		stopExonIndex += 1
	return truncatedExons

def mergeElements(elements):
	mergedElements = list(elements)
	while True:
		mergedElements, didMerge = _mergeElements(mergedElements)
		if not didMerge:
			break
	return mergedElements

def _mergeElements(elements):
	didMerge = False
	mergedElements = []
	for elementIndex, element in enumerate(elements):
		if elementIndex+1 < len(elements):
			nextElement = elements[elementIndex+1]
			#skip next element since it's same as current one
			if nextElement[0] == element[0] and nextElement[1] == element[1]:
				mergedElements.append(element)
				didMerge = True
				if elementIndex+2 < len(elements):
					mergedElements += elements[elementIndex+2:]
				break
			#merge next element with current one since they intersect
			elif nextElement[0] <= element[1]:
				didMerge = True
				if nextElement[1] < element[1]:
					mergedElements.append(element)
				else:
					#arbitrarily choose one of the elements to merge the last field from
					mergedElements.append([element[0], nextElement[1]] + element[2:])
				if elementIndex+2 < len(elements):
					mergedElements += elements[elementIndex+2:]
				break
			#no intersection
			else:
				mergedElements.append(element)
		#last element
		else:
			mergedElements.append(element)
	return mergedElements, didMerge

def getRejectionElementIntersectionData(exons, truncatedExons, GERPelements, GERPelementIndex, direction):
	rejectedElements = []

	elementIndex = GERPelementIndex

	#Make a list of elements that may be relevant to the truncated exons
	relevantElements = []
	while elementIndex >= 0 and elementIndex < len(GERPelements) and ((direction == '-' and GERPelements[elementIndex][1] >= truncatedExons[0][0]) or (direction == '+' and GERPelements[elementIndex][0] <= truncatedExons[-1][1])):
		relevantElements.append(GERPelements[elementIndex])
		if direction == '+':
			elementIndex += 1
		elif direction == '-':
			elementIndex -= 1

	#find all truncated exons and elements that intersect
	for truncatedExon in truncatedExons:
		for relevantElement in relevantElements:
			#check if they overlap in any way
			if truncatedExon[0] <= relevantElement[1] and truncatedExon[1] >= relevantElement[0]:
				distanceCovered = 0
				exonLength = truncatedExon[1] - truncatedExon[0] + 1
	
				#see if element completly contains exon
				if relevantElement[0] <= truncatedExon[0] and relevantElement[1] >= truncatedExon[1]:
					distanceCovered = exonLength
				#if exon completely contains element
				elif truncatedExon[0] <= relevantElement[0] and truncatedExon[1] >= relevantElement[1]:
					distanceCovered = (relevantElement[0] - truncatedExon[0] + 1) + (truncatedExon[1] - relevantElement[1] + 1)
				#otherwise find the partial overlap
				elif relevantElement[0] > truncatedExon[0]:
					distanceCovered = truncatedExon[1] - relevantElement[0] + 1
				elif relevantElement[1] < truncatedExon[1]:
					distanceCovered = relevantElement[1] - truncatedExon[0] + 1
	
				#append exon number (1 based), rejection score, distance element is covered inside exon, percentage element is 	inside exon
				formatArguments = (exons.index(truncatedExon)+1, relevantElement[2], distanceCovered, exonLength, 100.0 * distanceCovered / exonLength)
				rejectedElements.append(formatArguments)

	return rejectedElements

def getDisopredData(disopredSequencesPath, transcriptID, stopPosition):
	newData = "."
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
					if stopPosition is not None and residueNumber >= stopPosition:
						disorderedResiduesAfterPrematureStop += 1

				residueCount += 1
				if stopPosition is not None and residueNumber >= stopPosition:
					residueCountAfterPrematureStop += 1

		if stopPosition is not None:
			if residueCountAfterPrematureStop == 0:
				disorderedResiduesAfterPrematureStopPercentage = '.'
				printError("Residue count after stop position %d is zero for %s" % (stopPosition, transcriptID), False)
			else:
				disorderedResiduesAfterPrematureStopPercentage = "%.2f" % (100.0 * disorderedResiduesAfterPrematureStop / residueCountAfterPrematureStop)

			newData = "/".join([str(disorderedResidues), str(residueCount), str(disorderedResiduesAfterPrematureStop), str(residueCountAfterPrematureStop), "%.2f" % (100.0 * disorderedResidues / residueCount), disorderedResiduesAfterPrematureStopPercentage, str(stopPosition)])
		else:
			newData = "/".join([str(disorderedResidues), str(residueCount), "%.2f" % (100.0 * disorderedResidues / residueCount)])
	except IOError:
		#print("Skipping transcript %s" % (transcriptID))
		pass
	
	return newData

def getGERPelements(elementFile):
	return [(int(eline.split('\t')[0]),int(eline.split('\t')[1]), float(eline.split('\t')[3])) for eline in elementFile]

def getCodingExonIntervals(annotationIntervalPath):
	codingExonIntervals = {}
	
	for line in open(annotationIntervalPath):
		if not line.startswith("#"):
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
		#build cache file if it does not already exist using our gerprate utility
		if not os.path.exists(cachepath):
			programName = "./gerprate"
			exit_status = subprocess.check_call([programName, ratefilepath, cachepath])
			if exit_status != 0:
				printError("Exit status for following command was nonzero: %s %s %s - was aloft properly installed?" % (programName, ratefilepath, cachepath))
		return open(cachepath, "rb")
	except:
		print(ratefilepath + " could not be opened.")
		print("Exiting program.")
		sys.exit(1)

#the gerp cache files are 1-indexed
def getGerpScore(cacheFile, start, length):
	cacheFile.seek(start*4, 0)
	GERPrates = struct.unpack("<%df" % (length), cacheFile.read(length*4))
	return sum(GERPrates) / length

#binary search to find GERP element
def findGERPelementIndex(elements, start, end):
	low = 0
	high = len(elements) - 1
	while low<=high:
		mid = (low+high)//2
		if start > elements[mid][1]:
			low = mid+1
		elif end < elements[mid][0]:
			high = mid-1
		else:
			break
	if low <= high:
		return mid
	return -1
