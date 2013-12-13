import os
import re
import struct
import subprocess
import sys
from subprocess import Popen, PIPE
import platform
import glob

def getScriptDirectory():
	return os.path.dirname(os.path.realpath(__file__))

def platformName():
	return platform.system() + "_" + platform.machine()

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
				break
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
					mergedElements.append((element[0], nextElement[1]))
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

#not used by aloft anymore but may be useful to keep around
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
					#distanceCovered = (relevantElement[0] - truncatedExon[0] + 1) + (truncatedExon[1] - relevantElement[1] + 1)
					distanceCovered = relevantElement[1] - relevantElement[0] + 1
				#otherwise find the partial overlap
				elif relevantElement[0] > truncatedExon[0]:
					distanceCovered = truncatedExon[1] - relevantElement[0] + 1
				elif relevantElement[1] < truncatedExon[1]:
					distanceCovered = relevantElement[1] - truncatedExon[0] + 1
	
				#append exon number (1 based), rejection score, distance element is covered inside exon, percentage element is 	inside exon
				formatArguments = (exons.index(truncatedExon)+1, relevantElement[2], distanceCovered, exonLength, 100.0 * distanceCovered / exonLength)
				rejectedElements.append(formatArguments)

	return rejectedElements

def getRejectionElementIntersectionPercentage(exons, truncatedExons, GERPelements, GERPelementIndex, direction):
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
	distanceCovered = 0
	intersections = []
	for truncatedExon in truncatedExons:
		for relevantElement in relevantElements:
			#check if they overlap in any way
			if truncatedExon[0] <= relevantElement[1] and truncatedExon[1] >= relevantElement[0]:
				exonLength = truncatedExon[1] - truncatedExon[0] + 1
	
				#see if element completly contains exon
				if relevantElement[0] <= truncatedExon[0] and relevantElement[1] >= truncatedExon[1]:
					distanceCovered += truncatedExon[1] - truncatedExon[0] + 1
				#if exon completely contains element
				elif truncatedExon[0] <= relevantElement[0] and truncatedExon[1] >= relevantElement[1]:
					distanceCovered += relevantElement[1] - relevantElement[0] + 1
				#otherwise find the partial overlap
				elif relevantElement[0] > truncatedExon[0]:
					distanceCovered += truncatedExon[1] - relevantElement[0] + 1
				elif relevantElement[1] < truncatedExon[1]:
					distanceCovered += relevantElement[1] - truncatedExon[0] + 1

	truncatedExonsLength = sum([block[1] - block[0] + 1 for block in truncatedExons])

	return float(distanceCovered) / truncatedExonsLength * 100.0

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

			newData = "/".join(["%.2f" % (100.0 * disorderedResidues / residueCount), disorderedResiduesAfterPrematureStopPercentage])
		else:
			newData = "/".join([str(disorderedResidues), str(residueCount), "%.2f" % (100.0 * disorderedResidues / residueCount)])
	except IOError:
		#print("Skipping transcript %s" % (transcriptID))
		pass
	
	return newData

# Get a mapping of Transcript ID's (ENST) -> Proteins ID's (ENSP)
def getTranscriptToProteinHash(transcriptToProteinFilePath):
	try:
		inputFile = open(transcriptToProteinFilePath, "r")
	except:
		printError("Failed to open %s" % (transcriptToProteinFilePath))

	transcriptToProteinHash = {}
	firstLine = True
	for line in inputFile:
		if firstLine:
			firstLine = False
		else:
			components = line.split('\t')
			if components[1].strip() and components[2].strip():
				transcriptToProteinHash[components[1]] = components[2]

	inputFile.close()
	return transcriptToProteinHash

def getFilePathMatchingPattern(filePathPattern, abortOnFatalError):
	potentialPaths = glob.glob(filePathPattern)
	if len(potentialPaths) == 0:
		printError("Failed to find file matching %s" % filePathPattern, abortOnFatalError)
	elif len(potentialPaths) > 1:
		printError("Found more than 1 potential match for %s: %s" % (filePathPattern, str(potentialPaths)), False)
	return None if len(potentialPaths) == 0 else potentialPaths[0]

def getChromosomesPfamTable(chrs, pfamDirectory, strformat, domainTypeList, domainTypeColumn=0):
	# Get a mapping of Protein ID's -> Pfam information, for each chromosome
	chromosomesPFam = {i:{} for i in domainTypeList}
	for chromosome in chrs:
		for domainType in domainTypeList:
			chromosomesPFam[domainType][chromosome] = {}

		patternPath = os.path.join(pfamDirectory, strformat % (chromosome))
		path = getFilePathMatchingPattern(patternPath, False)
		if path is None:
			printError("Couldn't find path matching %s, skipping %s" % (patternPath, chromosome), False)
			continue

		#Get rid of duplicate lines
		try:
			pipe1 = Popen(['sort', path], stdout=PIPE)
			pipe2 = Popen(['uniq'], stdin=pipe1.stdout, stdout=PIPE)
			inputFile = pipe2.stdout
		except:
			printError("Couldn't read %s, skipping %s" % (path, chromosome), False)
			continue

		linesToSkip = 2
		for lineBytes in inputFile:
			line = lineBytes.decode()
			if linesToSkip > 0:
				linesToSkip -= 1
			else:
				if not line.startswith("#"):
					components = line.split("\t")
					digitmatch = re.search("\d", components[domainTypeColumn])
					if not digitmatch:
						domainType = components[domainTypeColumn].strip()
					else:
						domainType = components[domainTypeColumn][:digitmatch.start()]
					if domainType not in domainTypeList:
						continue
					if len(components) >= 3:
						translationID = components[2].replace('(', '').replace(')', '').strip()
						if translationID in chromosomesPFam[domainType][chromosome]:
							chromosomesPFam[domainType][chromosome][translationID].append(components)
						else:
							chromosomesPFam[domainType][chromosome][translationID] = [components]

		inputFile.close()

	return chromosomesPFam

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
