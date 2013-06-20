import gerprate
import os

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
def findGERPelement(elements, start, end):
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
		return elements[mid]
	return None
