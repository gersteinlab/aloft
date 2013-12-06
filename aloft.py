#!/usr/bin/env python

#Aloft is now compatible with python 2.7+ and python 3
#In order to keep compatibility, make sure to:
#do print("spam") not print "spam" since print is now a function
#do str(5) not `5` since backticks have been removed
#use // for floor division (aka integer division). 5/2 is 2.5 in python3 but 5/2 is 2 in python2 while 5//2 is 2 in both and float(5) / 2 is 2.5 in both.
#convert from binary to strings when necessary (like output from a pipe); opening files intended for binary access requires passing 'b' in file mode
#do not use xrange since it is removed
#consider that some functions (e.g, range, dict.items()) return iterators instead of lists, use list() to convert them to lists if needed
#test aloft with both python3 and python2.7
#see http://docs.python.org/3.0/whatsnew/3.0.html

import sys, os, re, string, array, datetime
from optparse import OptionParser
from subprocess import Popen, PIPE
from vat_run import *
from sequencing import *
from common import *
import argparse
import pickle
import distutils.spawn
import gzip
import vcf2bigwigbed

VERBOSE = None

def abortIfPathDoesNotExist(parser, path, shouldShowHelp=False):
    if path is not None and not os.path.exists(path):
        if shouldShowHelp:
            parser.print_help()
        printError("%s does not exist" % (path))

def abortIfCannotCreateDirectory(parser, directory):
    if not os.path.exists(directory):
        try:
            os.mkdir(directory)
        except:
            parser.print_help()
            printError("Failed to create directory %s" % (directory))

def abortIfCannotWriteFile(parser, filepath):
    try:
        newFile=open(filepath, 'w')
    except:
        parser.print_help()
        printError("%s could not be written to" % (filepath))
    return newFile

def parseCommandLineArguments():
    parser = argparse.ArgumentParser(description='Run aloft predictions. You must at least provide a VCF (via --vcf) or VAT (via --vat) input file. If you provide a VCF file, it will be ran through VAT and then through aloft. If you provide a VAT file instead, it must be sorted numerically (use vcf_sort.py for this).', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--version', action='version', version="aloft 1.0")

    parser.add_argument('--vcf', help='Path to VCF input file. This can be a compressed .gz file. If not specified, then --vat must be specified.')
    parser.add_argument('--vat', help='Path to VAT input file. If not specified, then --vcf must be specified. This file must be sorted numerically.')

    parser.add_argument('--annotation_interval', help='Path to annotation interval file for VAT', default='data/gencode.v16.pc.interval')
    parser.add_argument('--annotation_sequence', help='Path to annotation sequence file for VAT', default='data/gencode.v16.pc.fa')
    parser.add_argument('--genome', help='Path to directory containing chr*.fa files', default='data/genome/')
    
    parser.add_argument('--output', help='Path to output directory; directory is created if it does not exist', default='aloft_output/')

    parser.add_argument('--cache', help='Output to directory for cached files; directory is created if it does not exist. Only needed for networkx module.', default='cache/')

    parser.add_argument('--nmd_threshold', help='Distance from premature stop to last exon-exon junction; used to find NMD cause', type=int, default=50)

    parser.add_argument('--verbose', '-v', help='Verbose mode', action='store_true')

    parser.add_argument('--ensembl_table', help='Path to transcript to protein lookup table file', default='data/ensIDs.ens70.txt')
    parser.add_argument('--protein_features', help='Path to directory containing chr*.prot-features-ens70.txt files', default='data/ens73_all_domain_features/')
    parser.add_argument('--phosphorylation', help='Path to directory containing ptm.phosphorylation.chr*.txt files', default='data/ptm')
    parser.add_argument('--transmembrane', help='Path to directory containing transmembrane chr*.tmsigpcoilslc.ens70.txt', default='data/ens73_all_domain_features/')

    parser.add_argument('--thousandG', help='Path to 1000G file', default='data/1000G.vat')

    parser.add_argument('--ppi', help='Path to protein-protein interaction network file', default='data/BIOGRID-ORGANISM-Homo_sapiens-3.2.95.tab.txt')
    parser.add_argument('--dominant_genes', help='Path to list of dominant genes', default='data/dominantonly.list')
    parser.add_argument('--recessive_genes', help='Path to list of recessive genes', default='data/science_lofpaper_omim_recessive_filtered.list')

    parser.add_argument('--scores', help='Path to binary bw (bigwig) file containing GERP scores', default='data/All_hg19_RS.bw')
    parser.add_argument('--rates', help='Path to directory containing chr*.maf.rates files', default='data/bases/')
    parser.add_argument('--elements', help='Path to directory containing hg19_chr*_elems.txt files', default='data/elements/')

    parser.add_argument('--dNdS', help='Path to dNdS file', default='data/dNdS_avgs.txt')
    parser.add_argument('--paralogs', help='Path to paralogs file', default='data/within_species_geneparalogs.ens70')
    parser.add_argument('--LOF_score', help='Path to LOF disease scores', default='data/prob_recessive_disease_scores.txt')
    parser.add_argument('--ancestor', help='Path to directory containing homo_sapiens_ancestor_*.fa files', default='data/homo_sapiens_ancestor_GRCh37_e71/')
    parser.add_argument('--netSNP_score', help='Path to netSNP disease scores', default='data/Supplementary_Table8.20Jul2012.txt')
    parser.add_argument('--segdup', help='Path to segdup annotation file', default='data/hg19-segdup.txt')
    parser.add_argument('--annotation', help='Path to .gtf annotation file', default='data/gencode.v16.annotation.gtf')
    parser.add_argument('--exomes', help='Path to directory containing ESP6500.chr*.snps.vcf files', default='data/new_esp6500/')
    parser.add_argument('--pseudogenes', help='Path to pseudogenes file', default='data/gencode.v7.pgene.parents')
    parser.add_argument('--disopred_sequences', help='Path to disorder prediction sequences', default='data/disopred_sequences')

    args = parser.parse_args()

    global VERBOSE
    VERBOSE = args.verbose

    #safe way to test if args has an attribute named arg whose name is equal to key
    def testArgumentEquality(args, arg, key):
        if not hasattr(args, key):
            raise Exception("%s attribute does not exist!" % (key))
        return key == arg

    #Expand ~ to user's home directory for all argument paths
    for arg, path in vars(args).items():
        if path is not None and not any(map(lambda key: testArgumentEquality(args, arg, key), ['nmd_threshold', 'verbose'])):
            setattr(args, arg, os.path.expanduser(path))

    if not args.vcf and not args.vat:
        parser.print_help()
        printError("Neither a VCF or VAT file was specified. You must supply one of these as your input file")

    if args.vcf and args.vat:
        parser.print_help()
        printError("Both a VCF or VAT file were specified. You must supply only one of these as your input file, but not both")

    abortIfPathDoesNotExist(parser, args.vat)
    abortIfPathDoesNotExist(parser, args.vcf)

    abortIfCannotCreateDirectory(parser, args.output)
    try:
        import networkx
        abortIfCannotCreateDirectory(parser, args.cache)
    except:
        pass

    #Try to see if we can detect and open all input files
    for arg, path in vars(args).items():
        if not any(map(lambda key: testArgumentEquality(args, arg, key), ['vat', 'vcf', 'output', 'cache', 'nmd_threshold', 'verbose'])):
            abortIfPathDoesNotExist(parser, path, True)
            if not os.path.isdir(path):
                try:
                    f = open(path)
                    f.close()
                except:
                    printError("--%s: %s cannot be opened (insufficient read privileges?)" % (arg, path))
    return parser, args

def verifyUNIXUtility(utility):
    if distutils.spawn.find_executable(utility) is None:
        printError("Failed to find unix utility %s in PATH" % (utility))

def getAncestorData(ancespath, chromosome):
    individualAncestorPath = os.path.join(ancespath, "homo_sapiens_ancestor_%s.fa" % (chromosome))
    try:
        f=open(individualAncestorPath)
    except:
        printError("%s could not be opened... Exiting program" % (individualAncestorPath))

    f.readline()    ##first >**** line
    data = '0' + f.read().replace("\n", "")
    f.close()
    return data

def getGERPData(isSplice, GERPelements, exons, start, end, direction):
    truncatedExons = getTruncatedExons(exons, start, direction) if exons and not isSplice else None
    exonCountData = ":".join([str(len(truncatedExons)) if truncatedExons else ("." if not isSplice else "NA"), str(len(exons)) if exons else "."])

    rejectionData = "."
    elementIndex = findGERPelementIndex(GERPelements, start, end)
    if elementIndex == -1:
        elementData = "."
    else:
        elementData = str(GERPelements[elementIndex])

        if exons and truncatedExons and not isSplice:
            rejectedPercentage = getRejectionElementIntersectionPercentage(exons, truncatedExons, GERPelements, elementIndex, direction)
            rejectionData = "%.2f" % rejectedPercentage

    return elementData, rejectionData, exonCountData

def getSegDupData(vatFile, segdupPath, chrs):
    segdups={}
    segdupmax={}
    for i in chrs:
        segdups[i] = []
        segdupmax[i] = []
    if VERBOSE: print("Reading segdup information...")
    segdupfile = open(segdupPath)
    line = segdupfile.readline()
    while line.startswith("#") or line=="\n":
        line=segdupfile.readline()
    while line!="":
        data = line.split('\t')
        chr_num = data[0].split('chr')[-1]
        if '_' in chr_num:
            line = segdupfile.readline()
            continue
        segdups[chr_num].append((int(data[1]),int(data[2])))
        line = segdupfile.readline()
    for i in chrs:
        segdups[i] = sorted(segdups[i])
        maxsofar = 0
        for interval in segdups[i]:
            maxsofar = max(interval[1], maxsofar)
            segdupmax[i].append(maxsofar)
    segdupdata=[]
    vatFile.seek(0)
    line = vatFile.readline()
    while line.startswith("#") or line=="\n":
        segdupdata.append('')
        line=vatFile.readline()
    if VERBOSE: print('Calculating segdup overlaps...')
    while line!="":
        data = line.split('\t')
        chr_num = data[0].split('chr')[-1]
        start = int(data[1])
        length = len(data[3])
        end = start + length-1  ##inclusive endpoint

        if chr_num not in segdups:
            line = vatFile.readline()
            continue

        ##find right endpoint of interval search 
        low = 0; high = len(segdups[chr_num])-1
        while low<=high:
            mid = (low+high)//2
            if end<segdups[chr_num][mid][0]:
                high = mid-1
            elif mid == len(segdups[chr_num])-1 or end<segdups[chr_num][mid+1][0]:
                break
            else:
                low = mid+1
        right = mid
            
        ##find left endpoint of interval search
        low = 0; high = len(segdups[chr_num])-1
        while low<=high:
            mid = (low+high)//2
            if start>segdups[chr_num][mid][1] and start>segdupmax[chr_num][mid]:
                low = mid+1
            elif mid==0:
                break
            elif start>segdups[chr_num][mid-1][1] and start>segdupmax[chr_num][mid-1]:
                break
            else:
                high = mid-1
        left = mid
        
        overlaps = []
        for interval in segdups[chr_num][left:right+1]:
            ##compare bigger of left enpoints to smaller of right endpoints
            if max(start, interval[0]) <= min(end, interval[1]):
                overlaps.append(interval)
        segdupdata.append(str(overlaps))

        line = vatFile.readline()

    segdupfile.close()
    return segdupdata

#Returns a Pfam vcf-formatted description, and a verbose description
#The vcf-formatted description is in format Pfam_ID:domain_length:max_domain_percent_lost:number_pfams_in_domain:number_pfams_in_truncation
#The verbose description returned is a) a series of concatenated Pfam_ID:domain_length:percent_lost, and a series of concatenated domain_id_lost:domain_length 
##domain value=amino acid coordinate of premature stop
def getPfamDescription(transcriptToProteinHash, chromosome, transcriptID, domainValue, chromosomesPFam, domainType):
    pfamsMatched = []
    maxPercentageLostPfamIndex = -1
    pfamDescription = ""
    pfamVerboseDescription = None

    chromosomesPFam = chromosomesPFam[domainType]

    if transcriptID in transcriptToProteinHash and transcriptToProteinHash[transcriptID] in chromosomesPFam[chromosome]:
        pfamComponentsList = chromosomesPFam[chromosome][transcriptToProteinHash[transcriptID]]
        domainsLost = ""
        numberOfDomainsLost = 0
        for pfamComponents in pfamComponentsList:
            domainComponents = pfamComponents[1].split("-")
            domainStart = int(domainComponents[0])
            domainEnd = int(domainComponents[1])
            domainLength = domainEnd - domainStart + 1
            if domainValue >= domainStart and domainValue <= domainEnd:
                domainLengthLost = domainEnd - domainValue + 1
                domainPercentLost = domainLengthLost * 100.0 / domainLength

                if not pfamVerboseDescription:
                    pfamVerboseDescription = ""

                pfamVerboseDescription += ":%s:%d:%.2f" % (pfamComponents[0], domainLength, domainPercentLost)

                #Find largest percentage lost pfam
                if maxPercentageLostPfamIndex < 0 or domainPercentLost > float(pfamsMatched[maxPercentageLostPfamIndex].split(":")[3]):
                    maxPercentageLostPfamIndex = len(pfamsMatched)

                pfamsMatched.append(":%s:%d:%.2f" % (pfamComponents[0], domainLength, domainPercentLost))

            elif domainValue < domainStart:
                domainsLost += ":" + pfamComponents[0] + ":" + str(domainLength)
                numberOfDomainsLost += 1

        if maxPercentageLostPfamIndex >= 0:
            pfamDescription = pfamsMatched[maxPercentageLostPfamIndex] + ":" + str(len(pfamsMatched))

        if pfamDescription == "":
            pfamDescription = ":NO_"+domainType+":NA:NA:0"
            pfamVerboseDescription = "NO_"+domainType

        verboseDomainsLost = str(domainsLost)
        if verboseDomainsLost.startswith(":"):
            verboseDomainsLost = verboseDomainsLost[1:]

        if verboseDomainsLost == "":
            verboseDomainsLost = "NO_"+domainType

        if pfamVerboseDescription.startswith(":"):
            #Remove beginning colon
            pfamVerboseDescription = pfamVerboseDescription[1:]

        pfamVerboseDescription = [pfamVerboseDescription, verboseDomainsLost]

        #Add number of domains lost
        pfamDescription += ":" + str(numberOfDomainsLost)

    #If no ENST or ENSP ID could be matched
    if pfamDescription == "":
        pfamDescription = ":NA_"+domainType+":NA:NA:0:0"
    if pfamVerboseDescription is None:
        pfamVerboseDescription = ["NA_"+domainType, "NA_"+domainType]

    #pfamShortDescription = "NO" if pfamDescription.startswith(":NO") else "YES"
    if pfamDescription.startswith(":NO"):
        pfamShortDescription = 'NO'
    elif pfamDescription.startswith(":NA"):
        pfamShortDescription = 'NA'
    else:
        pfamShortDescription = 'YES'

    return pfamDescription, pfamShortDescription, pfamVerboseDescription

def getGenomeSequences(genomePath, chromosome):
    individualSequencePath = os.path.join(genomePath, "chr%s.fa" % (chromosome))
    try:
        f=open(individualSequencePath)
    except:
        printError("%s could not be opened" % (individualSequencePath))

    f.next() ##first >chr* line
    genomeSequences = '0' + ''.join([line.strip() for line in f])
    f.close()

    return genomeSequences

def getCDSAndExonDictionaries(annotationPath, chrs):
    CDS={}; exon={}; stop_codon={}  ##{chr_num: {transcript: [(a,b),(c,d)..] } }
    transcript_strand={}            ##{transcript_id:+ or -}
    for chr_num in chrs:
        CDS[chr_num]={}
        exon[chr_num]={}
        stop_codon[chr_num]={}
    CDS['M']={}
    exon['M']={}
    stop_codon['M']={}

    annotfile = open(annotationPath)

    ##begin going through actual annotation data
    oldtr = ""  ##last seen transcript
    oldchr = "" ##chr num of last seen transcript
    tlines = [] ##all split CDS lines in oldtr
    for line in annotfile:
        if line.startswith("#"):
            continue
        data = line.strip().split('\t')
        chr_num=data[0].split('chr')[-1]
        annottype = data[2]
        
        if annottype!='exon' and annottype!='transcript' and annottype!='CDS' and annottype!='stop_codon':
            continue
        
        if annottype=='transcript':
            transcript = data[8].split(';')[1].split('"')[1]
            transcript_strand[transcript]=data[6]
            if oldtr!="":
                if len(tlines)>0:
                    if transcript_strand[oldtr]=='+':
                        oldsort = sorted(tlines, key=lambda s: int(s[3]))
                        first = oldsort[0]
                        CDS[oldchr][oldtr].append((int(first[3])+int(first[7]), int(first[4])))
                        for CDSline in oldsort[1:]:
                            CDS[oldchr][oldtr].append((int(CDSline[3]),int(CDSline[4])))
                    else:
                        oldsort = sorted(tlines, key=lambda s: int(s[3]), reverse=True)
                        first = oldsort[0]
                        CDS[oldchr][oldtr].append((int(first[3]), int(first[4])-int(first[7])))
                        for CDSline in oldsort[1:]:
                            CDS[oldchr][oldtr].append((int(CDSline[3]),int(CDSline[4])))
            oldtr = transcript
            oldchr = chr_num
            tlines=[]
            exon[chr_num][transcript] = []
            CDS[chr_num][transcript] = []
        else:  ## then is either exon or CDS or stop codon
            begin = int(data[3])
            end = int(data[4])
            if data[2]=='exon':  
                exon[chr_num][transcript].append((begin, end))  ##
            elif data[2]=='CDS':
                tlines.append(data)  ##append data for reanalysis (mRNA_start_NF cases)
            else:  ##stop codon
                stop_codon[chr_num][transcript] = (begin,end)

    #this is really just copy of code above, todo: clean this up
    if len(tlines)>0 and oldtr != '':
        if transcript_strand[oldtr]=='+':
            oldsort = sorted(tlines, key=lambda s: int(s[3]))
            first = oldsort[0]
            CDS[oldchr][oldtr].append((int(first[3])+int(first[7]), int(first[4])))
            for CDSline in oldsort[1:]:
                CDS[oldchr][oldtr].append((int(CDSline[3]),int(CDSline[4])))
        else:
            oldsort = sorted(tlines, key=lambda s: int(s[3]), reverse=True)
            first = oldsort[0]
            CDS[oldchr][oldtr].append((int(first[3]), int(first[4])-int(first[7])))
            for CDSline in oldsort[1:]:
                CDS[oldchr][oldtr].append((int(CDSline[3]),int(CDSline[4])))

    annotfile.close()

    return transcript_strand, CDS, exon, stop_codon

def get1000GChromosomeInfo(thousandGPath):
    thousandGChromosomeInfo = {}

    if thousandGPath.endswith(".gz"):
        thousandGFile = gzip.open(thousandGPath)
    else:
        thousandGFile = open(thousandGPath)

    for thousandGLine in thousandGFile:
        if not thousandGLine.startswith("#"):
            thousandGLineComponents = thousandGLine.rstrip("\n").split("\t")
            thousandGChromosomeNumber = thousandGLineComponents[0].replace("chr", "")
            if not (thousandGChromosomeNumber in thousandGChromosomeInfo):
                thousandGChromosomeInfo[thousandGChromosomeNumber] = {}
            
            thousandGChromosomeInfo[thousandGChromosomeNumber][int(thousandGLineComponents[1])] = thousandGLineComponents[7]

    return thousandGChromosomeInfo

def getPPINetwork(networkx, ppiPath):
    ppi = networkx.Graph()
    ppifile = open(ppiPath)
    ppifile.readline()
    for line in ppifile:
        data = line.split('\t')
        ppi.add_edge(data[2], data[3])
    ppifile.close()
    return ppi

#This function assumes the key is the 0'th column
def getScores(filepath, scoreColumnIndex):
    scores = {}
    scorefile = open(filepath)
    scorefile.readline()
    for line in scorefile:
        data = line.strip().split("\t")
        scores[data[0].upper()] = data[scoreColumnIndex]
    scorefile.close()
    return scores

def getPseudogeneData(pseudogenesPath):
    numpseudogenes = {}     ##{parent transcript: # of assoc. pseudogenes}
    pseudogenesfile = open(pseudogenesPath)
    pseudogenesfile.readline()
    for line in pseudogenesfile:
        tx = line.split('\t')[6]
        if tx in numpseudogenes:
            numpseudogenes[tx] = numpseudogenes[tx]+1
        else:
            numpseudogenes[tx] = 1
    pseudogenesfile.close()
    return numpseudogenes

def getParalogData(paralogsPath):
    paralogs = {}     ##{ENSG ID (without . subclassifier): set(assoc. paralogs)}
    paralogsfile = open(paralogsPath)
    paralogsfile.readline()
    for line in paralogsfile:
        id1 = line.split('\t')[0]
        id2 = line.split('\t')[1]
        if id1 not in paralogs:
            paralogs[id1] = set()
        paralogs[id1].add(id2)
        if id2 not in paralogs:
            paralogs[id2] = set()
        paralogs[id2].add(id1)
    paralogsfile.close()
    return paralogs

def getdNdSData(dNdSPath):
    dNdSmacaque = {}      ##{ENST ID (without . subclassifier): dN/dS (STRING)}
    dNdSmouse = {}
    dNdSfile = open(dNdSPath)
    dNdSfile.readline()
    for line in dNdSfile:
        data = line.strip().split('\t')
        tx = data[1]
        dNdSmacaque[tx] = data[2]
        dNdSmouse[tx] = data[3]
    dNdSfile.close()
    return dNdSmacaque, dNdSmouse

def calculateExomeCoordinate(component):
    values = component.split("=")[1].split(",")
    if (int(values[0]) + int(values[1])) == 0:
        return 0.0
    return int(values[0]) * 1.0 / (int(values[0]) + int(values[1]))

def getESP6500ExomeChromosomeInfo(exomesPath, chromosome):
    exomesChromosomeInfo = {}
    #exomePath = os.path.join(exomesPath, 'ESP6500.chr%s.snps.vcf' % (chromosome))
    exomePath = os.path.join(exomesPath, 'ESP6500SI-V2-SSA137.updatedRsIds.chr%s.snps_indels.vcf' % (chromosome))
    try:
        exomeInputFile = open(exomePath)
    except:
        printError("Couldn't read %s, skipping.." % (exomePath), False)
        exomeInputFile = None
    if exomeInputFile:
        for exomeLine in exomeInputFile:
            if not exomeLine.startswith("#"):
                exomeLineComponents = exomeLine.split("\t")
                
                x = "NA"
                y = "NA"
                z = "NA"
                for component in exomeLineComponents[7].split(";"):
                    if component.startswith('EA_AC='):
                        x = "%.4f" % (calculateExomeCoordinate(component))
                    elif component.startswith('AA_AC='):
                        y = "%.4f" % (calculateExomeCoordinate(component))
                    elif component.startswith('TAC='):
                        z = "%.4f" % (calculateExomeCoordinate(component))
                        
                exomesChromosomeInfo[int(exomeLineComponents[1])] = ("%s,%s,%s" % (x, y, z))
        exomeInputFile.close()
    return exomesChromosomeInfo

def parsePPI(networkx, ppi, ppiHash, hashKey, gene_name, genes):
    dist = None
    if gene_name in ppiHash[hashKey]:
        dist = ppiHash[hashKey][gene_name]
    else:
        for gene in genes:
            if gene != gene_name and gene in ppi and networkx.has_path(ppi, gene_name, gene):
                shortestPathLength = networkx.shortest_path_length(ppi, gene_name, gene)
                if dist is None:
                    dist = shortestPathLength
                else:
                    dist = min(dist, shortestPathLength)
            ppiHash[hashKey][gene_name] = dist

    numberOfNeighbors = sum(1 for gene in genes if gene != gene_name and gene in ppi.neighbors(gene_name))
    return dist, numberOfNeighbors

def findNMDForIndelsAndPrematureStop(nmdThreshold, data, chr_num, transcript, start, end, exon, stop_codon, genomeSequences, CDS, subst, transcript_strand):
    nmdHash = {"NMD" : None, 'splice1' : None, 'splice2' : None, 'canonical' : None, 'newCDSpos' : None, 'stopCDS' : None, 'nextATG' : None, 'incrcodingpos' : None, 'issinglecodingexon' : None}

    l = sorted(CDS[chr_num][transcript])
    if len(l)==0:
        return nmdHash
    nmdHash['issinglecodingexon'] = "YES" if len(l)==1 else "NO"
    numberOfExonsHash = sorted(exon[chr_num][transcript])
    CDSseq = ''; exonseq = ''
    CDSprec = []; exonprec = []             ## prec holds # preceding nucleotides
    ispositivestr = transcript_strand[transcript]=='+'
    
    ## build spliced exon and CDS sequences and maintain coordinate wrt transcript
    tot = 0
    for j in range(0,len(l)):
        if ispositivestr:
            i=j
            CDSseq+=genomeSequences[l[i][0]:l[i][1]+1].upper()
        else:
            i=len(l)-j-1
            CDSseq+=compstr(genomeSequences[l[i][0]:l[i][1]+1].upper())
        CDSprec.append(tot)              ## stores in index i
        tot += l[i][1]+1-l[i][0]
    ## add on STOP sequence if annotated
    try:
        s = stop_codon[chr_num][transcript]
    except:
        s=(2,0)
    if ispositivestr:
        CDSseq+=genomeSequences[s[0]:s[1]+1].upper()
    else:
        CDSseq+=compstr(genomeSequences[s[0]:s[1]+1].upper())
    
    tot = 0
    for j in range(0,len(numberOfExonsHash)):
        if ispositivestr:
            i=j
            exonseq+=genomeSequences[numberOfExonsHash[i][0]:numberOfExonsHash[i][1]+1].upper()
        else:
            i=len(numberOfExonsHash)-j-1
            exonseq+=compstr(genomeSequences[numberOfExonsHash[i][0]:numberOfExonsHash[i][1]+1].upper())
        exonprec.append(tot)            ## stores in index i
        tot += numberOfExonsHash[i][1]+1-numberOfExonsHash[i][0]
    
    ##build coding exons IN ORDER OF TRANSLATION, i.e. start->stop
    coding_exons = []           ## flag coding exons (corresponds to exonpos)
    CDS2ex = {}                 ## maps CDSpos to exonpos    
    for i in range(0,len(numberOfExonsHash)):   ## i = exonpos
        k = i if ispositivestr else len(numberOfExonsHash)-i-1  ## k = exonindex
        coding_exons.append(0)
        for j in range(0,len(l)):   ## j = CDSindex
            if l[j][0]>=numberOfExonsHash[k][0] and l[j][1]<=numberOfExonsHash[k][1]:
                coding_exons[i] = 1
                j2 = j if ispositivestr else len(l)-j-1     ## j2 = CDSpos
                CDS2ex[j2]=i
                break
                
    ncodingexons = sum(coding_exons)    ## number of coding exons
    try:
        UTR=len(numberOfExonsHash)-(coding_exons.index(1)+ncodingexons) ## number of 3'UTR exons
    except:     ## no coding exons
        UTR=0
    
    ## find CDS and exon interval numbers
    flag1=0
    flag2=0
    CDSpos=-1        ## this gives the CDSpos 0-based: i.e. first CDS is index 0
    for i in range(0,len(l)):
        if start>=l[i][0] and start<=l[i][1]:
            CDSpos = i if ispositivestr else len(l)-i-1
        if start>=l[i][0] and end<=l[i][1]:
            flag1 = 1   ##indel is completely contained in CDS
            break
    
    exonpos=-1       ## this gives the exonpos also 0-based
    for i in range(0,len(numberOfExonsHash)):
        if start>=numberOfExonsHash[i][0] and start<=numberOfExonsHash[i][1]:
            exonpos = i if ispositivestr else len(numberOfExonsHash)-i-1
        if start>=numberOfExonsHash[i][0] and end<=numberOfExonsHash[i][1]:
            flag2 = 1   ##indel is completely contained in exon
            break
    if CDSpos==-1 or exonpos==-1:   ##start position of indel was not in ANY intervals
        nmdHash['NMD'] = "no exons or no CDS containing start of indel"
        return nmdHash
    
    exonindex= exonpos if ispositivestr else len(numberOfExonsHash)-exonpos-1
    CDSindex= CDSpos if ispositivestr else len(l)-CDSpos-1
    
    codingpos = exonpos-coding_exons.index(1)      ## this gives coding exon position 0-based
    
    diff = len(subst)-len(data[3])
    if ispositivestr:
        ## 1-based position of indel in CDS coordinates
        newCDSpos = CDSprec[CDSpos] + start - l[CDSindex][0] + 1
        ## 1-based position of indel in exon coordinates
        newexonpos = exonprec[exonpos] + start - numberOfExonsHash[exonindex][0] + 1
    else:
        newCDSpos = CDSprec[CDSpos] + l[CDSindex][1] - start + 1
        newexonpos = exonprec[exonpos] + numberOfExonsHash[exonindex][1] - start + 1
    ## # of exon nucleotides before e-e junction
    if newexonpos-1<exonprec[-1]:       ##indel being before last e-e junction shifts e-e position
        juncpos = exonprec[-1]+diff     ##WRONG IF START POSITION IS LAST NUCLEOTIDE BEFORE E-E AND
                                        ##IS DELETION (WOULD BE SPLICE OVERLAP)
    else:                               ##indel is after last junction, position unchanged
        juncpos = exonprec[-1]
    
    nmdHash['newCDSpos'] = newCDSpos
    
    lastindex = -1 if ispositivestr else 0
    indeltoend = exonprec[-1]+numberOfExonsHash[lastindex][1]-numberOfExonsHash[lastindex][0]+1 + diff - (newexonpos-1) ##CHECK THIS EXTRA +1
    if flag1==0:
        nmdHash['NMD'] = "no CDS regions completely containing variant"
        return nmdHash

    if flag2==0:
        nmdHash['NMD'] = "no exon regions completely containing variant"
        return nmdHash
        
    if ispositivestr:
        modCDSseq = CDSseq[0:newCDSpos-1]
        modCDSseq += subst
        modCDSseq += CDSseq[newCDSpos-1+len(data[3]):]
    else:
        modCDSseq = CDSseq[0:newCDSpos-len(data[3])]
        modCDSseq += compstr(subst)
        modCDSseq += CDSseq[newCDSpos:]
    ref_aa = translate_aa(CDSseq)
    alt_aa = translate_aa(modCDSseq)
    
    try:
        nextATG = str(3*(alt_aa[:1].index('M')+1))
    except:
        nextATG = 'NA'

    nmdHash['nextATG'] = nextATG
    
    ## # of CDS nucleotides before stop codon in alternate sequence
    try:
        stopCDS = 3*alt_aa.index('*')
        nmdHash['stopCDS'] = stopCDS
    except:
        nmdHash['NMD'] = "No stop codon found in alt_aa"
        return nmdHash
    
    ## stopexon is # of exon nucleotides preceding first nucleotide of stop codon
    ## increxon is the exon position (not exon index) where the new stop occurs
    
    ## if is in very last CDS or STOP is in current CDS
    if CDSpos==len(l)-1 or (stopCDS>=CDSprec[CDSpos] and stopCDS<CDSprec[CDSpos+1]+diff):
        increxon = exonpos
        if ispositivestr:
            stopexon = exonprec[exonpos]+l[CDSpos][0]-numberOfExonsHash[exonpos][0]+stopCDS-CDSprec[CDSpos]
        else:
            stopexon = exonprec[exonpos]+numberOfExonsHash[exonindex][1]-l[CDSindex][1]+stopCDS-CDSprec[CDSpos]
    else:
        incrCDS = CDSpos
        while incrCDS<len(l):
            increxon = CDS2ex[incrCDS]
            if incrCDS==len(l)-1 or (stopCDS>=CDSprec[incrCDS]+diff and stopCDS<CDSprec[incrCDS+1]+diff):
                if ispositivestr:
                    stopexon = exonprec[increxon]+l[incrCDS][0]-numberOfExonsHash[increxon][0]+stopCDS-(CDSprec[incrCDS]+diff)
                else:
                    stopexon = exonprec[increxon]+numberOfExonsHash[len(numberOfExonsHash)-increxon-1][1]-l[len(l)-incrCDS-1][1]+stopCDS-(CDSprec[incrCDS]+diff)
                break
            incrCDS+=1
    incrcoding = sum(coding_exons[:increxon])   ##incrcoding is the coding exon position where new stop occurs
    if ncodingexons == 1:
        incrcodingpos = 'single'
    elif incrcoding==0:
        incrcodingpos = 'first'
    elif incrcoding==ncodingexons-1:
        incrcodingpos = 'last'
    else:
        incrcodingpos = 'middle'

    nmdHash['incrcodingpos'] = incrcodingpos
    
    ## distances are calculated as follows: TAG_ _| is 5 from exon-exon junction/end of transcript etc.
    ## end of transcript is denoted as end of last exon.
    
    ##number of nucleotides in all exons
    transcriptend = exonprec[-1]+numberOfExonsHash[lastindex][1]-numberOfExonsHash[lastindex][0]+1 + diff    ## in exon coordinates
    stoptoend = transcriptend - stopexon
    stoptojunc = juncpos - stopexon
    indeltoend = transcriptend - (newexonpos-1)
    nmdHash['NMD'] = 'YES' if stoptojunc >= nmdThreshold else 'NO'
    
    ##exon index where new stop occurs
    increxonindex = increxon if ispositivestr else len(numberOfExonsHash)-increxon-1
    
    if increxon==0:
        splice1='.'     ## 5' flanking splice site (acceptor)
    else:
        if ispositivestr:
            splice1=genomeSequences[numberOfExonsHash[increxonindex][0]-2:numberOfExonsHash[increxonindex][0]].upper()
        else:
            splice1=compstr(genomeSequences[numberOfExonsHash[increxonindex][1]+1:numberOfExonsHash[increxonindex][1]+3].upper())                
    if increxon==len(numberOfExonsHash)-1:
        splice2='.'     ## 3' flanking splice site (donor)
    else:
        if ispositivestr:
            splice2=genomeSequences[numberOfExonsHash[increxon][1]+1:numberOfExonsHash[increxon][1]+3].upper()
        else:
            splice2=compstr(genomeSequences[numberOfExonsHash[increxonindex][0]-2:numberOfExonsHash[increxonindex][0]].upper())
        
    canonical = (splice1=='AG' or splice1=='.') and (splice2=='GT' or splice2=='.')
    canonical = 'YES' if canonical else 'NO'

    nmdHash['splice1'] = splice1
    nmdHash['splice2'] = splice2
    nmdHash['canonical'] = canonical

    return nmdHash

#not exactly sure what this function does so not exactly sure what to call it
def searchInSplices(chr_num, transcript, genomeSequences, ispositivestr, start, CDS, subst):
    newData = {'found' : None, 'new' : None, 'acceptor' : None, 'donor' : None, 'intronlength' : None}
    l = sorted(CDS[chr_num][transcript], reverse= not ispositivestr)
    found = False
    end = 0  ##0 for end toward smaller basepair number, 1 for other end
    for i in range(0,len(l)):
        r = l[i]
        if start-r[1] in [1,2]:
            end = 1
            found = True
            break
        elif r[0]-start in [1,2]:
            end = 0
            found = True
            break

    newData['found'] = found
    if not found:
        spliceOutputFile.write('\t'+'\t'.join(outdata[i] for i in ["shortest_path_to_recessive_gene", "recessive_neighbors"]))
        spliceOutputFile.write("\tCDS match not found: pos="+str(start)+' transcript='+transcript+'\n')
        return newData

    if ispositivestr:
        if (end==0 and i==0) or (end==1 and i==len(l)-1):
            return newData
        if end==0:
            acceptor = genomeSequences[l[i][0]-2:l[i][0]].upper()
            if start==l[i][0]-2:
                new = (1, subst+acceptor[1])
            else:
                new = (1, acceptor[0]+subst)
            donor = genomeSequences[l[i-1][1]+1:l[i-1][1]+3].upper()
            intronlength = l[i][0]-l[i-1][1]-1
        elif end==1:
            acceptor = genomeSequences[l[i+1][0]-2:l[i+1][0]].upper()
            donor = genomeSequences[l[i][1]+1:l[i][1]+3].upper()
            if start==l[i][1]+1:
                new = (0, subst+donor[1])
            else:
                new = (0, donor[0]+subst)
            intronlength = l[i+1][0]-l[i][1]-1
    else:   ##not ispositivestr
        if (end==1 and i==0) or (end==0 and i==len(l)-1):
            return newData
        if end==0:
            donor = genomeSequences[l[i][0]-2:l[i][0]].upper()
            acceptor = genomeSequences[l[i+1][1]+1:l[i+1][1]+3].upper()
            if start==l[i][0]-2:
                new = (0, subst+donor[1])
            else:
                new = (0, donor[0]+subst)
            intronlength = l[i][0]-l[i+1][1]-1
        elif end==1:
            donor = genomeSequences[l[i-1][0]-2:l[i-1][0]].upper()
            acceptor = genomeSequences[l[i][1]+1:l[i][1]+3].upper()
            if start==l[i][1]+1:
                new = (1, subst+acceptor[1])
            else:
                new = (1, acceptor[0]+subst)
            intronlength = l[i-1][0]-l[i][1]-1
        donor = compstr(donor.upper())
        acceptor = compstr(acceptor.upper())
        new = (new[0],compstr(new[1].upper()))

    newData['donor'] = donor
    newData['acceptor'] = acceptor
    newData['new'] = new
    newData['intronlength'] = intronlength

    return newData

def getMatchingNagnagnagPositions(genomeSequences, start, ispositivestr):
    #NAGN <snp>AG NAG
    nagnagSequence = genomeSequences[start-4:start+5].upper()
    if not ispositivestr:
        nagnagSequence = compstr(nagnagSequence)

    nagNagPositions = []

    if nagnagSequence[1:3] == 'AG':
        if ispositivestr:
            nagNagPositions.append(start-3)
        else:
            nagNagPositions.append(start+3)
    
    if nagnagSequence[7:9] == 'AG':
        if ispositivestr:
            nagNagPositions.append(start+3)
        else:
            nagNagPositions.append(start-3)

    return nagNagPositions

#find stop position to use for GERP calculations by using relative stop position in CDS and mapping it against the coding exon intervals for the transcript
def calculateAbsolutePosition(lofPosition, codingExonIntervals, direction):
    distanceLeft = lofPosition
    exonIntervals = codingExonIntervals if direction == '+' else reversed(codingExonIntervals)
    absoluteStopPosition = 0
    for exonInterval in exonIntervals:
        zeroBasedExonInterval = (exonInterval[0]-1, exonInterval[1])
        exonDistance = zeroBasedExonInterval[1] - zeroBasedExonInterval[0]
        if distanceLeft > exonDistance:
            distanceLeft -= exonDistance
        else:
            if direction == '+':
                absoluteStopPosition = distanceLeft + zeroBasedExonInterval[0]
            else:
                absoluteStopPosition = zeroBasedExonInterval[1] - distanceLeft + 1
            break
    return absoluteStopPosition

def main():
    startProgramExecutionTime = datetime.datetime.now()

    #getChromosomesPfamTable() function relies on these two dependencies
    verifyUNIXUtility('sort')
    verifyUNIXUtility('uniq')

    parser, args = parseCommandLineArguments()

    if args.vcf:
        #run VAT
        vatPath = os.path.join(args.output, os.path.basename(args.vcf) + ".vat")
        run_vat([sys.argv[0], args.vcf, vatPath, args.annotation_interval, args.annotation_sequence], VERBOSE)
    else:
        vatPath = args.vat
    
    if VERBOSE: print("Running ALoFT on %s" % (vatPath) + "\n")
    
    try:
        vatFile = open(vatPath)
    except:
        printError("Failed to read %s" % (vatPath))
    
    tabbedOutputLofPath = os.path.join(args.output, os.path.basename(vatPath) + ".aloft.lof")
    tabbedOutputSplicePath = os.path.join(args.output, os.path.basename(vatPath) + ".aloft.splice")
    vcfOutputPath = os.path.join(args.output, os.path.basename(vatPath) + ".aloft.vcf")

    lofOutputFile = abortIfCannotWriteFile(parser, tabbedOutputLofPath)
    spliceOutputFile = abortIfCannotWriteFile(parser, tabbedOutputSplicePath)
    vcfOutputFile = abortIfCannotWriteFile(parser, vcfOutputPath)
    
    chrs = [str(i) for i in range(1, 23)] + ['X', 'Y']

    gerpScoresHash = {}
    bigWigAverageOverBedPath = os.path.join(os.path.join('bigwig-bin', platform.system() + "_" + platform.machine()), 'bigWigAverageOverBed')
    bigWigTabOutputPath = os.path.join(args.output, 'bigwig.tab')
    bigWigBedInputPath = os.path.join(args.output, 'bigwig_input.bed')
    bigWigBedOutputPath = os.path.join(args.output, 'bigwig_output.bed')
    try:
        vcf2bigwigbed.writeBed(vatPath, bigWigBedInputPath)
        if VERBOSE: print("Running bigWigAverageOver...")
        subprocess.check_output([bigWigAverageOverBedPath, args.scores, bigWigBedInputPath, bigWigTabOutputPath, '-bedOut=%s' % bigWigBedOutputPath], stderr=subprocess.STDOUT)

        with open(bigWigBedOutputPath) as bigWigFile:
            for line in bigWigFile:
                data = line.strip().split("\t")
                chromosome = data[0]
                if chromosome not in gerpScoresHash:
                    gerpScoresHash[chromosome] = {}

                position = int(data[1])+1 #to 1 based coordinate
                score = float(data[4])
                gerpScoresHash[chromosome][position] = score

        os.remove(bigWigBedInputPath)
        os.remove(bigWigTabOutputPath)
        os.remove(bigWigBedOutputPath)
    except:
        printError("Failed to call bigWigAverageOverBed")
    
    #Load exon intervals from .interval file, used later for intersecting with gerp elements
    codingExonIntervals = getCodingExonIntervals(args.annotation_interval)
    
    segdupdata = getSegDupData(vatFile, args.segdup, chrs)
    
    if VERBOSE:
        print('Building CDS and exon dictionaries...')
        startTime = datetime.datetime.now()
    
    transcript_strand, CDS, exon, stop_codon = getCDSAndExonDictionaries(args.annotation, chrs)
    
    if VERBOSE:
        print(str((datetime.datetime.now() - startTime).seconds) + " seconds.")
        print('Begin ALoFT Calculations and Write-Out (this may take a while)...')
            
    transcriptToProteinHash = getTranscriptToProteinHash(args.ensembl_table)

    ##{'1':{'ENSP...':'PF...\t4-25\t(ENSP...)'}, '2':{...}, ...}
    #proteinFeaturesList = list(getChromosomesPfamTable(chrs, args.protein_features, r"chr%s.prot-features-ens70.txt", ["PF", "SSF", "SM"]).items())
    #phosphorylationTags = ["ACETYLATION", "DI-METHYLATION", "METHYLATION", "MONO-METHYLATION", "O-GlcNAc", "PHOSPHORYLATION", "SUMOYLATION", "TRI-METHYLATION", "UBIQUITINATION"]
    #phosphorylationFeaturesList = list(getChromosomesPfamTable(chrs, args.phosphorylation, r"ptm.phosphosite.chr%s.txt", phosphorylationTags, 3).items())
    #transmembraneFeaturesList = list(getChromosomesPfamTable(chrs, args.transmembrane, r"chr%s.tmsigpcoilslc.ens70.txt", ["Tmhmm", "Sigp"]).items())

    proteinFeaturesList = list(getChromosomesPfamTable(chrs, args.protein_features, r"%s.ens73.alldomainfeatures.txt", ["PF", "SSF", "SM"]).items())
    phosphorylationTags = ["ACETYLATION", "DI-METHYLATION", "METHYLATION", "MONO-METHYLATION", "O-GlcNAc", "PHOSPHORYLATION", "SUMOYLATION", "TRI-METHYLATION", "UBIQUITINATION"]
    phosphorylationFeaturesList = list(getChromosomesPfamTable(chrs, args.phosphorylation, r"ptm.phosphosite.chr%s.txt", phosphorylationTags, 3).items())
    transmembraneFeaturesList = list(getChromosomesPfamTable(chrs, args.transmembrane, r"%s.ens73.alldomainfeatures.txt", ["Tmhmm", "Sigp"]).items())

    chromosomesPFam = dict(proteinFeaturesList + phosphorylationFeaturesList + transmembraneFeaturesList)

    #Scan 1000G file
    if VERBOSE: print("Scanning 1000G file")
    thousandGChromosomeInfo = get1000GChromosomeInfo(args.thousandG)
    
    ppi = None
    ppiHashPath = None

    try:
        import networkx
        if VERBOSE: print("Reading PPI network")
        ppi = getPPINetwork(networkx, args.ppi)
        ppiHashPath = os.path.join(args.cache, "ppi")
    except:
        if VERBOSE: print("Skipping PPI network, since no networkx module found")

    ppiHash = None
    if ppiHashPath is not None:
        #ppiHash will contain cached values of shortest paths to genes
        if os.path.exists(ppiHashPath):
            ppiHash = pickle.load(open(ppiHashPath, "rb"))
        else:
            ppiHash = {"dgenes" : {}, "rgenes" : {}}
    
    if VERBOSE: print("Reading recessive genes list")
    rgenes = [line.strip() for line in open(args.recessive_genes)]
    
    if VERBOSE: print("Reading dominant genes list")
    dgenes = [line.strip() for line in open(args.dominant_genes)]
    
    if VERBOSE: print("Reading pseudogene data")
    numpseudogenes = getPseudogeneData(args.pseudogenes)
    
    if VERBOSE: print("Reading paralog data")
    paralogs = getParalogData(args.paralogs)
    
    if VERBOSE: print("Reading dNdS data")
    dNdSmacaque, dNdSmouse = getdNdSData(args.dNdS)

    ptmParams = ["ACETYLATION", "DI-METHYLATION", "METHYLATION", "MONO-METHYLATION", "O-GlcNAc","PHOSPHORYLATION", "SUMOYLATION", "TRI-METHYLATION", "UBIQUITINATION"]

    #params for PF, SSF, SM, etc
    #this variable could use a better name since it's not just PFAM, but not sure what to call it
    pfamParams = ["PF", "SSF", "SM", "Tmhmm", "Sigp"]

    pfamParamsWithTruncations = sum([[param, param + "truncated"] for param in pfamParams + ptmParams], []) #using sum to flatten the list

    ##list of output parameters for LOF and splice variants
    basicparams = ["gene", "gene_id", "partial/full", "transcript", "transcript_length", "longest_transcript?"]
    LOFparams = ["shortest_path_to_recessive_gene", "recessive_neighbors",\
                "shortest_path_to_dominant_gene", "dominant_neighbors",\
                "is_single_coding_exon?",\
                "indel_position_in_CDS", "stop_position_in_CDS",\
                "causes_NMD?", "5'_flanking_splice_site",\
                "3'_flanking_splice_site", "canonical?",\
                "#_failed_filters", "filters_failed",\
                "ancestral_allele", "GERP_score", "GERP_element", "percentage_gerp_elements_in_truncated_exons", "truncated_exons:total_exons",\
                "segmental_duplications", "disorder_prediction"] + pfamParamsWithTruncations +\
                ["1000GPhase1", "1000GPhase1_AF", "1000GPhase1_ASN_AF",\
                "1000GPhase1_AFR_AF", "1000GPhase1_EUR_AF",\
                "ESP6500", "ESP6500_AAF",\
                "#_pseudogenes_associated_to_transcript",\
                "#_paralogs_associated_to_gene",\
                "dN/dS_(macaque)", "dN/dS_(mouse)"]
    spliceparams = ["shortest_path_to_recessive_gene", "recessive_neighbors",\
                "shortest_path_to_dominant_gene", "dominant_neighbors",\
                "donor", "acceptor",\
                "SNP_in_canonical_site?", "other_splice_site_canonical?",\
                "SNP_location", "alt_donor", "alt_acceptor", "nagnag_positions",\
                "intron_length", "#_failed_filters", "filters_failed",\
                "GERP_score", "GERP_element", "percentage_gerp_elements_in_truncated_exons", "truncated_exons:total_exons",\
                "segmental_duplications", "1000GPhase1", "1000GPhase1_AF", "1000GPhase1_ASN_AF",\
                "1000GPhase1_AFR_AF", "1000GPhase1_EUR_AF",\
                "ESP6500", "ESP6500_AAF",\
                "#_pseudogenes_associated_to_transcript",\
                "#_paralogs_associated_to_gene",\
                "dN/dS_(macaque)", "dN/dS_(mouse)"]
    outdata = {i : "" for i in set(basicparams) | set(LOFparams) | set(spliceparams)}

    lofOutputFile.write('chr\tpos\trsID\tref\talt\tscore\tPASS?\tdetails\t')
    lofOutputFile.write('\t'.join(i for i in basicparams)+'\t')
    lofOutputFile.write('\t'.join(i for i in LOFparams)+'\n')
    
    spliceOutputFile.write('chr\tpos\trsID\tref\talt\tscore\tPASS?\tdetails\t')
    spliceOutputFile.write('\t'.join(i for i in basicparams)+'\t')
    spliceOutputFile.write('\t'.join(i for i in spliceparams)+'\n')

    ##scan through VCF file metadata
    counter = 0
    vatFile.seek(0)
    line = vatFile.readline()
    while line=="\n" or line.startswith("#"):
        vcfOutputFile.write(line)
        counter+=1
        line = vatFile.readline()

    currentLoadedChromosome = None

    while line!="":
        data = line.strip().split('\t')
        chr_num = data[0].split("chr")[-1]
        start = int(data[1])
        end = start+len(data[3])-1

        if chr_num not in chrs:
            line = vatFile.readline()
            continue

        if not currentLoadedChromosome or currentLoadedChromosome != chr_num:
            #This is where we get a chance to data that is unique to a chromosome
            if VERBOSE: print("Reading data from chromosome %s..." % (chr_num))
            ancestorData = getAncestorData(args.ancestor, chr_num)
            exomesChromosomeInfo = getESP6500ExomeChromosomeInfo(args.exomes, chr_num) #Scan ESP6500 (exome) fields
            genomeSequences = getGenomeSequences(args.genome, chr_num)
            GERPelements = mergeElements(getGERPelements(open(os.path.join(args.elements, "hg19_chr%s_elems.txt" % (chr_num)))))
            currentLoadedChromosome = chr_num
        
        #Filter lines
        if "deletionFS" in line or "insertionFS" in line or "premature" in line or "splice" in line:
            ancesdata = ancestorData[start:start+len(data[3])].upper()
            if data[3] == ancesdata:
                ancestral = "Ref"
            elif data[4] == ancesdata:
                ancestral = "Alt"
            else:
                ancestral = "Neither"

            GERPscore = gerpScoresHash[data[0]][start]
            
            ##screen for variant types here.  skip variant if it is not deletion(N)FS, insertion(N)FS, or premature SNP
            lineinfo = {'AA':'AA='+ancesdata,\
                        'Ancestral':'Ancestral='+ancestral,\
                        'GERPscore':'GERPscore='+"%.2f" % GERPscore,\
                        'SegDup':'SegDup='+str(segdupdata[counter].count('('))}
            infotypes = ['AA', 'Ancestral', 'GERPscore', 'SegDup']
    
            outdata["ancestral_allele"] = ancesdata
            outdata["GERP_score"] = "%.2f" % GERPscore
            outdata["segmental_duplications"] = '.' if segdupdata[counter].count('(') == '0' else segdupdata[counter]
    
            #Adding 1000G fields
            thousandGTags = ['1000GPhase1_AF', '1000GPhase1_ASN_AF', '1000GPhase1_AFR_AF', '1000GPhase1_EUR_AF']
            thousandGComponents = []
            for thousandGTag in thousandGTags:
                thousandGComponents.append(thousandGTag + "=NA")
            
            if chr_num in thousandGChromosomeInfo and start in thousandGChromosomeInfo[chr_num]:
                for info in thousandGChromosomeInfo[chr_num][start].split(";"):
                    infotype = info.split('=')[0]  
                    newComponent = "1000GPhase1_" + info
                    thousandGComponentIndex = -1
                    for findIndex in range(len(thousandGTags)):
                        if infotype == "_".join(thousandGTags[findIndex].split("_")[1:]):
                            thousandGComponentIndex = findIndex
                            break
                    
                    if thousandGComponentIndex >= 0:
                        thousandGComponents[thousandGComponentIndex] = newComponent
            
            infotypes += ['1000GPhase1'] + thousandGTags
            if chr_num in thousandGChromosomeInfo and start in thousandGChromosomeInfo[chr_num]:
                lineinfo['1000GPhase1'] = '1000GPhase1=Yes'
            else:
                lineinfo['1000GPhase1'] = '1000GPhase1=No'
            
            #Add 1000G entries to output
            for tagIndex in range(len(thousandGTags)):
                lineinfo[thousandGTags[tagIndex]] = thousandGComponents[tagIndex]
            
            #Add exomes info to output
            infotypes += ['ESP6500', 'ESP6500_AAF']
            if start in exomesChromosomeInfo:
                lineinfo['ESP6500'] = 'ESP6500=Yes'
                lineinfo['ESP6500_AAF'] = 'ESP6500_AAF=' + exomesChromosomeInfo[start]
            else:
                lineinfo['ESP6500'] = 'ESP6500=No'
                lineinfo['ESP6500_AAF'] = 'ESP6500_AAF=NA,NA,NA'
    
            for tag in ['1000GPhase1'] + thousandGTags + ['ESP6500', 'ESP6500_AAF']:
                outdata[tag] = lineinfo[tag]
            
            dataInfoComponents = data[7].split(';')
            found = 0
            for info in dataInfoComponents:
                infotype = info.split('=')[0]
                if infotype == 'VA':
                    variants = info.split('VA=')[-1].split(',')
                    found = 1
                if infotype!='AA' and infotype!='VA':
                    lineinfo[infotype]=info
                    infotypes.append(infotype)
            
            if found==1:
                lineinfo['VA']='VA='
                infotypes.append('VA')

            def insertAncestralField(outfile):
                ancestralInsertion = ";".join([lineinfo['Ancestral']] + data[7].split(";"))
                outfile.write("\t".join(data[0:7] + [ancestralInsertion] + data[8:]))
                outfile.write('\t'+ '\t'.join(outdata[i] for i in basicparams))
            
            LOFvariants = []
            splicevariants = []
            othervariants = []
            for variant in variants:
                ##alternate allele corresponding to variant
                subst = data[4].split(',')[int(variant.split(':')[0])-1]
                
                if "deletionFS" not in variant and "insertionFS" not in variant:
                    if "premature" not in variant and "splice" not in variant:
                        othervariants.append(variant)
                        continue
                details = variant.split(":")
    
                outdata["gene"], outdata["gene_id"] = details[1], details[2]
    
                if details[5].split("/")[0]==details[5].split("/")[1]:
                    pf = "full"
                else:
                    pf = "partial"
                outdata["partial/full"] = pf
                
                transcripts = []
    
                for i in range(6, len(details)-1, 3):
                    transcripts.append(details[i:i+3])
                longesttranscript = max([int(i[2].split('_')[0]) for i in transcripts])
    
                ##calculate distance to dominant and recessive genes
                gene_name = outdata["gene"]
                if ppi is not None and gene_name in ppi:
                    dominantdist, numberOfDominantNeighbors = parsePPI(networkx, ppi, ppiHash, "dgenes", gene_name, dgenes)
                    outdata["shortest_path_to_dominant_gene"] = 'NA' if dominantdist is None else str(dominantdist)
                    outdata["dominant_neighbors"] = str(numberOfDominantNeighbors)

                    recessdist, numberOfRecessiveNeighbors = parsePPI(networkx, ppi, ppiHash, "rgenes", gene_name, rgenes)
                    outdata["shortest_path_to_recessive_gene"] = 'NA' if recessdist is None else str(recessdist)
                    outdata["recessive_neighbors"] = str(numberOfRecessiveNeighbors)
                else:
                    outdata["shortest_path_to_recessive_gene"] = 'NA'
                    outdata["recessive_neighbors"] = 'NA'
    
                    outdata["shortest_path_to_dominant_gene"] = 'NA'
                    outdata["dominant_neighbors"] = 'NA'
    
                outdata["#_paralogs_associated_to_gene"] = str(len(paralogs[outdata["gene_id"].split('.')[0]])) if outdata["gene_id"].split('.')[0] in paralogs else "0"
    
                ##number of associated pseudogenes computation goes here
    
                if "splice" in variant:
                    ##check that is a SNP splice variant
                    if len(data[3])>1 or len(subst)>1:
                        splicevariants.append(variant)
                        continue
                    splicevariants.append(':'.join(details[:6]))
    
                    for entry in transcripts:
                        splicevariants[-1]+=':' + ':'.join(entry[0:1] + [pf] + entry[1:])
                        transcript = entry[1]
                        outdata["transcript"] = transcript
                        outdata["transcript_length"] = entry[2]
                        outdata["longest_transcript?"] = "YES" if int(outdata["transcript_length"])==longesttranscript else "NO"
                        ispositivestr = transcript_strand[transcript]=='+'

                        GERPelementdata, GERPrejectiondata, exonCountData = getGERPData(True, GERPelements, codingExonIntervals[chr_num][transcript] if transcript in codingExonIntervals[chr_num] else None, start, end, transcript_strand[transcript])

                        outdata['GERP_element'] = GERPelementdata
                        outdata['percentage_gerp_elements_in_truncated_exons'] = GERPrejectiondata
                        outdata['truncated_exons:total_exons'] = exonCountData
    
                        outdata["#_pseudogenes_associated_to_transcript"] = str(numpseudogenes[transcript]) if transcript in numpseudogenes else "0"

                        macaque = 'NA'
                        if transcript.split('.')[0] in dNdSmacaque:
                            if dNdSmacaque[transcript.split('.')[0]] != 'N/A':
                                macaque = "%.3f" % float(dNdSmacaque[transcript.split('.')[0]])

                        outdata["dN/dS_(macaque)"] =  macaque

                        mouse = 'NA'
                        if transcript.split('.')[0] in dNdSmouse:
                            if dNdSmouse[transcript.split('.')[0]] != 'N/A':
                                mouse = "%.3f" % float(dNdSmouse[transcript.split('.')[0]])

                        outdata["dN/dS_(mouse)"] = mouse
                        
                        insertAncestralField(spliceOutputFile)

                        spliceSearchData = searchInSplices(chr_num, transcript, genomeSequences, ispositivestr, start, CDS, subst)

                        def writeSpliceOutput(failure):
                            spliceOutputFile.write('\t'+'\t'.join(outdata[i] for i in ["shortest_path_to_recessive_gene", "recessive_neighbors"]))
                            spliceOutputFile.write("\t%s: pos=" % (failure) +str(start)+' transcript='+transcript+'\n')

                        if not spliceSearchData['found']:
                            writeSpliceOutput("CDS match not found")
                            continue

                        if not spliceSearchData['new']:
                            writeSpliceOutput("no donor/acceptor pair")
                            continue

                        new = spliceSearchData['new']
                        donor = spliceSearchData['donor']
                        acceptor = spliceSearchData['acceptor']
                        intronlength = spliceSearchData['intronlength']

                        outdata["donor"] = donor
                        outdata["acceptor"] = acceptor
                        outdata["intron_length"] = str(intronlength)
                        ##write to output
                        if new[0]==0:
                            isCanonical = 'YES' if (donor in ['GT', 'GC'] or (donor == 'AT' and acceptor == 'AC')) else 'NO'
                            otherCanonical = 'YES' if (acceptor=='AG' or (acceptor == 'AC' and donor == 'AT')) else 'NO'
                        elif new[0]==1:
                            isCanonical = 'YES' if (acceptor=='AG' or (acceptor == 'AC' and donor == 'AT')) else 'NO'
                            otherCanonical = 'YES' if (donor in ['GT', 'GC'] or (donor == 'AT' and acceptor == 'AC')) else 'NO'
                        outdata["SNP_in_canonical_site?"] = isCanonical
                        outdata["other_splice_site_canonical?"] = otherCanonical
                        
                        if new[0]==0:
                            outdata["SNP_location"] = "donor"
                            outdata["alt_donor"] = new[1].upper()
                            outdata["alt_acceptor"] = acceptor
                        else:
                            outdata["SNP_location"] = "acceptor"
                            outdata["alt_donor"] = donor
                            outdata["alt_acceptor"] = new[1].upper()

                        if new[0] == 1: #acceptor snp location
                            nagNagPositions = getMatchingNagnagnagPositions(genomeSequences, start, ispositivestr)
                            outdata['nagnag_positions'] = '/'.join(map(str, nagNagPositions)) if len(nagNagPositions) > 0 else '.'
                            alternateAcceptorSite = 'YES' if len(nagNagPositions) > 0 else 'NO'
                        else:
                            outdata['nagnag_positions'] = 'NA'
                            alternateAcceptorSite = 'NA'
    
                        #calculation of filters
                        failed_filters = []
                        if isCanonical == 'NO':
                            failed_filters.append('ref_noncanonical')
                            if new[0] == 0 and new[1].upper() != 'GT': #snp is in donor, and alt donor is not GT
                                failed_filters.append('alt_noncanonical')
                            elif new[0] != 0 and new[1].upper() != 'AG': #snp is in acceptor, and alt acceptor is not AG
                                failed_filters.append('alt_noncanonical')

                        if otherCanonical == 'NO':
                            failed_filters.append('other_noncanonical')
                        if intronlength < 15:
                            failed_filters.append('short_intron')
                            smallIntron = 'YES'
                        else:
                            smallIntron = 'NO'
                        if segdupdata[counter].count('(') > 3:
                            failed_filters.append('heavily_duplicated')
                            heavilyDuplicated = 'YES'
                        else:
                            heavilyDuplicated = 'NO'

                        isLofAnc = 'NO'
                        if ancesdata==subst:
                            failed_filters.append('lof_anc')
                            isLofAnc = 'YES'
    
                        outdata["#_failed_filters"] = str(len(failed_filters))
                        outdata["filters_failed"] = ','.join(failed_filters)

    ########################################################
                        spliceOutputFile.write("\t"+"\t".join(outdata[i] for i in spliceparams)+"\n")
    #########################################################
                        splicevariants[-1]+=':'+':'.join(['GERPelement='+("YES" if GERPelementdata != '.' else "NO"), 'exoncounts='+exonCountData, donor+'/'+acceptor, 'is_canonical=' + isCanonical, 'other_noncanonical=' + otherCanonical, 'intron_length=' + str(intronlength), 'small_intron=' + smallIntron, 'heavily_duplicated=' + heavilyDuplicated, 'lof_anc=' + isLofAnc, 'alternate_acceptor_site=' + alternateAcceptorSite])
                        
                else:   ##deletionFS, insertionFS, or prematureStop
                    LOFvariants.append(':'.join(details[:6]))
    
                    for entry in transcripts:
                        LOFvariants[-1]+=':'+':'.join(entry[0:1] + [pf] + entry[1:])
                        
                        tlength = entry[2].split('_')[0]
                        outdata["transcript_length"] = tlength
                        try:
                            LOFposition = entry[2].split('_')[1]
                        except:
                            LOFposition = '.'
                        outdata["longest_transcript?"] = "YES" if int(tlength)==longesttranscript else "NO"
                        transcript = entry[1]
                        outdata["transcript"]=transcript
                       
    		  #calculation of filters
                        filters_failed = 0
                        failed_filters = []

                        nearStart = 'NO'
                        nearEnd = 'NO'
                        try:	#since LOFposition may not be provided
                            if float(LOFposition)/float(tlength) <= 0.05:
                                filters_failed = filters_failed+1
                                failed_filters.append('near_start')
                                nearStart = 'YES'
                            if float(LOFposition)/float(tlength) >= 0.95:
                                filters_failed = filters_failed+1
                                failed_filters.append('near_stop')
                                nearEnd = 'YES'
                        except:
                            pass
                        
                        isLofAnc = 'NO'
                        if ancesdata==subst:
                            filters_failed = filters_failed+1
                            failed_filters.append('lof_anc')
                            isLofAnc = 'YES'
                        
                        heavilyDuplicated = 'NO'
                        if segdupdata[counter].count('(') > 3:
                            filters_failed = filters_failed+1
                            failed_filters.append('heavily_duplicated')
                            heavilyDuplicated = 'YES'

                        outdata["#_failed_filters"] = str(filters_failed)
                        outdata["filters_failed"] = ','.join(failed_filters)
    
                        outdata["indel_position_in_CDS"] = "NA"
                        outdata["stop_position_in_CDS"] = "NA"
                        outdata["5'_flanking_splice_site"] = "NA"
                        outdata["3'_flanking_splice _site"] = "NA"
                        outdata["canonical?"] = "NA"
                        outdata["#_pseudogenes_associated_to_transcript"] = str(numpseudogenes[transcript]) if transcript in numpseudogenes else "0"
                        outdata["dN/dS_(macaque)"] = dNdSmacaque[transcript.split('.')[0]] if transcript.split('.')[0] in dNdSmacaque else "NA"
                        outdata["dN/dS_(mouse)"] = dNdSmouse[transcript.split('.')[0]] if transcript.split('.')[0] in dNdSmouse else "NA"
                        
                        insertAncestralField(lofOutputFile)
                        
                        nmdData = findNMDForIndelsAndPrematureStop(args.nmd_threshold, data, chr_num, transcript, start, end, exon, stop_codon, genomeSequences, CDS, subst, transcript_strand)

                        if nmdData['NMD'] is None:
                            continue

                        outdata['causes_NMD?'] = nmdData['NMD']

                        if nmdData['issinglecodingexon']:
                            outdata["is_single_coding_exon?"] = nmdData['issinglecodingexon']

                        #NA for premature SNPs
                        if nmdData['newCDSpos'] and not ('prematureStop' in variant and (len(data[3])>1 or len(subst)>1)):
                            outdata['indel_position_in_CDS'] = str(nmdData['newCDSpos'])

                        if nmdData['NMD'] not in ['YES', 'NO']:
                            lofOutputFile.write('\t'+'\t'.join(outdata[i] for i in LOFparams) + '\n')
                            continue

                        lofPosition = nmdData['newCDSpos'] if "prematureStop" in variant else nmdData['stopCDS']

                        outdata["5'_flanking_splice_site"] = nmdData['splice1']
                        outdata["3'_flanking_splice_site"] = nmdData['splice2']
                        outdata["canonical?"] = nmdData['canonical']
                        outdata["stop_position_in_CDS"] = str(lofPosition)

                        vcfPfamDescriptions = {}
                        phosphorylationResults = {}

                        oneNA = False
                        oneNO = False
                        oneYES = False

                        stopPositionInAminoSpace = int(entry[2].split('_')[2]) if "prematureStop" in variant else (lofPosition - 1) // 3 + 1

                        if "prematureStop" not in variant:
                            stopPositionForGERP = calculateAbsolutePosition(lofPosition, codingExonIntervals[chr_num][transcript], transcript_strand[transcript])
                        else:
                            stopPositionForGERP = start
                        
                        GERPelementdata, GERPrejectiondata, exonCountData = getGERPData(False, GERPelements, codingExonIntervals[chr_num][transcript] if transcript in codingExonIntervals[chr_num] else None, stopPositionForGERP, stopPositionForGERP + len(data[3]) - 1, transcript_strand[transcript])

                        outdata['GERP_element'] = GERPelementdata
                        outdata['percentage_gerp_elements_in_truncated_exons'] = GERPrejectiondata
                        outdata['truncated_exons:total_exons'] = exonCountData

                        for paramKey in pfamParams + ptmParams:
                            newDescriptions = getPfamDescription(transcriptToProteinHash, chr_num, transcript.split(".")[0], stopPositionInAminoSpace, chromosomesPFam, paramKey)
                            #we just want a select few domain types for VCF output
                            if paramKey in pfamParams:
                                vcfPfamDescriptions[paramKey] = ":%s=%s" % (paramKey, newDescriptions[1])
                            elif paramKey in phosphorylationTags:
                                if newDescriptions[1] == "YES":
                                    phosphorylationResults[paramKey] = newDescriptions[1]

                                if newDescriptions[1] == 'YES':
                                    oneYES = True
                                elif newDescriptions[1] == 'NO':
                                    oneNO = True
                                elif newDescriptions[1] == 'NA':
                                    oneNA = True
                            else:
                                vcfPfamDescriptions[paramKey] = ''

                            outdata[paramKey] = newDescriptions[2][0]
                            outdata[pfamParamsWithTruncations[pfamParamsWithTruncations.index(paramKey)+1]] = newDescriptions[2][1]

                        if len(phosphorylationResults) == 0:
                            if oneNA and oneNO:
                                vcfPfamDescriptions[phosphorylationTags[0]] = ':PTM=NO/NA'
                            elif oneNA:
                                vcfPfamDescriptions[phosphorylationTags[0]] = ':PTM=NA'
                            else:
                                vcfPfamDescriptions[phosphorylationTags[0]] = ':PTM=NO'
                            #outdata["PTM"] = "PTM=NO"
                        else:
                            vcfPfamDescriptions[phosphorylationTags[0]] = ':PTM=' + '|'.join([key + "/" + value for key, value in phosphorylationResults.items()])
                            #outdata["PTM"] = "PTM=" + ','.join([key + "/" + value for key, value in phosphorylationResults.items()])

                        disorderPredictionData = getDisopredData(args.disopred_sequences, transcript, stopPositionInAminoSpace)
                        outdata["disorder_prediction"] = disorderPredictionData

    #########################################################
                        lofOutputFile.write('\t' + '\t'.join(outdata[i] for i in LOFparams)+'\n')
    #########################################################
                        LOFvariants[-1]+=':'+':'.join(['GERPelement='+("YES" if GERPelementdata != '.' else "NO"), 'exoncounts='+exonCountData, 'nearstart=' + nearStart, 'nearend=' + nearEnd, 'canonical='+nmdData['canonical'], nmdData['splice1']+'/'+nmdData['splice2'], str(nmdData['newCDSpos']), 'lofposition='+str(lofPosition), nmdData['nextATG'], 'nmd=' + nmdData['NMD'], nmdData['incrcodingpos'], 'lof_anc=' + isLofAnc, 'heavilyduplicated='+heavilyDuplicated, 'disorder_prediction='+disorderPredictionData]) + ''.join([vcfPfamDescriptions[param] for param in pfamParams + [phosphorylationTags[0]]])

            vcfOutputFile.write('\t'.join(data[k] for k in range(0,7))+'\t')
            allvariants = []
            for variant in LOFvariants:
                allvariants.append(variant)
            for variant in splicevariants:
                allvariants.append(variant)
            for variant in othervariants:
                allvariants.append(variant)
            lineinfo['VA']+=','.join(allvariants) 
            vcfOutputFile.write(';'.join(lineinfo[infotype] for infotype in infotypes))
            if len(data[8:]) > 0:
                vcfOutputFile.write('\t' + '\t'.join(data[8:]) + '\n')
            else:
                vcfOutputFile.write('\n')
        
        line=vatFile.readline()
        counter+=1
    
    vcfOutputFile.close()
    lofOutputFile.close()
    spliceOutputFile.close()
    vatFile.close()

    if ppiHash is not None:
        try:
            #save shortest path values to cache file
            pickle.dump(ppiHash, open(ppiHashPath, "wb"), protocol=2)
        except:
            printError("Failed to write PPI cache, skipping..", False)

    if VERBOSE: print("Finished execution in %d seconds" % ((datetime.datetime.now() - startProgramExecutionTime).seconds))

if __name__ == "__main__":
    main()
