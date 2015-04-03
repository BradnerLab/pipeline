#!/usr/bin/env python
'''
PROGRAM TO STITCH TOGETHER REGIONS TO FORM ENHANCERS, MAP READ DENSITY TO STITCHED REGIONS,
AND RANK ENHANCERS BY READ DENSITY TO DISCOVER SUPER-ENHANCERS
May2, 2014
VERSION 0.2
CONTACT: youngcomputation@wi.mit.edu
'''

from __future__ import absolute_import  # , division, print_function, unicode_literals
import utils

import sys
# import ROSE_utils
import time
import copy
import os
import numpy
import subprocess

from collections import defaultdict




#==================================================================
#=====================HELPER FUNCTIONS=============================
#==================================================================

def getBamChromList(bamFileList):

    '''
    gets the consensus list of chromosomes mapped by the bams
    '''
    
    #start w/ the first bam
    cmd = 'samtools idxstats %s' % (bamFileList[0])
    idxStats = subprocess.Popen(cmd,stdout=subprocess.PIPE,shell=True)
    idxStats= idxStats.communicate()
    finalChromList = [line.split('\t')[0] for line in idxStats[0].split('\n')[0:-2]]
    
    #now go through each additional bam
    for bamFile in bamFileList:
        cmd = 'samtools idxstats %s' % (bamFile)
        idxStats = subprocess.Popen(cmd,stdout=subprocess.PIPE,shell=True)
        idxStats= idxStats.communicate()
        chromList = [line.split('\t')[0] for line in idxStats[0].split('\n')[0:-2]]
        finalChromList = [chrom for chrom in finalChromList if chromList.count(chrom) != 0]

    return utils.uniquify(finalChromList)


def checkRefCollection(referenceCollection):

    '''
    makes sure the names of all loci in the reference collection are unique
    '''

    namesList = [locus.ID() for locus in referenceCollection.getLoci()]
    
    if len(namesList) != len(utils.uniquify(namesList)):
        print("ERROR: REGIONS HAVE NON-UNIQUE IDENTIFIERS")
        sys.exit()
    else:
        print("REFERENCE COLLECTION PASSES QC")
        return

    

def filterGFF(gffFile,chromList):

    '''
    takes in a gff and filters out all lines that don't belong to a chrom in the chromList
    '''
    gff = utils.parseTable(gffFile,'\t')
    filteredGFF = []
    excludeList=[]
    for line in gff:
        if chromList.count(line[0]) ==1:
            filteredGFF.append(line)
        else:
            excludeList.append(line[0])

    excludeList = utils.uniquify(excludeList)
    if len(excludeList) > 0:
        print("EXCLUDED GFF REGIONS FROM THE FALLING CHROMS: %s" % (','.join(excludeList)))

    return filteredGFF
             

          




#==================================================================
#=====================REGION STITCHING=============================
#==================================================================




def optimizeStitching(locusCollection, name, outFolder, stepSize=500):
    '''
    takes a locus collection and starts writing out stitching stats at step sized intervals
    '''
    maxStitch = 15000  # set a hard wired match stitching parameter

    stitchTable = [['STEP', 'NUM_REGIONS', 'TOTAL_CONSTIT', 'TOTAL_REGION', 'MEAN_CONSTIT', 'MEDIAN_CONSTIT', 'MEAN_REGION', 'MEDIAN_REGION', 'MEAN_STITCH_FRACTION', 'MEDIAN_STITCH_FRACTION']]
    # first consolidate the collection
    locusCollection = locusCollection.stitchCollection(stitchWindow=0)
    total_constit = sum([locus.len() for locus in locusCollection.getLoci()])
    step = 0
    while step <= maxStitch:

        print("Getting stitch stats for %s (bp)" % (step))
        stitchCollection = locusCollection.stitchCollection(stitchWindow=step)
        num_regions = len(stitchCollection)
        stitchLoci = stitchCollection.getLoci()
        regionLengths = [locus.len() for locus in stitchLoci]
        total_region = sum(regionLengths)
        constitLengths = []
        for locus in stitchLoci:

            constitLoci = locusCollection.getOverlap(locus)
            constitLengths.append(sum([locus.len() for locus in constitLoci]))

        meanConstit = round(numpy.mean(constitLengths), 2)
        medianConstit = round(numpy.median(constitLengths), 2)

        meanRegion = round(numpy.mean(regionLengths), 2)
        medianRegion = round(numpy.median(regionLengths), 2)

        stitchFractions = [float(constitLengths[i]) / float(regionLengths[i]) for i in range(len(regionLengths))]
        meanStitchFraction = round(numpy.mean(stitchFractions), 2)
        medianStitchFraction = round(numpy.median(stitchFractions), 2)

        newLine = [step, num_regions, total_constit, total_region, meanConstit, medianConstit, meanRegion, medianRegion, meanStitchFraction, medianStitchFraction]

        stitchTable.append(newLine)

        step += stepSize

    # write the stitch table to disk
    stitchParamFile = '%s%s_stitch_params.tmp' % (outFolder, name)
    utils.unParseTable(stitchTable, stitchParamFile, '\t')
    # call the rscript
    rCmd = 'Rscript ./ROSE2_stitchOpt.R %s %s %s' % (stitchParamFile, outFolder, name)
    print(rCmd)
    # get back the stitch parameter
    rOutput = subprocess.Popen(rCmd, stdout=subprocess.PIPE, shell=True)
    rOutputTest = rOutput.communicate()

    print(rOutputTest)

    stitchParam = rOutputTest[0].split('\n')[2]
    try:
        stitchParam = int(stitchParam)
    except ValueError:
        print("INVALID STITCHING PARAMETER. STITCHING OPTIMIZATION FAILED")
        sys.exit()

    # delete? the table
    # os.system('rm -f %s' % (stitchParamFile))
    return stitchParam


def regionStitching(referenceCollection, name, outFolder, stitchWindow, tssWindow, annotFile, removeTSS=True):
    print('PERFORMING REGION STITCHING')
    # first have to turn bound region file into a locus collection

    # need to make sure this names correctly... each region should have a unique name
    #referenceCollection 

    debugOutput = []
    # filter out all bound regions that overlap the TSS of an ACTIVE GENE
    if removeTSS:

        print('REMOVING TSS FROM REGIONS USING AN EXCLUSION WINDOW OF %sBP' % (tssWindow))
        # first make a locus collection of TSS

        startDict = utils.makeStartDict(annotFile)

        # now makeTSS loci for active genes
        removeTicker = 0
        # this loop makes a locus centered around +/- tssWindow of transcribed genes
        # then adds it to the list tssLoci
        tssLoci = []
        for geneID in startDict.keys():
            tssLoci.append(utils.makeTSSLocus(geneID, startDict, tssWindow, tssWindow))

        # this turns the tssLoci list into a LocusCollection
        # 50 is the internal parameter for LocusCollection and doesn't really matter
        tssCollection = utils.LocusCollection(tssLoci, 50)

        # gives all the loci in referenceCollection
        boundLoci = referenceCollection.getLoci()

        # this loop will check if each bound region is contained by the TSS exclusion zone
        # this will drop out a lot of the promoter only regions that are tiny
        # typical exclusion window is around 2kb
        for locus in boundLoci:
            if len(tssCollection.getContainers(locus, 'both')) > 0:

                # if true, the bound locus overlaps an active gene
                referenceCollection.remove(locus)
                debugOutput.append([locus.__str__(), locus.ID(), 'CONTAINED'])
                removeTicker += 1
        print('REMOVED %s LOCI BECAUSE THEY WERE CONTAINED BY A TSS' % (removeTicker))

    # referenceCollection is now all enriched region loci that don't overlap an active TSS

    if stitchWindow == '':
        print('DETERMINING OPTIMUM STITCHING PARAMTER')
        optCollection = copy.deepcopy(referenceCollection)
        stitchWindow = optimizeStitching(optCollection, name, outFolder, stepSize=500)
    print('USING A STITCHING PARAMETER OF %s' % stitchWindow)
    stitchedCollection = referenceCollection.stitchCollection(stitchWindow, 'both')

    if removeTSS:
        # now replace any stitched region that overlap 2 distinct genes
        # with the original loci that were there
        fixedLoci = []
        tssLoci = []
        for geneID in startDict.keys():
            tssLoci.append(utils.makeTSSLocus(geneID, startDict, 50, 50))

        # this turns the tssLoci list into a LocusCollection
        # 50 is the internal parameter for LocusCollection and doesn't really matter
        tssCollection = utils.LocusCollection(tssLoci, 50)
        removeTicker = 0
        originalTicker = 0
        for stitchedLocus in stitchedCollection.getLoci():
            overlappingTSSLoci = tssCollection.getOverlap(stitchedLocus, 'both')
            tssNames = [startDict[tssLocus.ID()]['name'] for tssLocus in overlappingTSSLoci]
            tssNames = utils.uniquify(tssNames)
            if len(tssNames) > 2:

                # stitchedCollection.remove(stitchedLocus)
                originalLoci = referenceCollection.getOverlap(stitchedLocus, 'both')
                originalTicker += len(originalLoci)
                fixedLoci += originalLoci
                debugOutput.append([stitchedLocus.__str__(), stitchedLocus.ID(), 'MULTIPLE_TSS'])
                removeTicker += 1
            else:
                fixedLoci.append(stitchedLocus)

        print('REMOVED %s STITCHED LOCI BECAUSE THEY OVERLAPPED MULTIPLE TSSs' % (removeTicker))
        print('ADDED BACK %s ORIGINAL LOCI' % (originalTicker))
        fixedCollection = utils.LocusCollection(fixedLoci, 50)
        return fixedCollection, debugOutput, stitchWindow
    else:
        return stitchedCollection, debugOutput, stitchWindow

#==================================================================
#=====================REGION LINKING MAPPING=======================
#==================================================================


def mapCollection(stitchedCollection, referenceCollection, bamFileList, mappedFolder, output, refName):
    '''
    makes a table of factor density in a stitched locus and ranks table by number of loci stitched together
    '''

    print('FORMATTING TABLE')
    loci = stitchedCollection.getLoci()

    locusTable = [['REGION_ID', 'CHROM', 'START', 'STOP', 'NUM_LOCI', 'CONSTITUENT_SIZE']]

    lociLenList = []

    # strip out any that are in chrY
    for locus in list(loci):
        if locus.chr() == 'chrY':
            loci.remove(locus)

    for locus in loci:
        # numLociList.append(int(stitchLocus.ID().split('_')[1]))
        lociLenList.append(locus.len())
        # numOrder = order(numLociList,decreasing=True)
    lenOrder = utils.order(lociLenList, decreasing=True)
    ticker = 0
    for i in lenOrder:
        ticker += 1
        if ticker % 1000 == 0:
            print(ticker)
        locus = loci[i]

        # First get the size of the enriched regions within the stitched locus
        refEnrichSize = 0
        refOverlappingLoci = referenceCollection.getOverlap(locus, 'both')
        for refLocus in refOverlappingLoci:
            refEnrichSize += refLocus.len()

        try:
            stitchCount = int(locus.ID().split('_')[0])
        except ValueError:
            stitchCount = 1
        coords = [int(x) for x in locus.coords()]

        locusTable.append([locus.ID(), locus.chr(), min(coords), max(coords), stitchCount, refEnrichSize])

    print('GETTING MAPPED DATA')
    print("USING A BAMFILE LIST:")
    print(bamFileList)
    for bamFile in bamFileList:

        bamFileName = bamFile.split('/')[-1]

        print('GETTING MAPPING DATA FOR  %s' % bamFile)
        # assumes standard convention for naming enriched region gffs

        # opening up the mapped GFF
        print('OPENING %s%s_%s_MAPPED/matrix.txt' % (mappedFolder, refName, bamFileName))

        mappedGFF = utils.parseTable('%s%s_%s_MAPPED/matrix.txt' % (mappedFolder, refName, bamFileName), '\t')

        signalDict = defaultdict(float)
        print('MAKING SIGNAL DICT FOR %s' % (bamFile))
        mappedLoci = []
        for line in mappedGFF[1:]:

            chrom = line[1].split('(')[0]
            start = int(line[1].split(':')[-1].split('-')[0])
            end = int(line[1].split(':')[-1].split('-')[1])
            mappedLoci.append(utils.Locus(chrom, start, end, '.', line[0]))
            try:
                signalDict[line[0]] = float(line[2]) * (abs(end - start))
            except ValueError:
                print('WARNING NO SIGNAL FOR LINE:')
                print(line)
                continue

        mappedCollection = utils.LocusCollection(mappedLoci, 500)
        locusTable[0].append(bamFileName)

        for i in range(1, len(locusTable)):
            signal = 0.0
            line = locusTable[i]
            lineLocus = utils.Locus(line[1], line[2], line[3], '.')
            overlappingRegions = mappedCollection.getOverlap(lineLocus, sense='both')
            for region in overlappingRegions:
                signal += signalDict[region.ID()]
            locusTable[i].append(signal)

    utils.unParseTable(locusTable, output, '\t')


#==================================================================
#=========================MAIN METHOD==============================
#==================================================================
def main():
    '''
    main run call
    '''
    debug = False

    from optparse import OptionParser
    usage = "usage: %prog [options] -g [GENOME] -i [INPUT_REGION_GFF] -r [RANKBY_BAM_FILE] -o [OUTPUT_FOLDER] [OPTIONAL_FLAGS]"
    parser = OptionParser(usage=usage)
    # required flags
    parser.add_option("-i", "--i", dest="input", nargs=1, default=None,
                      help="Enter a .gff or .bed file of binding sites used to make enhancers")
    parser.add_option("-r", "--rankby", dest="rankby", nargs=1, default=None,
                      help="bamfile to rank enhancer by")
    parser.add_option("-o", "--out", dest="out", nargs=1, default=None,
                      help="Enter an output folder")
    parser.add_option("-g", "--genome", dest="genome", nargs=1, default=None,
                      help="Enter the genome build (MM9,MM8,HG18,HG19)")

    # optional flags
    parser.add_option("-b", "--bams", dest="bams", nargs=1, default=None,
                      help="Enter a comma separated list of additional bam files to map to")
    parser.add_option("-c", "--control", dest="control", nargs=1, default=None,
                      help="bamfile to rank enhancer by")
    parser.add_option("-s", "--stitch", dest="stitch", nargs=1, default='',
                      help="Enter a max linking distance for stitching. Default will determine optimal stitching parameter")
    parser.add_option("-t", "--tss", dest="tss", nargs=1, default=0,
                      help="Enter a distance from TSS to exclude. 0 = no TSS exclusion")

    parser.add_option("--mask", dest="mask", nargs=1, default=None,
                      help="Mask a set of regions from analysis.  Provide a .bed or .gff of masking regions")

    # RETRIEVING FLAGS
    (options, args) = parser.parse_args()

    if not options.input or not options.rankby or not options.out or not options.genome:
        print('hi there')
        parser.print_help()
        exit()

    # making the out folder if it doesn't exist
    outFolder = utils.formatFolder(options.out, True)

    # figuring out folder schema
    gffFolder = utils.formatFolder(outFolder + 'gff/', True)
    mappedFolder = utils.formatFolder(outFolder + 'mappedGFF/', True)

    # GETTING INPUT FILE
    if options.input.split('.')[-1] == 'bed':
        # CONVERTING A BED TO GFF
        inputGFFName = options.input.split('/')[-1][0:-4]
        inputGFFFile = '%s%s.gff' % (gffFolder, inputGFFName)
        utils.bedToGFF(options.input, inputGFFFile)
    elif options.input.split('.')[-1] == 'gff':
        # COPY THE INPUT GFF TO THE GFF FOLDER
        inputGFFFile = options.input
        os.system('cp %s %s' % (inputGFFFile, gffFolder))

    else:
        print('WARNING: INPUT FILE DOES NOT END IN .gff or .bed. ASSUMING .gff FILE FORMAT')
        # COPY THE INPUT GFF TO THE GFF FOLDER
        inputGFFFile = options.input
        os.system('cp %s %s' % (inputGFFFile, gffFolder))

    # GETTING THE LIST OF BAMFILES TO PROCESS
    if options.control:
        bamFileList = [options.rankby, options.control]

    else:
        bamFileList = [options.rankby]

    if options.bams:
        bamFileList += options.bams.split(',')
        #bamFileList = utils.uniquify(bamFileList) # makes sad when you have the same control bam over and over again
    # optional args

    # Stitch parameter
    if options.stitch == '':
        stitchWindow = ''
    else:
        stitchWindow = int(options.stitch)

    # tss options
    tssWindow = int(options.tss)
    if tssWindow != 0:
        removeTSS = True
    else:
        removeTSS = False

    # GETTING THE BOUND REGION FILE USED TO DEFINE ENHANCERS
    print('USING %s AS THE INPUT GFF' % (inputGFFFile))
    inputName = inputGFFFile.split('/')[-1].split('.')[0]

    # GETTING THE GENOME
    genome = options.genome
    print('USING %s AS THE GENOME' % genome)

    # GETTING THE CORRECT ANNOT FILE
    cwd = os.getcwd()
    genomeDict = {
        'HG18': '%s/annotation/hg18_refseq.ucsc' % (cwd),
        'MM9': '%s/annotation/mm9_refseq.ucsc' % (cwd),
        'HG19': '%s/annotation/hg19_refseq.ucsc' % (cwd),
        'MM8': '%s/annotation/mm8_refseq.ucsc' % (cwd),
        'MM10': '%s/annotation/mm10_refseq.ucsc' % (cwd),
        'RN4': '%s/annotation/rn4_refseq.ucsc' % (cwd),
    }

    annotFile = genomeDict[genome.upper()]

    # MAKING THE START DICT
    print('MAKING START DICT')
    startDict = utils.makeStartDict(annotFile)


    #GET CHROMS FOUND IN THE BAMS
    print('GETTING CHROMS IN BAMFILES')
    bamChromList = getBamChromList(bamFileList)
    print("USING THE FOLLOWING CHROMS")
    print(bamChromList)

    #LOADING IN THE GFF AND FILTERING BY CHROM
    print('LOADING AND FILTERING THE GFF')
    inputGFF = filterGFF(inputGFFFile,bamChromList)
    # LOADING IN THE BOUND REGION REFERENCE COLLECTION
    print('LOADING IN GFF REGIONS')
    referenceCollection = utils.gffToLocusCollection(inputGFF)

    print('CHECKING REFERENCE COLLECTION:')
    checkRefCollection(referenceCollection)
        

    # MASKING REFERENCE COLLECTION
    # see if there's a mask
    if options.mask:
        maskFile = options.mask
        # if it's a bed file
        if maskFile.split('.')[-1].upper() == 'BED':
            maskGFF = utils.bedToGFF(maskFile)
        elif maskFile.split('.')[-1].upper() == 'GFF':
            maskGFF = utils.parseTable(maskFile, '\t')
        else:
            print("MASK MUST BE A .gff or .bed FILE")
            sys.exit()
        maskCollection = utils.gffToLocusCollection(maskGFF)

        # now mask the reference loci
        referenceLoci = referenceCollection.getLoci()
        filteredLoci = [locus for locus in referenceLoci if len(maskCollection.getOverlap(locus, 'both')) == 0]
        print("FILTERED OUT %s LOCI THAT WERE MASKED IN %s" % (len(referenceLoci) - len(filteredLoci), maskFile))
        referenceCollection = utils.LocusCollection(filteredLoci, 50)

    # NOW STITCH REGIONS
    print('STITCHING REGIONS TOGETHER')
    stitchedCollection, debugOutput, stitchWindow = regionStitching(referenceCollection, inputName, outFolder, stitchWindow, tssWindow, annotFile, removeTSS)

    # NOW MAKE A STITCHED COLLECTION GFF
    print('MAKING GFF FROM STITCHED COLLECTION')
    stitchedGFF = utils.locusCollectionToGFF(stitchedCollection)
    # making sure start/stop ordering are correct
    for i in range(len(stitchedGFF)):

        line = stitchedGFF[i]
        start = int(line[3])
        stop = int(line[4])
        if start > stop:
            line[3] = stop
            line[4] = start

    print(stitchWindow)
    print(type(stitchWindow))
    if not removeTSS:
        stitchedGFFFile = '%s%s_%sKB_STITCHED.gff' % (gffFolder, inputName, str(stitchWindow / 1000))
        stitchedGFFName = '%s_%sKB_STITCHED' % (inputName, str(stitchWindow / 1000))
        debugOutFile = '%s%s_%sKB_STITCHED.debug' % (gffFolder, inputName, str(stitchWindow / 1000))
    else:
        stitchedGFFFile = '%s%s_%sKB_STITCHED_TSS_DISTAL.gff' % (gffFolder, inputName, str(stitchWindow / 1000))
        stitchedGFFName = '%s_%sKB_STITCHED_TSS_DISTAL' % (inputName, str(stitchWindow / 1000))
        debugOutFile = '%s%s_%sKB_STITCHED_TSS_DISTAL.debug' % (gffFolder, inputName, str(stitchWindow / 1000))

    # WRITING DEBUG OUTPUT TO DISK

    if debug:
        print('WRITING DEBUG OUTPUT TO DISK AS %s' % (debugOutFile))
        utils.unParseTable(debugOutput, debugOutFile, '\t')

    # WRITE THE GFF TO DISK
    print('WRITING STITCHED GFF TO DISK AS %s' % (stitchedGFFFile))
    utils.unParseTable(stitchedGFF, stitchedGFFFile, '\t')

    # SETTING UP THE OVERALL OUTPUT FILE
    outputFile1 = outFolder + stitchedGFFName + '_ENHANCER_REGION_MAP.txt'
    print('OUTPUT WILL BE WRITTEN TO  %s' % (outputFile1))

    # MAPPING TO THE NON STITCHED (ORIGINAL GFF)
    # MAPPING TO THE STITCHED GFF

    # Try to use the bamliquidatior_path.py script on cluster, otherwise, failover to local (in path), otherwise fail.
    bamliquidator_path = 'bamliquidator_batch'


    bamFileListUnique = list(bamFileList)
    bamFileListUnique = utils.uniquify(bamFileListUnique)
    #prevent redundant mapping
    print("MAPPING TO THE FOLLOWING BAMS:")
    print(bamFileListUnique)
    for bamFile in bamFileListUnique:

        bamFileName = bamFile.split('/')[-1]

        # MAPPING TO THE STITCHED GFF
        mappedOut1Folder = '%s%s_%s_MAPPED' % (mappedFolder, stitchedGFFName, bamFileName)
        mappedOut1File = '%s%s_%s_MAPPED/matrix.txt' % (mappedFolder, stitchedGFFName, bamFileName)
        if utils.checkOutput(mappedOut1File, 0.2, 0.2):
            print("FOUND %s MAPPING DATA FOR BAM: %s" % (stitchedGFFFile, mappedOut1File))
        else:
            cmd1 = bamliquidator_path + " --sense . -e 200 --match_bamToGFF -r %s -o %s %s" % (stitchedGFFFile, mappedOut1Folder, bamFile)
            print(cmd1)

            os.system(cmd1)
            if utils.checkOutput(mappedOut1File,0.2,5):
                print("SUCCESSFULLY MAPPED TO %s FROM BAM: %s" % (stitchedGFFFile, bamFileName))
            else:
                print("ERROR: FAILED TO MAP %s FROM BAM: %s" % (stitchedGFFFile, bamFileName))
                sys.exit()

    print('BAM MAPPING COMPLETED NOW MAPPING DATA TO REGIONS')
    # CALCULATE DENSITY BY REGION
    # NEED TO FIX THIS FUNCTION TO ACCOUNT FOR DIFFERENT OUTPUTS OF LIQUIDATOR
    mapCollection(stitchedCollection, referenceCollection, bamFileList, mappedFolder, outputFile1, refName=stitchedGFFName)

    print('CALLING AND PLOTTING SUPER-ENHANCERS')

    if options.control:

        rankbyName = options.rankby.split('/')[-1]
        controlName = options.control.split('/')[-1]
        cmd = 'R --no-save %s %s %s %s < ROSE2_callSuper.R' % (outFolder, outputFile1, inputName, controlName)

    else:
        rankbyName = options.rankby.split('/')[-1]
        controlName = 'NONE'
        cmd = 'R --no-save %s %s %s %s < ROSE2_callSuper.R' % (outFolder, outputFile1, inputName, controlName)
    print(cmd)

    os.system(cmd)

    # calling the gene mapper
    time.sleep(20)
    superTableFile = "%s_SuperEnhancers.table.txt" % (inputName)
    if options.control:
        cmd = "python ROSE2_geneMapper.py -g %s -r %s -c %s -i %s%s &" % (genome, options.rankby, options.control, outFolder, superTableFile)
    else:
        cmd = "python ROSE2_geneMapper.py -g %s -r %s -i %s%s &" % (genome, options.rankby, outFolder, superTableFile)
    os.system(cmd)


    stretchTableFile = "%s_StretchEnhancers.table.txt" % (inputName)
    if options.control:
        cmd = "python ROSE2_geneMapper.py -g %s -r %s -c %s -i %s%s &" % (genome, options.rankby, options.control, outFolder, stretchTableFile)
    else:
        cmd = "python ROSE2_geneMapper.py -g %s -r %s -i %s%s &" % (genome, options.rankby, outFolder, stretchTableFile)
    os.system(cmd)


    superStretchTableFile = "%s_SuperStretchEnhancers.table.txt" % (inputName)
    if options.control:
        cmd = "python ROSE2_geneMapper.py -g %s -r %s -c %s -i %s%s &" % (genome, options.rankby, options.control, outFolder, superStretchTableFile)
    else:
        cmd = "python ROSE2_geneMapper.py -g %s -r %s -i %s%s &" % (genome, options.rankby, outFolder, superStretchTableFile)
    os.system(cmd)





if __name__ == "__main__":
    main()
