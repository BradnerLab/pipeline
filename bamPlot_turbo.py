#!/usr/bin/env python

# bamPlot.py

'''
The MIT License (MIT)

Copyright (c) 2013 Charles Lin

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
'''


#==========================================================================
#=======================DEPENDENCIES=======================================
#==========================================================================
import argparse
import cPickle
import sys
import utils
import pipeline_dfci
import subprocess
import os
import string
import tempfile
import zlib
from distutils.spawn import find_executable

# Try to use the bamliquidatior script on cluster, otherwise, failover to local default, otherwise fail.
bamliquidatorString = '/ark/home/cl512/pipeline/bamliquidator'
if not os.path.isfile(bamliquidatorString):
    bamliquidatorString = find_executable('bamliquidator')
    if bamliquidatorString is None:
        raise ValueError('bamliquidator not found in path')

# as of now the number of bins to sample the space is hard wired
nBins = 200

# Get the script's full local path
whereAmI = os.path.dirname(os.path.realpath(__file__))

# script that takes in a list of bams, makes an intermediate table file, and then calls R to make the plot

# updated for batching with bamliquidator magic


#==========================================================================
#======================HELPER FUNCTIONS====================================
#==========================================================================


def loadAnnotFile(genome, skip_cache=False):
    """
    load in the annotation and create a geneDict and transcription collection
    """
    genomeDict = {
        'HG18': 'annotation/hg18_refseq.ucsc',
        'MM9': 'annotation/mm9_refseq.ucsc',
        'MM10': 'annotation/mm10_refseq.ucsc',
        'hg18': 'annotation/hg18_refseq.ucsc',
        'mm9': 'annotation/mm9_refseq.ucsc',
        'HG19': 'annotation/hg19_refseq.ucsc',
        'hg19': 'annotation/hg19_refseq.ucsc',
        'hg19_ribo': 'annotation/hg19_refseq.ucsc',
        'HG19_RIBO': 'annotation/hg19_refseq.ucsc',
        'rn4': 'annotation/rn4_refseq.ucsc',
        'RN4': 'annotation/rn4_refseq.ucsc',
        }

    annotFile = whereAmI + '/' + genomeDict[genome]


    if not skip_cache:
        # Try loading from a cache, if the crc32 matches
        annotPathHash = zlib.crc32(annotFile) & 0xFFFFFFFF  # hash the entire location of this script
        annotFileHash = zlib.crc32(open(annotFile, "rb").read()) & 0xFFFFFFFF

        cache_file_name = "%s.%s.%s.cache" % (genome, annotPathHash, annotFileHash)

        cache_file_path = '%s/%s' % (tempfile.gettempdir(), cache_file_name)

        if os.path.isfile(cache_file_path):
            # Cache exists! Load it!
            try:
                print('\tLoading genome data from cache.')
                with open(cache_file_path, 'rb') as cache_fh:
                    cached_data = cPickle.load(cache_fh)
                    print('\tCache loaded.')
                return cached_data
            except (IOError, cPickle.UnpicklingError):
                # Pickle corrupt? Let's get rid of it.
                print('\tWARNING: Cache corrupt or unreadable. Ignoring.')
        else:
            print('\tNo cache exists: Loading annotation (slow).')


    # We're still here, so either caching was disabled, or the cache doesn't exist

    # geneList =['NM_002460','NM_020185']
    geneList = []
    geneDict = utils.makeGenes(annotFile, geneList, True)
    txCollection = utils.makeTranscriptCollection(annotFile, 0, 0, 500, geneList)

    if not skip_cache:
        print('Writing cache for the first time.')
        with open(cache_file_path, 'wb') as cache_fh:
            cPickle.dump((geneDict, txCollection), cache_fh, cPickle.HIGHEST_PROTOCOL)

    return geneDict, txCollection


def tasteTheRainbow(n):
    '''
    samples rainbow color space
    '''
    from colorsys import hsv_to_rgb

    colorList = []
    nRange = [x / float(n) for x in range(0, n)]
    for i in nRange:
        color = [int(255 * x) for x in list(hsv_to_rgb(i, .9, .9))]

        colorList.append(color)
    return colorList


def mapGFFLineToAnnot(gffLine, outFolder, nBins, geneDict, txCollection, sense='both', header=''):
    '''
    for every line produces a file with all of the rectangles to draw
    '''

    if len(header) == 0:
        gffString = '%s_%s_%s_%s' % (gffLine[0], gffLine[6], gffLine[3], gffLine[4])
    else:
        gffString = header
    diagramTable = [[0, 0, 0, 0]]
    nameTable = [['', 0, 0]]
    gffLocus = utils.Locus(gffLine[0], int(gffLine[3]), int(gffLine[4]), gffLine[6], gffLine[1])

    scaleFactor = float(nBins) / gffLocus.len()
    # plotting buffer for diagrams
    plotBuffer = int(gffLocus.len() / float(nBins) * 20)

    overlapLoci = txCollection.getOverlap(gffLocus, sense='both')
    geneList = [locus.ID() for locus in overlapLoci]

    if gffLine[6] == '-':
        refPoint = int(gffLine[4])
    else:
        refPoint = int(gffLine[3])
    offsetCollection = utils.LocusCollection([], 500)
    for geneID in geneList:

        gene = geneDict[geneID]

        print(gene.commonName())
        if len(gene.commonName()) > 1:
            name = gene.commonName()
        else:
            name = geneID
        offset = 4 * len(offsetCollection.getOverlap(gene.txLocus()))
        offsetCollection.append(utils.makeSearchLocus(gene.txLocus(), plotBuffer, plotBuffer))
        # write the name of the gene down
        if gene.sense() == '+':
            geneStart = gene.txLocus().start()
        else:
            geneStart = gene.txLocus().end()
        geneStart = abs(geneStart - refPoint) * scaleFactor
        nameTable.append([name, geneStart, -2 - offset])
        # draw a line across the entire txLocus

        [start, stop] = [abs(x - refPoint) * scaleFactor for x in gene.txLocus().coords()]
        diagramTable.append([start, -0.01 - offset, stop, 0.01 - offset])

        # now draw thin boxes for all txExons
        if len(gene.txExons()) > 0:
            for txExon in gene.txExons():

                [start, stop] = [abs(x - refPoint) * scaleFactor for x in txExon.coords()]

                diagramTable.append([start, -0.5 - offset, stop, 0.5 - offset])

        # now draw fatty boxes for the coding exons if any
        if len(gene.cdExons()) > 0:
            for cdExon in gene.cdExons():

                [start, stop] = [abs(x - refPoint) * scaleFactor for x in cdExon.coords()]

                diagramTable.append([start, -1 - offset, stop, 1 - offset])

    utils.unParseTable(diagramTable, outFolder + gffString + '_diagramTemp.txt', '\t')
    utils.unParseTable(nameTable, outFolder + gffString + '_nameTemp.txt', '\t')


def makeBedCollection(bedFileList):
    '''
    takes in a list of bedFiles and makes a single huge collection
    each locus has as its ID the name of the bed file
    '''

    bedLoci = []
    print("MAKING BED COLLECTION FOR:")
    for bedFile in bedFileList:

        bedName = bedFile.split('/')[-1].split('.')[0]
        print(bedName)
        bed = utils.parseTable(bedFile, '\t')
        for line in bed:
            if len(line) >= 3:
                #check that line[0]
                if line[0][0:3] == 'chr':
                    try:
                        coords = [int(line[1]),int(line[2])]
                        bedLocus = utils.Locus(line[0], min(coords), max(coords), '.', bedName)
                        bedLoci.append(bedLocus)

                    except ValueError:
                        pass

        print("IDENTIFIED %s BED REGIONS" % (len(bedLoci)))

    return utils.LocusCollection(bedLoci, 50)


def mapGFFLineToBed(gffLine, outFolder, nBins, bedCollection, header=''):
    '''
    for every line produces a file with all of the rectangles to draw
    '''

    if len(header) == 0:
        gffString = '%s_%s_%s_%s' % (gffLine[0], gffLine[6], gffLine[3], gffLine[4])
    else:
        gffString = header
    diagramTable = [[0, 0, 0, 0]]
    nameTable = [['', 0, 0]]
    gffLocus = utils.Locus(gffLine[0], int(gffLine[3]), int(gffLine[4]), gffLine[6], gffLine[1])

    scaleFactor = float(nBins) / gffLocus.len()
    # plotting buffer for diagrams
    # plotBuffer = int(gffLocus.len() / float(nBins) * 20) # UNUSED (?)

    overlapLoci = bedCollection.getOverlap(gffLocus, sense='both')
    print("IDENTIFIED %s OVERLAPPING BED LOCI FOR REGION %s" % (len(overlapLoci),gffLine))

    # since beds come from multiple sources, we want to figure out how to offset them
    offsetDict = {}  # this will store each ID name
    bedNamesList = utils.uniquify([locus.ID() for locus in overlapLoci])
    bedNamesList.sort()
    for i in range(len(bedNamesList)):
        offsetDict[bedNamesList[i]] = 2 * i  # offsets different categories of bed regions

    if gffLine[6] == '-':
        refPoint = int(gffLine[4])
    else:
        refPoint = int(gffLine[3])

    # fill out the name table
    for name in bedNamesList:
        offset = offsetDict[name]
        nameTable.append([name, 0, 0.0 - offset])

    for bedLocus in overlapLoci:

        offset = offsetDict[bedLocus.ID()]

        [start, stop] = [abs(x - refPoint) * scaleFactor for x in bedLocus.coords()]

        diagramTable.append([start, -0.5 - offset, stop, 0.5 - offset])

    utils.unParseTable(diagramTable, outFolder + gffString + '_bedDiagramTemp.txt', '\t')
    utils.unParseTable(nameTable, outFolder + gffString + '_bedNameTemp.txt', '\t')


def mapBamToGFFLine(bamFile, MMR, name, gffLine, color, nBins, sense='both', extension=200):
    '''maps reads from a bam to a gff'''

    print('using a MMR/scaling denominator value of %s' % (MMR))

    line = gffLine[0:9]
    gffLocus = utils.Locus(line[0], int(line[3]), int(line[4]), line[6], line[1])

    # setting up the output clusterline
    colorLine = color
    bamName = bamFile.split('/')[-1]
    clusterLine = [bamName, gffLocus.ID(), name, gffLocus.__str__()] + colorLine

    binSize = gffLocus.len() / nBins
    # some regions will be too short to get info on
    # we just kick these back and abandon them
    if binSize == 0:
        clusterLine += ['NA'] * int(nBins)
        return clusterLine

    # flippy flip if sense is negative
    senseTrans = string.maketrans('-+.', '+-+')
    if sense == '-':
        bamSense = string.translate(gffLocus.sense(), senseTrans)
    elif sense == '+':
        bamSense = gffLocus.sense()
    else:
        bamSense = '.'
    # using the bamLiquidator to get the readstring
    # print('using nBin of %s' % nBin)

    bamCommand = "%s %s %s %s %s %s %s %s" % (bamliquidatorString, bamFile, gffLocus.chr(), gffLocus.start(), gffLocus.end(), bamSense, nBins, extension)
    # print(bamCommand)
    getReads = subprocess.Popen(bamCommand, stdin=subprocess.PIPE, stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
    readString = getReads.communicate()
    denList = readString[0].split('\n')[:-1]

    # flip the denList if the actual gff region is -
    if gffLocus.sense() == '-':
        denList = denList[::-1]

    # converting from units of total bp of read sequence per bin to rpm/bp
    denList = [round(float(x) / binSize / MMR, 4) for x in denList]

    clusterLine += denList

    return clusterLine


def callRPlot(summaryFile, outFile, yScale, plotStyle,multi):
    '''
    calls the R plotting thingy
    '''
    if multi == True:
        pageFlag = 'MULTIPLE_PAGE'
    else:
        pageFlag = 'SINGLE_PAGE'

    cmd = 'R --no-save %s %s %s %s %s < %s/bamPlot_turbo.R' % (summaryFile, outFile, yScale, plotStyle, pageFlag,whereAmI)
    print('calling command %s' % (cmd))
    return cmd


def makeBamPlotTables(gff, genome, bamFileList, colorList, nBins, sense, extension, rpm, outFolder, names, title, bedCollection,scale=''):
    '''
    makes a plot table for each line of the gff mapped against all the bams in the bamList
    '''

    # load in the gff
    if type(gff) == str:
        gff = utils.parseTable(gff, '\t')

    # load in the annotation
    print('loading in annotation for %s' % (genome))
    geneDict, txCollection = loadAnnotFile(genome)

    # make an MMR dict so MMRs are only computed once
    print('Getting information about read depth in bams')
    mmrDict = {}

    if len(scale) >0:
        print("Applying scaling factors")
        scaleList = [float(x) for x in scale.split(',')]
    else:
        scaleList = [1]*len(bamFileList)

    #now iterate through the bam files
    for i,bamFile in enumerate(bamFileList):
        # millionMappedReads
        idxCmd = 'samtools idxstats %s' % (bamFile)
        idxPipe = subprocess.Popen(idxCmd, stdin=subprocess.PIPE, stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
        idxStats = idxPipe.communicate()
        idxStats = idxStats[0].split('\n')
        idxStats = [line.split('\t') for line in idxStats]

        rawCount = sum([int(line[2]) for line in idxStats[:-1]])

        #implement scaling
        readScaleFactor = scaleList[i]



        if rpm:
            MMR = round(float(rawCount) / 1000000 / readScaleFactor, 4)
        else:
            MMR = round(1/float(readScaleFactor),4)
        mmrDict[bamFile] = MMR

    ticker = 1
    # go line by line in the gff
    summaryTable = [['DIAGRAM_TABLE', 'NAME_TABLE', 'BED_DIAGRAM_TABLE', 'BED_NAME_TABLE', 'PLOT_TABLE', 'CHROM', 'ID', 'SENSE', 'START', 'END']]
    for gffLine in gff:
        gffString = 'line_%s_%s_%s_%s_%s_%s' % (ticker, gffLine[0], gffLine[1], gffLine[6], gffLine[3], gffLine[4])
        ticker += 1
        print('writing the gene diagram table for region %s' % (gffLine[1]))
        mapGFFLineToAnnot(gffLine, outFolder, nBins, geneDict, txCollection, sense='both', header=gffString)
        mapGFFLineToBed(gffLine, outFolder, nBins, bedCollection, header=gffString)
        outTable = []

        outTable.append(['BAM', 'GENE_ID', 'NAME', 'LOCUSLINE', 'COLOR1', 'COLOR2', 'COLOR3'] + ['bin_' + str(n) for n in range(1, int(nBins) + 1, 1)])

        for i in range(0, len(bamFileList), 1):
            bamFile = bamFileList[i]
            name = names[i]

            color = colorList[i]
            print('getting data for location %s in dataset %s' % (gffLine[1], bamFile))
            mmr = mmrDict[bamFile]
            newLine = mapBamToGFFLine(bamFile, mmr, name, gffLine, color, nBins, sense, extension)

            outTable.append(newLine)

        # get the gene name
        if geneDict.has_key(gffLine[1]):
            geneName = geneDict[gffLine[1]].commonName()
        else:
            geneName = gffLine[1]
        utils.unParseTable(outTable, outFolder + gffString + '_plotTemp.txt', '\t')
        diagramTable = outFolder + gffString + '_diagramTemp.txt'
        plotTable = outFolder + gffString + '_plotTemp.txt'
        nameTable = outFolder + gffString + '_nameTemp.txt'
        bedNameTable = outFolder + gffString + '_bedNameTemp.txt'
        bedDiagramTable = outFolder + gffString + '_bedDiagramTemp.txt'
        summaryTable.append([diagramTable, nameTable, bedDiagramTable, bedNameTable, plotTable, gffLine[0], geneName, gffLine[6], gffLine[3], gffLine[4]])
    summaryTableFileName = "%s%s_summary.txt" % (outFolder, title)
    utils.unParseTable(summaryTable, summaryTableFileName, '\t')
    return summaryTableFileName


def main():
    """
    main run function
    """

    #usage = "usage: %prog [options] -g [GENOME] -b [SORTED BAMFILE(S)] -i [INPUTFILE] -o [OUTPUTFOLDER]"
    parser = argparse.ArgumentParser(usage='%(prog)s [options]')

    # required flags
    parser.add_argument("-b", "--bam", dest="bam", nargs='*',
                        help="Enter a comma separated list of .bam files to be processed.", required=True)
    parser.add_argument("-i", "--input", dest="input", type=str,
                        help="Enter .gff or genomic region e.g. chr1:+:1-1000.", required=True)
    parser.add_argument("-g", "--genome", dest="genome", type=str,
                        help="specify a genome, HG18,HG19,MM8,MM9,MM10 are currently supported", required=True)

    # output flag
    parser.add_argument("-o", "--output", dest="output", type=str,
                        help="Enter the output folder.", required=True)
    # additional options
    parser.add_argument("--stretch-input", dest="stretch_input", default=None, type=int,
                        help="Stretch the input regions to a minimum length in bp, e.g. 10000 (for 10kb)")
    parser.add_argument("-c", "--color", dest="color", default=None,
                        help="Enter a colon separated list of colors e.g. 255,0,0:255,125,0, default samples the rainbow")
    parser.add_argument("-s", "--sense", dest="sense", default='both',
                        help="Map to '+','-' or 'both' strands. Default maps to both.")
    parser.add_argument("-e", "--extension", dest="extension", default=200,
                        help="Extends reads by n bp. Default value is 200bp")
    parser.add_argument("-r", "--rpm", dest="rpm", action='store_true', default=False,
                        help="Normalizes density to reads per million (rpm) Default is False")
    parser.add_argument("-y", "--yScale", dest="yScale", default="relative",
                        help="Choose either relative or uniform y axis scaling. options = 'relative,uniform' Default is relative scaling")
    parser.add_argument("-n", "--names", dest="names", default=None,
                        help="Enter a comma separated list of names for your bams")
    parser.add_argument("-p", "--plot", dest="plot", default="MULTIPLE",
                        help="Choose either all lines on a single plot or multiple plots. options = 'SINGLE,MULTIPLE,MERGE'")
    parser.add_argument("-t", "--title", dest="title", default='',
                        help="Specify a title for the output plot(s), default will be the coordinate region")

    # DEBUG OPTION TO SAVE TEMP FILES
    parser.add_argument("--scale", dest="scale", default='',
                        help="Enter a comma separated list of scaling factors for your bams. Default is none")
    parser.add_argument("--save-temp", dest="save", action='store_true', default=False,
                        help="If flagged will save temporary files made by bamPlot")
    parser.add_argument("--bed", dest="bed",
                        help="Add a space-delimited list of bed files to plot")
    parser.add_argument("--multi-page", dest="multi", action='store_true', default=False,
                        help="If flagged will create a new pdf for each region")

    args = parser.parse_args()

    print(args)

    if args.bam and args.input and args.genome and args.output:

        # Support a legacy mode where a ',' delimited multiple files
        bamFileList = args.bam
        if len(args.bam) == 1:
            bamFileList = args.bam[0].split(',')

        # Make sure these are actually files & readable (!)
        for filename in bamFileList:
            assert(os.access(filename, os.R_OK))

        # bringing in any beds
        if args.bed:
            bedFileList = args.bed
            if type(bedFileList) == str:
                bedFileList = args.bed.split(',')
            print(bedFileList)
            bedCollection = makeBedCollection(bedFileList)
        else:
            bedCollection = utils.LocusCollection([], 50)

        # Load the input for graphing. One of:
        # - A .gff
        # - A .bed
        # - a specific input region (e.g. chr10:.:93150000-93180000)

        valid_sense_options = {'+', '-', '.'}
        if os.access(args.input, os.R_OK):
            if args.input.endswith('.bed'):
                # Uniquely graph every input of this bed
                parsed_input_bed = utils.parseTable(args.input, '\t')
                gffName = os.path.basename(args.input)  # Graph title
                gff = None
                try:
                    if parsed_input_bed[0][5] in valid_sense_options:
                        # This .bed might have a sense parameter
                        gff = [[e[0], '', args.input, e[1], e[2], '', e[5], '', ''] for e in parsed_input_bed]
                except IndexError:
                    pass

                if gff is None:
                    print("Your bed doesn't have a valid senese parameter. Defaulting to both strands, '.'")
                    # We only take chr/start/stop and ignore everything else.
                    gff = [[e[0], '', args.input, e[1], e[2], '', '.', '', ''] for e in parsed_input_bed]
            else:
                # Default to .gff, since that's the original behavior
                gff = utils.parseTable(args.input, '\t')
                gffName = args.input.split('/')[-1].split('.')[0]
        else:
            # means a coordinate line has been given e.g. chr1:+:1-100
            chromLine = args.input.split(':')
            try:
                chrom = chromLine[0]
                sense = chromLine[1]
            except IndexError:
                print('Invalid input line or inaccessible file. Try: chr1:.:1-5000')
                exit()
            assert(sense in valid_sense_options)
            [start, end] = chromLine[2].split('-')
            if chrom[0:3] != 'chr':
                print('ERROR: UNRECOGNIZED GFF OR CHROMOSOME LINE INPUT')
                exit()
            gffLine = [chrom, '', args.input, start, end, '', sense, '', '']
            gffName = "%s_%s_%s_%s" % (chrom, sense, start, end)
            gff = [gffLine]

        # Consider stretching the regions to a fixed minimum size
        if args.stretch_input:
            print('Stretching inputs to a minimum of: %d bp' % (args.stretch_input))
            minLength = args.stretch_input
            stretchGff = []
            for e in gff:
                difference = int(e[4]) - int(e[3])
                if difference < minLength:
                    pad = int((minLength - difference) / 2)
                    stretchGff.append([e[0], e[1], e[2], int(e[3])-pad, int(e[4])+pad, e[5], e[6], e[7], e[8]])
                else:
                    stretchGff.append(e)

            gff = stretchGff

        # Sanity test the gff object
        assert(all([e[6] in valid_sense_options for e in gff]))  # All strands are sane
        assert(all([int(e[3]) < int(e[4]) for e in gff]))  # All start/stops are ordered

        # bring in the genome
        genome = args.genome.upper()
        if ['HG18', 'HG19', 'HG19_RIBO','MM9', 'MM10', 'RN4'].count(genome) == 0:
            print('ERROR: UNSUPPORTED GENOME TYPE %s. USE HG19,HG18, RN4, MM9, or MM10' % (genome))
            parser.print_help()
            exit()

        # bring in the rest of the options

        # output
        rootFolder = args.output
        if rootFolder[-1] != '/':
            rootFolder += '/'
        try:
            os.listdir(rootFolder)
        except OSError:
            print('ERROR: UNABLE TO FIND OUTPUT DIRECTORY %s' % (rootFolder))
            exit()

        # Get analysis title
        if len(args.title) == 0:
            title = gffName
        else:
            title = args.title

        # make a temp folder
        tempFolder = rootFolder + title + '/'
        print("CREATING TEMP FOLDER %s" % (tempFolder))
        pipeline_dfci.formatFolder(tempFolder, create=True)

        # colors
        if args.color:
            colorList = args.color.split(':')
            colorList = [x.split(',') for x in colorList]
            if len(colorList) < len(bamFileList):
                print('WARNING: FEWER COLORS THAN BAMS SPECIFIED. COLORS WILL BE RECYCLED')
                # recycling the color list
                colorList += colorList * (len(bamFileList) / len(colorList))
                colorList = colorList[0:len(bamFileList)]

        else:
            # cycles through the colors of the rainbow
            colorList = tasteTheRainbow(len(bamFileList))

        # sense
        sense = args.sense

        extension = int(args.extension)

        rpm = args.rpm

        scale = args.scale

        yScale = args.yScale.upper()

        # names
        if args.names:
            names = args.names.split(',')

            if len(names) != len(bamFileList):
                print('ERROR: NUMBER OF NAMES AND NUMBER OF BAMS DO NOT CORRESPOND')
                parser.print_help()
                exit()
        else:
            names = [x.split('/')[-1] for x in bamFileList]

        # plot style
        plotStyle = args.plot.upper()
        if ['SINGLE', 'MULTIPLE','MERGE'].count(plotStyle) == 0:
            print('ERROR: PLOT STYLE %s NOT AN OPTION' % (plotStyle))
            parser.print_help()
            exit()

        # now run!
        summaryTableFileName = makeBamPlotTables(gff, genome, bamFileList, colorList, nBins, sense, extension, rpm, tempFolder, names, title, bedCollection,scale)
        print ("%s is the summary table" % (summaryTableFileName))

        #running the R command to plot
        multi = args.multi
        outFile = "%s%s_plots.pdf" % (rootFolder, title)
        rCmd = callRPlot(summaryTableFileName, outFile, yScale, plotStyle,multi)

        # open a bash file
        bashFileName = "%s%s_Rcmd.sh" % (tempFolder, title)
        bashFile = open(bashFileName, 'w')
        bashFile.write('#!/usr/bin/bash\n')
        bashFile.write(rCmd)
        bashFile.close()
        print("Wrote R command to %s" % (bashFileName))
        os.system("bash %s" % (bashFileName))

        # delete temp files
        if not args.save:
            if utils.checkOutput(outFile, 1, 10):
                # This is super dangerous (!). Add some sanity checks.
                assert(" " not in tempFolder)
                assert(tempFolder is not "/")
                removeCommand = "rm -rf %s" % (tempFolder)
                print(removeCommand)
                os.system(removeCommand)
            else:
                print("ERROR: NO OUTPUT FILE %s DETECTED" % (outFile))

    else:
        parser.print_help()
        sys.exit()

if __name__ == "__main__":
    main()
