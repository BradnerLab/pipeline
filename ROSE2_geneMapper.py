#!/usr/bin/python

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


# 140306_singleGeneMapper.py
# main method wrapped script to take the enhancer region table output of ROSE_Main and map genes to it
# will create two outputs a gene mapped region table where each row is an enhancer
# and a gene table where each row is a gene
# does this by default for super-enhancers only
# update to the gene mapper that finds nearest gene w/ highest signal
# also switching to using the pipeline utils module as opposed to the
# stripped down ROSE_utils module
import sys

import utils
# import pipeline_dfci

import os
import subprocess
from string import join

from collections import defaultdict



#==================================================================
#=========================GLOBAL===================================
#==================================================================

# Get the script's full local path
whereAmI = os.path.dirname(os.path.realpath(__file__))
print(whereAmI)
# Get the script folder
codeFolder = utils.formatFolder(whereAmI,False)

print('RUNNING ROSE2_META.py FROM %s' % (whereAmI))


#==================================================================
#===========MAPPING GENES TO ENHANCERS WITHOUT BAM RANKING=========
#==================================================================

#this is the traditional way of running gene mapper


def mapEnhancerToGene(annotFile,enhancerFile,transcribedFile='',uniqueGenes=True,searchWindow =50000,noFormatTable = False):
    
    '''
    maps genes to enhancers. if uniqueGenes, reduces to gene name only. Otherwise, gives for each refseq
    '''
    startDict = utils.makeStartDict(annotFile)
    enhancerTable = utils.parseTable(enhancerFile,'\t')

    #internal parameter for debugging
    byRefseq = False


    if len(transcribedFile) > 0:
        transcribedTable = utils.parseTable(transcribedFile,'\t')
        transcribedGenes = [line[1] for line in transcribedTable]
    else:
        transcribedGenes = startDict.keys()

    print('MAKING TRANSCRIPT COLLECTION')
    transcribedCollection = utils.makeTranscriptCollection(annotFile,0,0,500,transcribedGenes)


    print('MAKING TSS COLLECTION')
    tssLoci = []
    for geneID in transcribedGenes:
        tssLoci.append(utils.makeTSSLocus(geneID,startDict,0,0))


    #this turns the tssLoci list into a LocusCollection
    #50 is the internal parameter for LocusCollection and doesn't really matter
    tssCollection = utils.LocusCollection(tssLoci,50)

    

    geneDict = {'overlapping':defaultdict(list),'proximal':defaultdict(list)}

    #dictionaries to hold ranks and superstatus of gene nearby enhancers
    rankDict = defaultdict(list)
    superDict= defaultdict(list)

    #list of all genes that appear in this analysis
    overallGeneList = []

    if noFormatTable:
        #set up the output tables
        #first by enhancer
        enhancerToGeneTable = [enhancerTable[0]+['OVERLAP_GENES','PROXIMAL_GENES','CLOSEST_GENE']]

        
    else:
        #set up the output tables
        #first by enhancer
        enhancerToGeneTable = [enhancerTable[0][0:9]+['OVERLAP_GENES','PROXIMAL_GENES','CLOSEST_GENE'] + enhancerTable[5][-2:]]

        #next by gene
        geneToEnhancerTable = [['GENE_NAME','REFSEQ_ID','PROXIMAL_ENHANCERS']]

    #next make the gene to enhancer table
    geneToEnhancerTable = [['GENE_NAME','REFSEQ_ID','PROXIMAL_ENHANCERS','ENHANCER_RANKS','IS_SUPER']]

        


    for line in enhancerTable:
        if line[0][0] =='#' or line[0][0] == 'R':
            continue

        enhancerString = '%s:%s-%s' % (line[1],line[2],line[3])
        
        enhancerLocus = utils.Locus(line[1],line[2],line[3],'.',line[0])

        #overlapping genes are transcribed genes whose transcript is directly in the stitchedLocus         
        overlappingLoci = transcribedCollection.getOverlap(enhancerLocus,'both')           
        overlappingGenes =[]
        for overlapLocus in overlappingLoci:                
            overlappingGenes.append(overlapLocus.ID())

        #proximalGenes are transcribed genes where the tss is within 50kb of the boundary of the stitched loci
        proximalLoci = tssCollection.getOverlap(utils.makeSearchLocus(enhancerLocus,searchWindow,searchWindow),'both')           
        proximalGenes =[]
        for proxLocus in proximalLoci:
            proximalGenes.append(proxLocus.ID())


        distalLoci = tssCollection.getOverlap(utils.makeSearchLocus(enhancerLocus,1000000,1000000),'both')           
        distalGenes =[]
        for proxLocus in distalLoci:
            distalGenes.append(proxLocus.ID())

            
            
        overlappingGenes = utils.uniquify(overlappingGenes)
        proximalGenes = utils.uniquify(proximalGenes)
        distalGenes = utils.uniquify(distalGenes)
        allEnhancerGenes = overlappingGenes + proximalGenes + distalGenes
        #these checks make sure each gene list is unique.
        #technically it is possible for a gene to be overlapping, but not proximal since the
        #gene could be longer than the 50kb window, but we'll let that slide here
        for refID in overlappingGenes:
            if proximalGenes.count(refID) == 1:
                proximalGenes.remove(refID)

        for refID in proximalGenes:
            if distalGenes.count(refID) == 1:
                distalGenes.remove(refID)


        #Now find the closest gene
        if len(allEnhancerGenes) == 0:
            closestGene = ''
        else:
            #get enhancerCenter
            enhancerCenter = (int(line[2]) + int(line[3]))/2

            #get absolute distance to enhancer center
            distList = [abs(enhancerCenter - startDict[geneID]['start'][0]) for geneID in allEnhancerGenes]
            #get the ID and convert to name
            closestGene = startDict[allEnhancerGenes[distList.index(min(distList))]]['name']

        #NOW WRITE THE ROW FOR THE ENHANCER TABLE
        if noFormatTable:

            newEnhancerLine = list(line)
            newEnhancerLine.append(join(utils.uniquify([startDict[x]['name'] for x in overlappingGenes]),','))
            newEnhancerLine.append(join(utils.uniquify([startDict[x]['name'] for x in proximalGenes]),','))
            newEnhancerLine.append(closestGene)

        else:
            newEnhancerLine = line[0:9]
            newEnhancerLine.append(join(utils.uniquify([startDict[x]['name'] for x in overlappingGenes]),','))
            newEnhancerLine.append(join(utils.uniquify([startDict[x]['name'] for x in proximalGenes]),','))
            newEnhancerLine.append(closestGene)
            newEnhancerLine += line[-2:]

        enhancerToGeneTable.append(newEnhancerLine)
        #Now grab all overlapping and proximal genes for the gene ordered table

        overallGeneList +=overlappingGenes
        for refID in overlappingGenes:
            geneDict['overlapping'][refID].append(enhancerString)
            rankDict[refID].append(int(line[-2]))
            superDict[refID].append(int(line[-1]))
            
        overallGeneList+=proximalGenes
        for refID in proximalGenes:
            geneDict['proximal'][refID].append(enhancerString)
            rankDict[refID].append(int(line[-2]))
            superDict[refID].append(int(line[-1]))



    #End loop through
    
    #Make table by gene
    overallGeneList = utils.uniquify(overallGeneList)  

    #use enhancer rank to order
    rankOrder = utils.order([min(rankDict[x]) for x in overallGeneList])
        
    usedNames = []
    for i in rankOrder:
        refID = overallGeneList[i]
        geneName = startDict[refID]['name']
        if usedNames.count(geneName) > 0 and uniqueGenes == True:

            continue
        else:
            usedNames.append(geneName)
        
        proxEnhancers = geneDict['overlapping'][refID]+geneDict['proximal'][refID]
        
        superStatus = max(superDict[refID])
        enhancerRanks = join([str(x) for x in rankDict[refID]],',')
    
        newLine = [geneName,refID,join(proxEnhancers,','),enhancerRanks,superStatus]
        geneToEnhancerTable.append(newLine)

    #resort enhancerToGeneTable
    if noFormatTable:
        return enhancerToGeneTable,geneToEnhancerTable
    else:
        enhancerOrder = utils.order([int(line[-2]) for line in enhancerToGeneTable[1:]])
        sortedTable = [enhancerToGeneTable[0]]
        for i in enhancerOrder:
            sortedTable.append(enhancerToGeneTable[(i+1)])

        return sortedTable,geneToEnhancerTable





#==================================================================
#===========MAPPING GENES TO ENHANCERS WITH BAM RANKING============
#==================================================================


def makeSignalDict(mappedGFFFile, controlMappedGFFFile=''):
    '''
    makes a signal dict
    '''
    print('\t called makeSignalDict on %s (ctrl: %s)' % (mappedGFFFile, controlMappedGFFFile))
    signalDict = defaultdict(float)

    mappedGFF = utils.parseTable(mappedGFFFile, '\t')
    if len(controlMappedGFFFile) > 0:
        controlGFF = utils.parseTable(controlMappedGFFFile, '\t')

        for i in range(1, len(mappedGFF)):

            signal = float(mappedGFF[i][2]) - float(controlGFF[i][2])
            if signal < 0:
                signal = 0.0
            signalDict[mappedGFF[i][0]] = signal
    else:
        for i in range(1, len(mappedGFF)):
            signal = float(mappedGFF[i][2])
            signalDict[mappedGFF[i][0]] = signal

    return signalDict

#makeSignalDict('../sshfs/x_rose/mm9_TSS_ENHANCER_GENES_-5000_+5000_CONV3_CD4.nomito.rmdup.bam.gff')

def mapEnhancerToGeneTop(rankByBamFile, controlBamFile, genome, annotFile, enhancerFile, transcribedFile='', uniqueGenes=True, searchWindow=50000, noFormatTable=False):
    '''
    maps genes to enhancers. if uniqueGenes, reduces to gene name only. Otherwise, gives for each refseq
    '''
    startDict = utils.makeStartDict(annotFile)
    enhancerName = enhancerFile.split('/')[-1].split('.')[0]
    enhancerTable = utils.parseTable(enhancerFile, '\t')

    # internal parameter for debugging
    byRefseq = False

    if len(transcribedFile) > 0:
        transcribedTable = utils.parseTable(transcribedFile, '\t')
        transcribedGenes = [line[1] for line in transcribedTable]
    else:
        transcribedGenes = startDict.keys()

    print('MAKING TRANSCRIPT COLLECTION')
    transcribedCollection = utils.makeTranscriptCollection(
        annotFile, 0, 0, 500, transcribedGenes)

    print('MAKING TSS COLLECTION')
    tssLoci = []
    for geneID in transcribedGenes:
        tssLoci.append(utils.makeTSSLocus(geneID, startDict, 0, 0))

    # this turns the tssLoci list into a LocusCollection
    # 50 is the internal parameter for LocusCollection and doesn't really
    # matter
    tssCollection = utils.LocusCollection(tssLoci, 50)

    geneDict = {'overlapping': defaultdict(
        list), 'proximal': defaultdict(list)}

    # dictionaries to hold ranks and superstatus of gene nearby enhancers
    rankDict = defaultdict(list)
    superDict = defaultdict(list)

    # list of all genes that appear in this analysis
    overallGeneList = []

    # find the damn header
    for line in enhancerTable:
        if line[0][0] == '#':
            continue
        else:
            header = line
            break

    if noFormatTable:
        # set up the output tables
        # first by enhancer
        enhancerToGeneTable = [
            header + ['OVERLAP_GENES', 'PROXIMAL_GENES', 'CLOSEST_GENE']]

    else:
        # set up the output tables
        # first by enhancer
        enhancerToGeneTable = [
            header[0:9] + ['OVERLAP_GENES', 'PROXIMAL_GENES', 'CLOSEST_GENE'] + header[-2:]]

        # next by gene
        geneToEnhancerTable = [
            ['GENE_NAME', 'REFSEQ_ID', 'PROXIMAL_ENHANCERS']]

    # next make the gene to enhancer table
    geneToEnhancerTable = [
        ['GENE_NAME', 'REFSEQ_ID', 'PROXIMAL_ENHANCERS', 'ENHANCER_RANKS', 'IS_SUPER', 'ENHANCER_SIGNAL']]

    for line in enhancerTable:
        if line[0][0] == '#' or line[0][0] == 'R':
            continue

        enhancerString = '%s:%s-%s' % (line[1], line[2], line[3])

        enhancerLocus = utils.Locus(line[1], line[2], line[3], '.', line[0])

        # overlapping genes are transcribed genes whose transcript is directly
        # in the stitchedLocus
        overlappingLoci = transcribedCollection.getOverlap(
            enhancerLocus, 'both')
        overlappingGenes = []
        for overlapLocus in overlappingLoci:
            overlappingGenes.append(overlapLocus.ID())

        # proximalGenes are transcribed genes where the tss is within 50kb of
        # the boundary of the stitched loci
        proximalLoci = tssCollection.getOverlap(
            utils.makeSearchLocus(enhancerLocus, searchWindow, searchWindow), 'both')
        proximalGenes = []
        for proxLocus in proximalLoci:
            proximalGenes.append(proxLocus.ID())

        distalLoci = tssCollection.getOverlap(
            utils.makeSearchLocus(enhancerLocus, 1000000, 1000000), 'both')
        distalGenes = []
        for proxLocus in distalLoci:
            distalGenes.append(proxLocus.ID())

        overlappingGenes = utils.uniquify(overlappingGenes)
        proximalGenes = utils.uniquify(proximalGenes)
        distalGenes = utils.uniquify(distalGenes)
        allEnhancerGenes = overlappingGenes + proximalGenes + distalGenes
        # these checks make sure each gene list is unique.
        # technically it is possible for a gene to be overlapping, but not proximal since the
        # gene could be longer than the 50kb window, but we'll let that slide
        # here
        for refID in overlappingGenes:
            if proximalGenes.count(refID) == 1:
                proximalGenes.remove(refID)

        for refID in proximalGenes:
            if distalGenes.count(refID) == 1:
                distalGenes.remove(refID)

        # Now find the closest gene
        if len(allEnhancerGenes) == 0:
            closestGene = ''
        else:
            # get enhancerCenter
            enhancerCenter = (int(line[2]) + int(line[3])) / 2

            # get absolute distance to enhancer center
            distList = [abs(enhancerCenter - startDict[geneID]['start'][0])
                        for geneID in allEnhancerGenes]
            # get the ID and convert to name
            closestGene = startDict[
                allEnhancerGenes[distList.index(min(distList))]]['name']

        # NOW WRITE THE ROW FOR THE ENHANCER TABLE
        if noFormatTable:

            newEnhancerLine = list(line)
            newEnhancerLine.append(
                join(utils.uniquify([startDict[x]['name'] for x in overlappingGenes]), ','))
            newEnhancerLine.append(
                join(utils.uniquify([startDict[x]['name'] for x in proximalGenes]), ','))
            newEnhancerLine.append(closestGene)

        else:
            newEnhancerLine = line[0:9]
            newEnhancerLine.append(
                join(utils.uniquify([startDict[x]['name'] for x in overlappingGenes]), ','))
            newEnhancerLine.append(
                join(utils.uniquify([startDict[x]['name'] for x in proximalGenes]), ','))
            newEnhancerLine.append(closestGene)
            newEnhancerLine += line[-2:]

        enhancerToGeneTable.append(newEnhancerLine)
        # Now grab all overlapping and proximal genes for the gene ordered
        # table

        overallGeneList += overlappingGenes
        for refID in overlappingGenes:
            geneDict['overlapping'][refID].append(enhancerString)
            rankDict[refID].append(int(line[-2]))
            superDict[refID].append(int(line[-1]))

        overallGeneList += proximalGenes
        for refID in proximalGenes:
            geneDict['proximal'][refID].append(enhancerString)
            rankDict[refID].append(int(line[-2]))
            superDict[refID].append(int(line[-1]))

    # End loop through
    # Make table by gene
    print('MAKING ENHANCER ASSOCIATED GENE TSS COLLECTION')
    overallGeneList = utils.uniquify(overallGeneList)

    #get the chromLists from the various bams here
    cmd = 'samtools idxstats %s' % (rankByBamFile)
    idxStats = subprocess.Popen(cmd,stdout=subprocess.PIPE,shell=True)
    idxStats= idxStats.communicate()
    bamChromList = [line.split('\t')[0] for line in idxStats[0].split('\n')[0:-2]]
    
    if len(controlBamFile) > 0:
        cmd = 'samtools idxstats %s' % (controlBamFile)
        idxStats = subprocess.Popen(cmd,stdout=subprocess.PIPE,shell=True)
        idxStats= idxStats.communicate()
        bamChromListControl = [line.split('\t')[0] for line in idxStats[0].split('\n')[0:-2]]
        bamChromList = [chrom for chrom in bamChromList if bamChromListControl.count(chrom) != 0]



    #now make sure no genes have a bad chrom 
    overallGeneList = [gene for gene in overallGeneList if bamChromList.count(startDict[gene]['chr']) != 0]

    
    #now make an enhancer collection of all transcripts    
    enhancerGeneCollection = utils.makeTranscriptCollection(
        annotFile, 5000, 5000, 500, overallGeneList)

    enhancerGeneGFF = utils.locusCollectionToGFF(enhancerGeneCollection)

    # dump the gff to file
    enhancerFolder = utils.getParentFolder(enhancerFile)
    gffRootName = "%s_TSS_ENHANCER_GENES_-5000_+5000" % (genome)
    enhancerGeneGFFFile = "%s%s_%s.gff" % (enhancerFolder, enhancerName,gffRootName)
    utils.unParseTable(enhancerGeneGFF, enhancerGeneGFFFile, '\t')

    # now we need to run bamToGFF

    # Try to use the bamliquidatior_path.py script on cluster, otherwise, failover to local (in path), otherwise fail.
    bamliquidator_path = 'bamliquidator_batch'


    print('MAPPING SIGNAL AT ENHANCER ASSOCIATED GENE TSS')
    # map density at genes in the +/- 5kb tss region
    # first on the rankBy bam
    bamName = rankByBamFile.split('/')[-1]
    mappedRankByFolder = "%s%s_%s_%s/" % (enhancerFolder, enhancerName,gffRootName, bamName)
    mappedRankByFile = "%s%s_%s_%s/matrix.txt" % (enhancerFolder,enhancerName, gffRootName, bamName)
    cmd = bamliquidator_path + ' --sense . -e 200 --match_bamToGFF -r %s -o %s %s' % (enhancerGeneGFFFile, mappedRankByFolder,rankByBamFile)
    print("Mapping rankby bam %s" % (rankByBamFile))
    print(cmd)
    os.system(cmd)

    #check for completion
    if utils.checkOutput(mappedRankByFile,0.2,5):
        print("SUCCESSFULLY MAPPED TO %s FROM BAM: %s" % (enhancerGeneGFFFile, rankByBamFile))
    else:
        print("ERROR: FAILED TO MAP %s FROM BAM: %s" % (enhancerGeneGFFFile, rankByBamFile))
        sys.exit()

    # next on the control bam if it exists
    if len(controlBamFile) > 0:
        controlName = controlBamFile.split('/')[-1]
        mappedControlFolder = "%s%s_%s_%s/" % (
            enhancerFolder, enhancerName,gffRootName, controlName)
        mappedControlFile = "%s%s_%s_%s/matrix.txt" % (
            enhancerFolder, enhancerName,gffRootName, controlName)
        cmd = bamliquidator_path + ' --sense . -e 200 --match_bamToGFF -r %s -o %s %s' % (enhancerGeneGFFFile, mappedControlFolder,controlBamFile)
        print("Mapping control bam %s" % (controlBamFile))
        print(cmd)
        os.system(cmd)

        #check for completion
        if utils.checkOutput(mappedControlFile,0.2,5):
            print("SUCCESSFULLY MAPPED TO %s FROM BAM: %s" % (enhancerGeneGFFFile, controlBamFile))
        else:
            print("ERROR: FAILED TO MAP %s FROM BAM: %s" % (enhancerGeneGFFFile, controlBamFile))
            sys.exit()

    # now get the appropriate output files
    if len(controlBamFile) > 0:
        print("CHECKING FOR MAPPED OUTPUT AT %s AND %s" %
              (mappedRankByFile, mappedControlFile))
        if utils.checkOutput(mappedRankByFile, 1, 1) and utils.checkOutput(mappedControlFile, 1, 1):
            print('MAKING ENHANCER ASSOCIATED GENE TSS SIGNAL DICTIONARIES')
            signalDict = makeSignalDict(mappedRankByFile, mappedControlFile)
        else:
            print("NO MAPPING OUTPUT DETECTED")
            sys.exit()
    else:
        print("CHECKING FOR MAPPED OUTPUT AT %s" % (mappedRankByFile))
        if utils.checkOutput(mappedRankByFile, 1, 30):
            print('MAKING ENHANCER ASSOCIATED GENE TSS SIGNAL DICTIONARIES')
            signalDict = makeSignalDict(mappedRankByFile)
        else:
            print("NO MAPPING OUTPUT DETECTED")
            sys.exit()

    # use enhancer rank to order

    rankOrder = utils.order([min(rankDict[x]) for x in overallGeneList])

    usedNames = []

    # make a new dict to hold TSS signal by max per geneName
    geneNameSigDict = defaultdict(list)
    print('MAKING GENE TABLE')
    for i in rankOrder:
        refID = overallGeneList[i]
        geneName = startDict[refID]['name']
        if usedNames.count(geneName) > 0 and uniqueGenes == True:
            continue
        else:
            usedNames.append(geneName)

        proxEnhancers = geneDict['overlapping'][
            refID] + geneDict['proximal'][refID]

        superStatus = max(superDict[refID])
        enhancerRanks = join([str(x) for x in rankDict[refID]], ',')

        enhancerSignal = signalDict[refID]
        geneNameSigDict[geneName].append(enhancerSignal)

        newLine = [geneName, refID, join(
            proxEnhancers, ','), enhancerRanks, superStatus, enhancerSignal]
        geneToEnhancerTable.append(newLine)
    #utils.unParseTable(geneToEnhancerTable,'/grail/projects/newRose/geneMapper/foo.txt','\t')
    print('MAKING ENHANCER TO TOP GENE TABLE')

    if noFormatTable:
        enhancerToTopGeneTable = [
            enhancerToGeneTable[0] + ['TOP_GENE', 'TSS_SIGNAL']]
    else:
        enhancerToTopGeneTable = [enhancerToGeneTable[0][0:12] + [
            'TOP_GENE', 'TSS_SIGNAL'] + enhancerToGeneTable[0][-2:]]

    for line in enhancerToGeneTable[1:]:

        geneList = []
        if noFormatTable:
            geneList += line[-3].split(',')
            geneList += line[-2].split(',')

        else:
            geneList += line[10].split(',')
            geneList += line[11].split(',')

        geneList = utils.uniquify([x for x in geneList if len(x) > 0])
        if len(geneList) > 0:
            try:
                sigVector = [max(geneNameSigDict[x]) for x in geneList]
                maxIndex = sigVector.index(max(sigVector))
                maxGene = geneList[maxIndex]
                maxSig = sigVector[maxIndex]
                if maxSig == 0.0:
                    maxGene = 'NONE'
                    maxSig = 'NONE'
            except ValueError:
                if len(geneList) == 1:
                    maxGene = geneList[0]
                    maxSig = 'NONE'    
                else:
                    maxGene = 'NONE'
                    maxSig = 'NONE'    
        else:
            maxGene = 'NONE'
            maxSig = 'NONE'
        if noFormatTable:
            newLine = line + [maxGene, maxSig]
        else:
            newLine = line[0:12] + [maxGene, maxSig] + line[-2:]
        enhancerToTopGeneTable.append(newLine)

    # resort enhancerToGeneTable
    if noFormatTable:
        return enhancerToGeneTable, enhancerToTopGeneTable, geneToEnhancerTable
    else:
        enhancerOrder = utils.order([int(line[-2])
                                    for line in enhancerToGeneTable[1:]])
        sortedTable = [enhancerToGeneTable[0]]
        sortedTopGeneTable = [enhancerToTopGeneTable[0]]
        for i in enhancerOrder:
            sortedTable.append(enhancerToGeneTable[(i + 1)])
            sortedTopGeneTable.append(enhancerToTopGeneTable[(i + 1)])

        return sortedTable, sortedTopGeneTable, geneToEnhancerTable


#==================================================================
#=========================MAIN METHOD==============================
#==================================================================
def main():
    '''
    main run call
    '''

    from optparse import OptionParser
    usage = "usage: %prog [options] -g [GENOME] -i [INPUT_ENHANCER_FILE]"
    parser = OptionParser(usage=usage)
    # required flags
    parser.add_option("-i", "--i", dest="input", nargs=1, default=None,
                      help="Enter a ROSE ranked enhancer or super-enhancer file")
    parser.add_option("-g", "--genome", dest="genome", nargs=1, default=None,
                      help="Enter the genome build (MM9,MM8,HG18,HG19)")

    # optional flags
    parser.add_option("-r", "--rankby", dest="rankby", nargs=1, default=None,
                      help="Enter the bam used to rank enhancers")
    parser.add_option("-c", "--control", dest="control", nargs=1, default='',
                      help="Enter a background bam for background correction")

    parser.add_option("-l", "--list", dest="geneList", nargs=1, default=None,
                      help="Enter a gene list to filter through")
    parser.add_option("-o", "--out", dest="out", nargs=1, default=None,
                      help="Enter an output folder. Default will be same folder as input file")
    parser.add_option(
        "-w", "--window", dest="window", nargs=1, default=50000,
        help="Enter a search distance for genes. Default is 50,000bp")
    parser.add_option(
        "-f", "--format", dest="formatTable", action="store_true", default=False,
        help="If flagged, maintains original formatting of input table")

    # RETRIEVING FLAGS
    (options, args) = parser.parse_args()

    if not options.input or not options.genome:

        parser.print_help()
        exit()

    print(options)

    # GETTING THE GENOME
    genome = options.genome
    print('USING %s AS THE GENOME' % genome)

    # GETTING THE CORRECT ANNOT FILE

    genomeDict = {
        'HG18': '%s/annotation/hg18_refseq.ucsc' % (codeFolder),
        'MM9': '%s/annotation/mm9_refseq.ucsc' % (codeFolder),
        'HG19': '%s/annotation/hg19_refseq.ucsc' % (codeFolder),
        'MM8': '%s/annotation/mm8_refseq.ucsc' % (codeFolder),
        'MM10': '%s/annotation/mm10_refseq.ucsc' % (codeFolder),
        'RN4': '%s/annotation/rn4_refseq.ucsc' % (codeFolder),
    }

    annotFile = genomeDict[genome.upper()]

    # GETTING THE INPUT
    enhancerFile = options.input
    window = int(options.window)

    # making the out folder if it doesn't exist
    if options.out:
        outFolder = utils.formatFolder(options.out, True)
    else:
        outFolder = join(enhancerFile.split('/')[0:-1], '/') + '/'

    # GETTING BAM INFO
    rankByBamFile = options.rankby
    controlBamFile = options.control

    # CHECK FORMATTING FLAG
    if options.formatTable:
        noFormatTable = True
    else:
        noFormatTable = False

    # GETTING THE TRANSCRIBED LIST
    if options.geneList:

        transcribedFile = options.geneList
    else:
        transcribedFile = ''

    if options.rankby:
        print("MAPPING GENES TO ENHANCERS USING CLOSEST ACTIVE GENE")
        enhancerToGeneTable, enhancerToTopGeneTable, geneToEnhancerTable = mapEnhancerToGeneTop(
            rankByBamFile, controlBamFile, genome, annotFile, enhancerFile, transcribedFile, True, window, noFormatTable)

        # Writing enhancer output
        enhancerFileName = enhancerFile.split('/')[-1].split('.')[0]

        if window != 50000:
            # writing the enhancer table

            out1 = '%s%s_ENHANCER_TO_GENE_%sKB.txt' % (
                outFolder, enhancerFileName, window / 1000)
            print("writing output to %s" % (out1))
            utils.unParseTable(enhancerToGeneTable, out1, '\t')

            # writing enhancer top gene table
            out2 = '%s%s_ENHANCER_TO_TOP_GENE_%sKB.txt' % (
                outFolder, enhancerFileName, window / 1000)
            utils.unParseTable(enhancerToTopGeneTable, out2, '\t')

            # writing the gene table
            out3 = '%s%s_GENE_TO_ENHANCER_%sKB.txt' % (
                outFolder, enhancerFileName, window / 1000)
            utils.unParseTable(geneToEnhancerTable, out3, '\t')
        else:
            # writing the enhancer table
            out1 = '%s%s_ENHANCER_TO_GENE.txt' % (outFolder, enhancerFileName)
            utils.unParseTable(enhancerToGeneTable, out1, '\t')

            # writing the enhancer table
            out2 = '%s%s_ENHANCER_TO_TOP_GENE.txt' % (outFolder, enhancerFileName)
            utils.unParseTable(enhancerToTopGeneTable, out2, '\t')

            # writing the gene table
            out3 = '%s%s_GENE_TO_ENHANCER.txt' % (outFolder, enhancerFileName)
            utils.unParseTable(geneToEnhancerTable, out3, '\t')
    else:
        #do traditional mapping
        print("MAPPING GENES TO ENHANCERS USING PROXIMITY RULE")
        enhancerToGeneTable,geneToEnhancerTable = mapEnhancerToGene(annotFile,enhancerFile,transcribedFile,True,window,noFormatTable)

        #Writing enhancer output
        enhancerFileName = enhancerFile.split('/')[-1].split('.')[0]

        if window != 50000:
            #writing the enhancer table
            out1 = '%s%s_ENHANCER_TO_GENE_%sKB.txt' % (outFolder,enhancerFileName,window/1000)
            utils.unParseTable(enhancerToGeneTable,out1,'\t')

            #writing the gene table
            out2 = '%s%s_GENE_TO_ENHANCER_%sKB.txt' % (outFolder,enhancerFileName,window/1000)
            utils.unParseTable(geneToEnhancerTable,out2,'\t')
        else:
            #writing the enhancer table
            out1 = '%s%s_ENHANCER_TO_GENE.txt' % (outFolder,enhancerFileName)
            utils.unParseTable(enhancerToGeneTable,out1,'\t')

            #writing the gene table
            out2 = '%s%s_GENE_TO_ENHANCER.txt' % (outFolder,enhancerFileName)
            utils.unParseTable(geneToEnhancerTable,out2,'\t')
        
        

if __name__ == "__main__":
    main()
