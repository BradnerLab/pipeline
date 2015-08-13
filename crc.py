#!/user/bin/env python

######################
#
# Transcriptional Regulatory Network Prediction
# James Bradner Lab
# Dana Farber Cancer Institute
# Version 1.0
# 150801
#
######################

'''
The MIT License (MIT)

Copyright (c) 2015 Alexander Federation

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

######################           
# Dependencies       
######################

import sys
import utils
import os
import subprocess

import numpy
import networkx as nx
from networkx.algorithms.clique import find_cliques_recursive

from string import upper
from collections import defaultdict
import pickle

import gzip
import urllib2

######################
# Functions
######################

def downloadDataFiles():
    path = os.path.dirname(__file__)
    annotationLocation = "%s/data/annotation/"%path
    genomeLocation = "%s/data/genome/"%path
    if not os.path.isdir(annotationLocation):
        os.makedirs(annotationLocation)
        _downloadAnnotation(annotationLocation)
    if not os.path.isdir(genomeLocation):
        os.makedirs(genomeLocation)
        _downloadGenomes(genomeLocation)
        

def downloadFile(location,downloadDir):
    print "Downloading data file:%s"%location 
    data_file = urllib2.urlopen(location)
    out_name = "%s%s"%(downloadDir, location.split("/")[-1])
    print "Downloading to:%s"%(out_name)
    with open(out_name,'w') as out_file:
        out_file.write(data_file.read())
        out_file.close()

def _downloadAnnotation(annotationLocation):
    gitPipeline = "https://raw.githubusercontent.com/BradnerLab/pipeline/master/"
    dataLocations = ["%sannotation/hg18_refseq.ucsc"%(gitPipeline),
                     "%sannotation/hg19_refseq.ucsc"%(gitPipeline),
                     "%sannotation/mm9_refseq.ucsc"%(gitPipeline)]
    for location in dataLocations:
        downloadFile(location,annotationLocation)

def _downloadGenomes(genomeLocation):
    genomes = ['hg19','hg18','mm9']
    for genome in genomes: 
        downloadDir = utils.formatFolder("%s/%s"%(genomeLocation,genome),True)
        filePath = 'http://hgdownload.soe.ucsc.edu/goldenPath/%s/chromosomes/'%genome
        numChromosomes = 23 if 'hg' in genome else 20
        chromosomes = range(1,numChromosomes)+["X"]+["Y"]
        for i in chromosomes:
            downloadName = "chr%s.fa.gz"%i
            downloadFile("%s/%s"%(filePath,downloadName), downloadDir)
            zipName = "%s/%s"%(downloadDir,downloadName)
            inZip = gzip.open(zipName, 'rb')
            with open(zipName.replace(".gz",""), 'w') as outFile:
                outFile.write(inZip.read())
                outFile.close()
            inZip.close()


def calculatePromoterActivity(annotationFile, bamFile, projectName, projectFolder, refseqToNameDict):
    '''
    calculates the level of acetylation at each TF promoter
    by calling bamliquidator
    '''

    print 'GENERATING AN ACTIVITY TABLE USING CHIP DATA'

    annotTable = utils.parseTable(annotationFile, '\t')
    output = []
    counter = 0

    bam = utils.Bam(bamFile)

    startDict = utils.makeStartDict(annotationFile)

    tssLoci = []
    for gene in startDict:
        tssLoci.append(utils.makeTSSLocus(gene,startDict,2500,2500))
    tssCollection = utils.LocusCollection(tssLoci,50)

    gff = utils.locusCollectionToGFF(tssCollection)

    outputname = projectFolder + projectName + '_TSS.gff'
    utils.unParseTable(gff, outputname, '\t')


    mappingCmd = 'bamliquidator_batch'
    mappingCmd += ' -r ' + outputname
    mappingCmd += ' -o ' + projectFolder + 'bamliquidator'
    mappingCmd += ' -m -e 200 '
    mappingCmd += bamFile

    subprocess.call(mappingCmd, shell=True)

    print mappingCmd      
    return outputname

def createEnhancerLoci(enhancerTable, Enumber='super'):
    '''
    input a rose SuperEnhancer table 
    output a table of Loci of the super enhancers
    '''

    print 'CREATING SUPER LOCUS COLLECTION'

    output = []

    if Enumber == 'super':
        for line in enhancerTable[6:]:
            if line[-1] == '1':
                locus = utils.Locus(line[1], line[2], line[3], '.', 
                                    line[0], (float(line[6])-float(line[7])))
                # the final argument is the total signal at that enhancer
                output.append(locus)
    else:
        end = 6+int(Enumber) # skip 6 line header
        for line in enhancerTable[6:end]:
            locus = utils.Locus(line[1], line[2], line[3], '.', 
                                line[0], (float(line[6])-float(line[7])))
            output.append(locus)

    return output

def createExpressionDict(annotationFile, projectFolder, projectName, refseqToNameDict, expCutoff):
    '''
    input: an activity table with refseq in first column and expression or promoter
    acetylation in second column
    output: a dictionary keyed by refseq that points to activity
    '''

    print 'CREATING EXPRESSION DICTIONARY'

    expresionFilename = projectFolder + 'bamliquidator/matrix.txt'
    expressionTable = utils.parseTable(expresionFilename, '\t')

    expressionDictNM = {}
    expressionDictGene = {}

    for line in expressionTable[1:]:
        trid = line[0]
        geneName = refseqToNameDict[trid]
        exp = float(line[2])

        # Save the expression value of each NMid in a dict, keep higher value if multiple
        if trid in expressionDictNM and exp > expressionDictNM[trid]:
            expressionDictNM[trid] = exp
        elif trid not in expressionDictNM:
            expressionDictNM[trid] = exp

        # Save the value of the expression if it's the highest for that gene
        if geneName in expressionDictGene and exp > expressionDictGene[geneName]:
            expressionDictGene[geneName] = exp
        elif geneName not in expressionDictGene:
            expressionDictGene[geneName] = exp

    cutoff = numpy.percentile(expressionDictGene.values(), expCutoff)
    print 'Expression cutoff: ' + str(cutoff)

    expressedGenes = []
    expressedNM = []

    for nmid in expressionDictNM:
        if float(expressionDictNM[nmid]) > cutoff:
            expressedGenes.append(refseqToNameDict[nmid])
            expressedNM.append(nmid)

    expressedGenes = utils.uniquify(expressedGenes)
    Genefilename = projectFolder + projectName + '_EXPRESSED_GENES.txt'
    utils.unParseTable(expressedGenes, Genefilename, '')

    expressedNM = utils.uniquify(expressedNM)
    NMfilename = projectFolder + projectName + '_EXPRESSED_NM.txt'
    utils.unParseTable(expressedNM, NMfilename, '')
      
    return expressedNM, expressionDictNM

def findCanidateTFs(annotationFile, enhancerLoci, expressedNM, expressionDictNM,
                    bamFile, TFlist, refseqToNameDict, projectFolder, projectName, promoter):
    '''                                                           
    Assign each Super-Enhancer to the closest active TSS to its center
    Return a dictionary keyed by TF that points to a list of loci
    '''

    print 'FINDING CANIDATE TFs'

    enhancerAssignment = []
    TFtoEnhancerDict = defaultdict(list)

    # generate the locations of all the TSS
    startDict = utils.makeStartDict(annotationFile)    
    tssLoci = []
    for gene in expressedNM:
        tssLoci.append(utils.makeTSSLocus(gene,startDict,1000,1000))
    tssCollection = utils.LocusCollection(tssLoci,50)    


    # Loop through enhancers
    for enhancer in enhancerLoci:
        
        # If the enhancer overlaps a TSS, save it
        overlappingLoci = tssCollection.getOverlap(enhancer, 'both')
        overlappingGenes =[]
        for overlapLocus in overlappingLoci:
            overlappingGenes.append(overlapLocus.ID())

        # Find all gene TSS within 100 kb
        proximalLoci = tssCollection.getOverlap(utils.makeSearchLocus(enhancer,100000,100000),'both')
        proximalGenes =[]
        for proxLocus in proximalLoci:
            proximalGenes.append(proxLocus.ID())
        
        # If no genes are within 100 kb, find the closest active gene
        closestGene = ''
        if len(overlappingGenes) == 0 and len(proximalGenes) == 0:
        
            distalLoci = tssCollection.getOverlap(utils.makeSearchLocus(enhancer,1000000,1000000),'both')
            distalGenes =[]
            for distalLocus in distalLoci:
                distalGenes.append(distalLocus.ID())

            enhancerCenter = (int(enhancer.start()) + int(enhancer.end())) / 2
            distList = [abs(enhancerCenter - startDict[geneID]['start'][0])
                        for geneID in distalGenes]
            if distList:
                closestGene = distalGenes[distList.index(min(distList))]


        overlappingGenes = utils.uniquify(overlappingGenes)
        proximalGenes = utils.uniquify(proximalGenes)
        for refID in overlappingGenes:
            if proximalGenes.count(refID) == 1:
                proximalGenes.remove(refID)
 

        # If a TSS overlaps an enhancer, assign them together
        if overlappingGenes:
            for gene in overlappingGenes:
                if gene in TFlist:
                    TFtoEnhancerDict[gene].append(enhancer)
                    enhancerAssignment.append([gene, enhancer.chr(), enhancer.start(), enhancer.end(), enhancer.ID()])
                
        # Otherwise, assign the enhancer to the most active gene in 100 kb
        elif not overlappingGenes and proximalGenes:
            highestGene = ''
            highestActivity = 0
            for gene in proximalGenes:
                if expressionDictNM[gene] > highestActivity:
                    highestActivity = expressionDictNM[gene]
                    highestGene = gene
            if highestGene in TFlist:
                TFtoEnhancerDict[gene].append(enhancer)
                enhancerAssignment.append([gene, enhancer.chr(), enhancer.start(), enhancer.end(), enhancer.ID()])
            
        elif not overlappingGenes and not proximalGenes and closestGene:
            if closestGene in TFlist:
                gene = closestGene
                TFtoEnhancerDict[gene].append(enhancer)
                enhancerAssignment.append([gene, enhancer.chr(), enhancer.start(), enhancer.end(), enhancer.ID()])

    # Add promoter is it's not contained in the super
    if promoter:
        for gene in TFtoEnhancerDict.keys():
            promoter = utils.Locus(startDict[gene]['chr'], int(startDict[gene]['start'][0]) - 2000, 
                                   int(startDict[gene]['start'][0]) + 2000, startDict[gene]['sense'])
            overlapBool = False
            for enhancer in TFtoEnhancerDict[gene]:
                if promoter.overlaps(enhancer):
                    overlapBool = True
            if not overlapBool:
                TFtoEnhancerDict[gene].append(promoter)

    seAssignmentFile = projectFolder + projectName + '_ENHANCER_ASSIGNMENT.txt'
    utils.unParseTable(enhancerAssignment, seAssignmentFile, '\t')

    return TFtoEnhancerDict

def formatEnhancerAssignment(TFtoEnhancerDict, refseqToNameDict, projectName, projectFolder):

    '''                                                                             
    takes in the dict mapping TFs to all proximal supers                                     
    returns a file that lists each canidate TFs                                                     
    and gives the coordinates of the super enhancers around them                                    
    '''

    output = [['TF_refseq', 'TF_name', 'chr', 'start', 'stop', 'SuperID', 'Super_Load' ]]

    used = []
 
    for gene in TFtoEnhancerDict.keys():
        for superEnh in TFtoEnhancerDict[gene]:

            check = (refseqToNameDict[gene], superEnh.chr(), superEnh.start(), superEnh.end())
            
            if check not in used:
                newline = [gene, refseqToNameDict[gene]]
                newline.append(superEnh.chr())
                newline.append(superEnh.start())
                newline.append(superEnh.end())
                newline.append(superEnh.ID())
                newline.append(superEnh.score())
                output.append(newline)

                used.append(check)

    outputname = projectFolder + projectName + '_CANIDATE_TF_AND_SUPER_TABLE.txt'

    utils.unParseTable(output, outputname, '\t')

    return 1

def gaussianSmooth(readList, degree=5):
    '''
    Smoothing function for raw bamliquidator output
    '''

    window=degree*2-1
    weight=numpy.array([1.0]*window)
    weightGauss=[]

    for i in range(window):

        i = i-degree+1
        frac = i/float(window)
        gauss = 1/(numpy.exp((4*(frac))**2))
        weightGauss.append(gauss)

    weight=numpy.array(weightGauss)*weight
    smoothed=[0.0]*(len(readList)-window)

    for i in range(len(smoothed)):
        smoothed[i]=sum(numpy.array(readList[i:i+window])*weight)/sum(weight)

    smoothed = [0,0,0,0,0] + smoothed + [0,0,0,0] # return an array of the same length

    return smoothed

def scoreValley(locus, bamFile, projectName, projectFolder, paramater=0.3):
    '''
    calculate valley scores for a locus
    based on this refernce:
    http://bioinformatics.oxfordjournals.org/content/26/17/2071.full
    '''

    nbins = locus.len()/10

    #call bamliquidator on the region and store in a temp file
    os.system('bamliquidator ' + bamFile + ' ' + locus.chr() + ' ' + str(locus.start()) + ' '
              + str(locus.end()) + ' . ' + str(nbins) + ' 0 > ' + projectFolder + 'tempBamliquidator_'
              + projectName + '.txt')

    x = utils.parseTable(projectFolder + 'tempBamliquidator_' + projectName + '.txt', '\t')
    density = [int(y[0]) for y in x]
    smoothDensity =  gaussianSmooth(density, 5)

    scoreArray = []
    regionMax = max(smoothDensity)

    #Now take the smooth reads and calaculate a valley score

    for i in range(len(smoothDensity)):
        score = 0
        try:
            leftmax = max(smoothDensity[i-25:i-10])
        except:
            leftmax = 'edge'
        try:
            rightmax = max(smoothDensity[i+10:i+25])
        except:
            rightmax = 'edge'

        if rightmax == 'edge' and leftmax == 'edge':
            shoulderHeightMin = 0
            shoulderHeightMax = 0
        elif leftmax == 'edge':
            shoulderHeightMin = rightmax
            shoulderHeightMax = rightmax
        elif rightmax == 'edge':
            shoulderHeightMin = leftmax
            shoulderHeightMax = leftmax
        else:
            shoulderHeightMin = min(leftmax, rightmax)
            shoulderHeightMax = max(leftmax, rightmax)

        ratio = (shoulderHeightMax-float(smoothDensity[i]))/regionMax
        if ratio > paramater:
            score = 1
        else:
            score = 0

        scoreArray.append(score)

    return scoreArray

def stitchValleys(valleyList):
    '''
    takes a list of valley loci
    returns a stitched list of valleys to extract seq from
    '''

    valleyCollection = utils.LocusCollection(valleyList,1)
    stitchedValleyCollection = valleyCollection.stitchCollection()
    loci = []
    regions = []
    for valley in stitchedValleyCollection.getLoci():
        if [valley.chr(), valley.start(), valley.end()] not in regions:
            loci.append(valley)
            regions.append([valley.chr(), valley.start(), valley.end()])
    return loci

def findValleys(TFtoEnhancerDict, bamFile, projectName, projectFolder, valleyParamater=0.3, cutoff = 0.2):
    '''
    takes in the super dict
    returns a dictionary of refseqs with all valley loci that are associated
    '''

    print 'IDENTIFYING VALLEYS IN SUPER ENHANCERS'

    valleyBED = []
    valleyDict = {}

    for gene in TFtoEnhancerDict.keys():
        valleyDict[gene] = []
        print gene
        for region in TFtoEnhancerDict[gene]:
            scoreArray = scoreValley(region, bamFile, projectName, projectFolder, valleyParamater)
            for index,score in enumerate(scoreArray):
                if score > cutoff:
                    valley = utils.Locus(region.chr(), region.start() + index*10,
                                         region.start() + (index+1)*10, '.')
                    valleyDict[gene].append(valley)

        stitchedValleys = stitchValleys(valleyDict[gene])
        for valley in stitchedValleys:
            valleyBED.append([valley.chr(), valley.start(), valley.end()])
            valleyDict[gene] = stitchedValleys

    bedfilename = projectFolder + projectName + '_valleys.bed'
    utils.unParseTable(valleyBED, bedfilename, '\t')
    print bedfilename

    return bedfilename

def generateSubpeakFASTA(TFtoEnhancerDict, subpeaks, genomeDirectory, projectName, projectFolder, constExtension):
    '''
    from a BED file of constituents
    generate a FASTA for the consituients contained within the canidate supers
    '''

    subpeakDict = {}
    subpeakBED = [['track name=' + projectName + ' color=204,0,204']]
    subpeakTable = utils.parseTable(subpeaks, '\t')

    subpeakLoci = [utils.Locus(l[0], int(l[1]), int(l[2]), '.') for l in subpeakTable]
    subpeakCollection = utils.LocusCollection(subpeakLoci, 50)

    for gene in TFtoEnhancerDict.keys():
        subpeakDict[gene] = []
        for region in TFtoEnhancerDict[gene]:
            overlaps = subpeakCollection.getOverlap(region)
            extendedOverlaps = [utils.makeSearchLocus(x, constExtension, constExtension) for x in overlaps]

            overlapCollectionTemp = utils.LocusCollection(extendedOverlaps, 50)
            overlapCollection = overlapCollectionTemp.stitchCollection()
            for overlap in overlapCollection.getLoci():
                subpeakBED.append([overlap.chr(), overlap.start(), overlap.end()])
                subpeakDict[gene].append(overlap)

    bedfilename = projectFolder + projectName + '_subpeaks.bed'
    utils.unParseTable(subpeakBED, bedfilename, '\t')

    fasta = []

    for gene in subpeakDict:
        for subpeak in subpeakDict[gene]:

            fastaTitle = gene + '|'  + subpeak.chr() + '|' + str(subpeak.start()) + '|' + str(subpeak.end())
            fastaLine = utils.fetchSeq(genomeDirectory, subpeak.chr(), int(subpeak.start()+1), 
                                       int(subpeak.end()+1))

            fasta.append('>' + fastaTitle)
            fasta.append(upper(fastaLine))

    outname = projectFolder + projectName + '_SUBPEAKS.fa'

    utils.unParseTable(fasta, outname, '')

def findMotifs(canidateGenes, projectFolder, projectName, motifConvertFile, motifDatabaseFile):
    '''
    takes the refseq to subpeak seq dict
    returns the networkx object with all connections
    '''

    # Create a dictionary to call motif names keyed on gene names
    motifDatabase = utils.parseTable(motifConvertFile, '\t')
    motifDatabaseDict = {}
    motifNames = [line[1] for line in motifDatabase]

    # motifDatabaseDict points from gene name to motif code, 
    # allows multiple motifs per gene
    for line in motifDatabase:
        motifDatabaseDict[line[1]] = []
    for line in motifDatabase:
        motifDatabaseDict[line[1]].append(line[0])

    print 'GENERATING TF NETWORK'

    # select the TF candidates that have motifs
    canidateMotifs = []
    for gene in canidateGenes:
        if gene in motifNames:
            canidateMotifs.append(gene)

    print 'Number of annotated canidate TFs that have motifs: ' + str(len(canidateMotifs))
    canidateMotifs = sorted(canidateMotifs)

    # generate FIMO background file
    bgCmd = 'fasta-get-markov -m 1 < ' + projectFolder + projectName + '_SUBPEAKS.fa > ' + projectFolder + projectName + '_bg.meme'
    subprocess.call(bgCmd, shell=True)

    utils.formatFolder(projectFolder + 'FIMO/', True)

    # begin FIMO command here
    fimoCmd = 'fimo'
    for TF in canidateMotifs:
        print TF
        for x in motifDatabaseDict[TF]:
            fimoCmd += ' --motif ' + "'%s'" % (str(x))

    fimoCmd += ' -verbosity 1' 
    fimoCmd += ' -text'
    fimoCmd += ' -oc ' + projectFolder + 'FIMO'
    fimoCmd += ' --bgfile ' + projectFolder + projectName + '_bg.meme'
    fimoCmd += ' ' + motifDatabaseFile + ' '
    fimoCmd += projectFolder + projectName + '_SUBPEAKS.fa'
    fimoCmd += ' > '+ projectFolder + 'FIMO/fimo.txt' 
    print fimoCmd

    fimoOutput = subprocess.call(fimoCmd, shell=True)  #will wait that fimo is done to go on

    return fimoCmd

def buildGraph(projectFolder, projectName, motifConvertFile, refseqToNameDict, canidateGenes):
    '''
    import the FIMO output once it's finished
    build the networkX directed graph
    '''

    motifDatabase = utils.parseTable(motifConvertFile, '\t')
    motifDatabaseDict = {}
    motifNames = [line[1] for line in motifDatabase]

    # The reverse of the other dict, from motif name to gene name
    for line in motifDatabase:
        motifDatabaseDict[line[0]] = line[1]

    fimoFile =  projectFolder + 'FIMO/fimo.txt'
    fimoTable = utils.parseTable(fimoFile, '\t')

    # Use a directed, unweighted graph
    graph = nx.DiGraph(name=projectName)
    graph.add_nodes_from(canidateGenes)

    # MotifDict points from the TF to the genomic location of the motif
    motifDict = defaultdict(list)
    for line in fimoTable[1:]:

        source = motifDatabaseDict[line[0]]   #motifId
        region = line[1].split('|')
        target = refseqToNameDict[region[0]]   #gene name corresponding to the NMid
        graph.add_edge(source, target)
        
        motifDict[source].append((region[1], int(region[2]) + int(line[2]), 
                                  int(region[2]) + int(line[3])))

    utils.formatFolder(projectFolder + 'motifBED/', True)

    # Make a bed file from the motif dict for each factor
    for gene in motifDict.keys():
        if motifDict[gene]:
            bed = []
            for loc in motifDict[gene]:
                bed.append([loc[0], loc[1], loc[2]])

            filename = projectFolder + 'motifBED/' + gene + '_' + projectName + '_motifs.bed'
            utils.unParseTable(bed, filename, '\t')

    return graph

def formatNetworkOutput(graph, projectFolder, projectName, canidateGenes):
    '''
    takes the networkx graph
    returns all figures, tables, etc
    '''

    # output the network as a .ntx dictionary of lists

    networkFilename = projectFolder + projectName + '.ntx'
    networkFile = open(networkFilename, 'w')
    networkDictOfLists = nx.to_dict_of_lists(graph)
    pickle.dump(networkDictOfLists, networkFile)

    # output the adjacency list and nodelist

    nodeFile = projectFolder + projectName + '_NODELIST.txt'
    nodeList = [ [n] for n in graph.nodes_iter()]
    utils.unParseTable(nodeList, nodeFile, '\t')

    adjFile = projectFolder + projectName + '_ADJ_LIST.txt'
    adjList = graph.adjacency_list()
    utils.unParseTable(adjList, adjFile, '\t')

    edgesTable = [['From', 'To']]
    targetList = []
    for i,gene in enumerate(nodeList):
        for j in adjList[i]:
            newline = [gene[0],j]
            edgesTable.append(newline)
            TFname = gene[0]

    edgeFile = projectFolder + projectName + '_EDGE_LIST.txt'
    utils.unParseTable(edgesTable, edgeFile, '\t')


    # Make the degree table
    
    degTable = [['Tf', 'In_Degree', 'Out_Degree', 'Total_Connections' ]]
    degFile = projectFolder + projectName + '_DEGREE_TABLE.txt'

    for node in graph.nodes(): 
        newline = [node, graph.in_degree()[node], graph.out_degree()[node], graph.degree()[node]]
        degTable.append(newline)

    utils.unParseTable(degTable, degFile, '\t')


    # TFs regulating their own enhancer
    autoreg = graph.selfloop_edges()
    selfLoops = [x for x,y in autoreg]
    selfLoopFile = projectFolder + projectName + '_SELF_LOOPS.txt'
    utils.unParseTable(selfLoops, selfLoopFile, '')

    #recover bidirectional edges
    pairs = []
    for n in selfLoops:
        for m in selfLoops:
            if n != m:
                if graph.has_edge(n,m) and graph.has_edge(m,n):
                    pairs.append([n,m])

    # Find cliques
    unDirGraph = nx.from_edgelist(pairs)
    cliqueGen = find_cliques_recursive(unDirGraph)
    cliqueList = list(cliqueGen)

    utils.unParseTable(cliqueList, projectFolder + projectName + '_CLIQUES_ALL.txt', '\t')

    # Find frequency of TF in a clique
    cliqueRanking = []
    factorEnrichmentDict = {}
    outDegreeDict = graph.out_degree()

    for c in cliqueList:
        score = 0
        for gene in c:
            score += outDegreeDict[gene]
            score = score/len(c)
        if score > 0 and len(c) > 2:
            cliqueRanking.append((c, score))

    sortCliqueRanking = sorted(cliqueRanking, reverse=True, key=lambda x:x[1])
    cliqueFile = projectFolder + projectName + '_CLIQUES_DEGREE.txt'
    utils.unParseTable(sortCliqueRanking, cliqueFile, '\t')

    for factor in selfLoops:
        factorEnrichmentDict[factor] = 0
    for pair in cliqueRanking:
        c = pair[0]
        for factor in c:
            factorEnrichmentDict[factor] += 1

    factorRankingTable = []
    for factor in selfLoops:
        newline = [factor, factorEnrichmentDict[factor]/float(len(cliqueRanking))]
        factorRankingTable.append(newline)

    factorRankingFile = projectFolder + projectName + '_ENRICHED_CLIQUE_FACTORS.txt'
    utils.unParseTable(factorRankingTable, factorRankingFile, '\t')

    # Clique ranking algorithm
    print 'Number of cliques:'
    print len(cliqueList)

    #count the occurences of the TFs accross the loops

    dicoTFinloopsCounts={}
    
    for clique in cliqueList:
        for TF in clique:
            if dicoTFinloopsCounts.has_key(TF):
                dicoTFinloopsCounts[TF]+=1
            else:
                dicoTFinloopsCounts[TF]=1

    #calculate a score by loop

    cliqueRanking = []
    cliqueNub = 0
    for clique in cliqueList:
        cliqueScore=0
        for TF in clique:
            cliqueScore = (float(cliqueScore) + (float(dicoTFinloopsCounts[TF])))
        cliqueRanking.append((clique, cliqueScore/len(clique), len(clique)))

    sortCliqueRanking = sorted(cliqueRanking, reverse=True, key=lambda x:x[1])
    cliqueFile = projectFolder + projectName + '_CLIQUES_RANKED.txt'
    utils.unParseTable(sortCliqueRanking, cliqueFile, '\t')

    # Visualizations

    sizeFile = projectFolder + projectName + '_CANIDATE_TF_AND_SUPER_TABLE.txt'
    os.system('Rscript networkScatter.R ' + degFile + ' ' + sizeFile + ' ' +
              projectFolder + projectName + '_NETWORK_SCATTER.pdf')



###################### 
#
# Main Method
#
######################

def main():

    from optparse import OptionParser

    usage = "usage: %prog [options] -e [ENHANCER_FILE] -b [BAM_FILE] -g [GENOME] -o [OUTPUTFOLDER] -n [NAME]" 
    parser = OptionParser(usage = usage)

    #required flags                                                                                                                                        
    parser.add_option("-e","--enhancer_file", dest="enhancers",nargs = 1, default=None,
                      help = "Provide a ROSE generated enhancer table (_AllEnhancers.table.txt)")
    parser.add_option("-b","--bam",dest="bam",nargs =1, default = None,
                      help = "Provide a bam that corresponds to the super enhancer table")
    parser.add_option("-g","--genome",dest="genome",nargs =1, default = None,
                      help = "Provide the build of the genome to be used for the analysis. Currently supports HG19, HG18 and MM9")
    parser.add_option("-o","--output",dest="output",nargs =1, default = None,
                      help = "Enter an output folder")
    parser.add_option("-n","--name",dest="name",nargs =1, default = None,
                      help = "Provide a name for the job")

    #additional options                                                                                                               

    parser.add_option("-s","--subpeaks", dest="subpeaks",nargs=1,default=None,
                      help = "Enter a BED file of regions to search for motifs")                     
    parser.add_option("-x","--expCutoff", dest="expCutoff",nargs=1,default=50,
                      help = "Enter the expression cutoff to be used to define canidate TFs")
    parser.add_option("-l","--extension-length", dest="extension",nargs = 1, default=0,
                      help = "Enter the length to extend subpeak regions for motif finding")
    parser.add_option("-a","--activity", dest="activity",nargs = 1, default=None,
                      help = "A table with refseq in the first column and activity (expression or promoter acetylation) in second")
    parser.add_option("-E","--enhancer_number", dest="Enumber",nargs = 1, default='super',
                      help = "Enter the number of top ranked enhancers to include in the anlaysis. Default is all super-enhancers")
    parser.add_option("--promoter", dest="promoter",nargs = 1, default=False,
                      help = "Enter True if the promoters should be included in the analysis")
    parser.add_option("--motifs", dest="motifs",nargs = 1, default=False,
                      help = "Enter an alternative PWM file for the analysis")
    parser.add_option("-t","--tfs", dest="tfs",nargs=1,default=None,
                      help = "Enter additional TFs (comma separated) to be used in the bindinf analysis")
    parser.add_option("-v","--valley", dest="valley",nargs=1,default=0.3,
                      help = "Paramater for valley calling threshold, default = 0.3")

    (options,args) = parser.parse_args()

    print(options)

    if options.enhancers and options.genome and options.output and options.name:

        ###
        # Define all global file names
        ###

        dataDir = "%s/data/"%(os.path.dirname(__file__))
        annotationDir = "%s/annotation/"%(dataDir)

        
        if options.motifs:
            motifDatabaseFile = options.motifs
        else:
            motifConvertFile = '%s/MotifDictionary.txt'%(annotationDir)
            motifDatabaseFile = '%s/VertebratePWMs.txt'%(annotationDir)

        # User input files
        enhancerFile = options.enhancers
        enhancerTable = utils.parseTable(enhancerFile, '\t')

        if options.bam:
            bamFile = options.bam
            bam = utils.Bam(bamFile)

        genome = options.genome
        genome = upper(genome)

        if genome == 'HG19':
            genomeDirectory = '/grail/genomes/Homo_sapiens/UCSC/hg19/Sequence/Chromosomes/'
            annotationFile = '%s/hg19_refseq.ucsc'%(annotationDir)
            TFfile = '%s/TFlist_NMid_hg19.txt'%(annotationDir)

        if genome == 'HG18':
            genomeDirectory = '/grail/genomes/Homo_sapiens/human_gp_mar_06_no_random/fasta/'
            annotationFile = '%s/hg18_refseq.ucsc'%(annotationDir)
            TFfile = '%s/TFlist_NMid_hg19.txt'%(annotationDir)

        if genome == 'MM9':
            genomeDirectory = '/grail/genomes/Mus_musculus/UCSC/mm9/Sequence/Chromosomes/'
            annotationFile = '%s/mm9_refseq.ucsc'%(annotationDir)
            TFfile = '%s/TFlist_NMid_mm9.txt'%(annotationDir)


        TFtable = utils.parseTable(TFfile, '\t')
        TFlist = [line[0] for line in TFtable]
        TFlistGene = [line[1] for line in TFtable]

        projectFolder = options.output
        projectName = options.name

        if options.subpeaks:
            subpeakFile = options.subpeaks
        else: subpeakFile = None

        refseqToNameDict = {}
        expressionFile = options.activity
        
        if expressionFile:
            expressionTable = utils.parseTable(expressionFile, '\t')

        else:
            expressionTable = calculatePromoterActivity(annotationFile, bamFile, projectName, projectFolder, refseqToNameDict) 

        # Collect paramaters
        expCutoff = int(options.expCutoff)
        constExtension = int(options.extension)
        valleyParamater = float(options.valley)
        enhancerNumber = options.Enumber
        if options.Enumber != 'super':
            enhancerNumber = options.Enumber
        else:
            enhancerNumber = 'super'

        promoter = options.promoter
        additionalTFs = options.tfs

        annotTable = utils.parseTable(annotationFile, '\t')
        for line in annotTable:
            gid = line[1]
            genename = upper(line[12])
            refseqToNameDict[gid] = genename    

        ###
        # Now run all the functions
        ###

        enhancerLoci = createEnhancerLoci(enhancerTable, enhancerNumber)
        expressedNM, expressionDictNM = createExpressionDict(annotationFile, projectFolder, projectName, refseqToNameDict, expCutoff)
        TFtoEnhancerDict = findCanidateTFs(annotationFile, enhancerLoci, expressedNM, expressionDictNM,  bamFile, TFlist, refseqToNameDict, projectFolder, projectName, promoter)
        formatEnhancerAssignment(TFtoEnhancerDict, refseqToNameDict, projectName, projectFolder)
        canidateGenes = [upper(refseqToNameDict[x]) for x in TFtoEnhancerDict.keys()]
        
        if additionalTFs:
            for tf in additionalTFs.split(','):
                canidateGenes.append(tf)
        canidateGenes = utils.uniquify(canidateGenes)

        print canidateGenes

        if subpeakFile == None:
            subpeakFile = findValleys(TFtoEnhancerDict, bamFile, projectName, projectFolder, 
                                      valleyParamater, cutoff = 0.2)
            
        generateSubpeakFASTA(TFtoEnhancerDict, subpeakFile, genomeDirectory, projectName, projectFolder, constExtension)        
        subpeakFile = projectFolder + projectName + '_SUBPEAKS.fa'

        findMotifs(canidateGenes, projectFolder, projectName, motifConvertFile, motifDatabaseFile)
        graph = buildGraph(projectFolder, projectName, motifConvertFile, refseqToNameDict, canidateGenes)
        formatNetworkOutput(graph, projectFolder, projectName, canidateGenes)
        

    else:
        parser.print_help()
        sys.exit()



if __name__ == '__main__':
    main()

