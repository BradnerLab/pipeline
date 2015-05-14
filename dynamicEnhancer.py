#!/usr/bin/python

#131108_dynamicEnhancer.py
#131108
#Charles Lin


#Description:

'''
pipeline to run dynamic enhancer analysis


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



#================================================================================
#=============================DEPENDENCIES=======================================
#================================================================================

import sys

print "Using python version %s" % sys.version
    

#importing utils package
sys.path.append('/ark/home/cl512/pipeline/')
import utils
import pipeline_dfci
import os
import time
import string
import numpy

#================================================================================
#============================GLOBAL PARAMETERS===================================
#================================================================================

#add locations of files and global parameters in this section

pipelineDir = '/ark/home/cl512/pipeline/'


#dataFile = '/ark/home/cl512/projects/athero/EC_TABLE_FINAL.txt'
#genome = 'hg18'

#dataDict = pipeline_dfci.loadDataTable(dataFile)

#================================================================================
#===================================CLASSES======================================
#================================================================================

#user defined classes here

#================================================================================
#=================================FUNCTIONS======================================
#================================================================================

def getFile(fileString,fileList,parentFolder):
    '''
    returns full path of file from fileList containing the fileString
    returns an error if multiple files match
    '''
    if not utils.formatFolder(parentFolder,False):
        print "ERROR: Folder %s does not exist" % (parentFolder)
        sys.exit()
    parentFolder = utils.formatFolder(parentFolder,False)
    matchFiles = [fileName for fileName in fileList if fileName.count(fileString) == 1]
    if len(matchFiles) == 0:
        print "WARNING: No files found in %s with %s in title" % (parentFolder,fileString)
        return ''
    if len(matchFiles) > 1:
        print "ERROR: Multiple files found in %s with %s in title" % (parentFolder,fileString)
        sys.exit()
    matchFilePath  = "%s%s" % (parentFolder,matchFiles[0])
    return matchFilePath



def makeRoseDict(roseFolder):

    '''
    analyzes a rose folder to try to find all of the various necessary files
    creates a dictionary with their full paths
    '''
    if not utils.formatFolder(roseFolder,False):
        
        print "Folder %s does not exist" % (roseFolder)
        sys.exit()

    roseFolder = utils.formatFolder(roseFolder,False)
    roseFileList = [x for x in os.listdir(roseFolder) if x[0] != '.'] #no hidden files
    if len(roseFileList) == 0:
        print "No files found in %s" % (roseFolder)
        sys.exit()

    #create a dictionary to store stuff
    roseDict = {}
    #there are 5 files that we're interested in
    #REGION_MAP, AllEnhancers.table.txt, SuperEnhancers.table.txt, ENHANCER_TO_GENE, Enhancers_withSuper.bed

    #sequentially find each one and add the full path to the roseDict
    roseDict['AllEnhancer'] = getFile('AllEnhancers.table.txt',roseFileList,roseFolder)
    roseDict['super'] = getFile('SuperEnhancers.table.txt',roseFileList,roseFolder)
    roseDict['stretch'] = getFile('_StretchEnhancers.table.txt',roseFileList,roseFolder)
    roseDict['superstretch'] = getFile('SuperStretchEnhancers.table.txt',roseFileList,roseFolder)

    roseDict['EnhancerToGene'] = getFile('_SuperEnhancers_ENHANCER_TO_GENE',roseFileList,roseFolder)
    roseDict['RegionMap'] = getFile('REGION_MAP',roseFileList,roseFolder)
    roseDict['bed'] = getFile('Enhancers_withSuper.bed',roseFileList,roseFolder)

    return roseDict


def getMedianSignal(enhancerFile,name,dataFile):

    '''
    returns the median enhancer signal of a file
    '''
    dataDict = pipeline_dfci.loadDataTable(dataFile)
    enhancerTable = utils.parseTable(enhancerFile,'\t')

    backgroundName = dataDict[name]['background']
    if dataDict.has_key(backgroundName):
        enhancerVector = [float(line[6]) - float(line[7]) for line in enhancerTable[6:]]
    else:
        enhancerVector = [float(line[6]) for line in enhancerTable[6:]]
            

    median= numpy.median(enhancerVector)

    return median




def makeSECollection(enhancerFile,name,top=0):
    '''
    returns a locus collection from a super table
    top gives the number of rows
    '''
    enhancerTable = utils.parseTable(enhancerFile,'\t')
    superLoci = []

    ticker = 0
    for line in enhancerTable:
        if line[0][0] == '#' or line[0][0] == 'R':
            continue
        else:
            ticker+=1

            superLoci.append(utils.Locus(line[1],line[2],line[3],'.',name+'_'+line[0]))

            if ticker == top:
                break
    return utils.LocusCollection(superLoci,50)

def makeSEDict(enhancerFile,name,superOnly = True):

    '''
    makes an attribute dict for enhancers keyed by uniqueID
    '''

    seDict = {}
    enhancerTable = utils.parseTable(enhancerFile,'\t')

    superLoci = []
    for line in enhancerTable:
        if line[0][0] == '#':
            continue
        if line[0][0] == 'R':
            header = line
            supColumn = header.index('isSuper')
            continue
        if superOnly:
            if int(line[supColumn]) == 1:
                
            
                rank = int(line[-2])
                enhancerID = name+'_'+line[0]
                seDict[enhancerID] = {'rank':rank}

        else:

            signal = float(line[6]) - float(line[7])
            rank = int(line[-2])
            enhancerID = name+'_'+line[0]
            seDict[enhancerID] = {'rank':rank}

    return seDict


def mergeCollections(superFile1,superFile2,name1,name2,output=''):

    '''
    merges them collections
    '''

    conSuperCollection = makeSECollection(superFile1,name1)

    tnfSuperCollection = makeSECollection(superFile2,name2)


    #now merge them
    mergedLoci = conSuperCollection.getLoci() + tnfSuperCollection.getLoci()

    mergedCollection = utils.LocusCollection(mergedLoci,50)

    #stitch the collection together
    stitchedCollection = mergedCollection.stitchCollection()

    stitchedLoci = stitchedCollection.getLoci()
    
    #loci that are in both get renamed with a new unique identifier

    renamedLoci =[]
    ticker = 1
    for locus in stitchedLoci:

        if len(conSuperCollection.getOverlap(locus)) > 0 and len(tnfSuperCollection.getOverlap(locus)):

            newID = 'CONSERVED_%s' % (str(ticker))
            ticker +=1
            locus._ID = newID
        else:
            locus._ID = locus.ID()[2:]
        renamedLoci.append(locus)

    #now we turn this into a gff and write it out
    gff = utils.locusCollectionToGFF(utils.LocusCollection(renamedLoci,50))

    if len(output) == 0:
        return gff
    else:
        print "writing merged gff to %s" % (output)
        utils.unParseTable(gff,output,'\t')
        return output






#call rose on the mergies

def callRoseMerged(dataFile,mergedGFFFile,name1,name2,parentFolder):

    '''
    makes a rose call for the merged supers
    '''

    dataDict = pipeline_dfci.loadDataTable(dataFile)

    
    backgroundName1 = dataDict[name1]['background']
    backgroundName2 = dataDict[name2]['background']
    if dataDict.has_key(backgroundName1) and dataDict.has_key(backgroundName2):
        hasBackground = True
    elif not dataDict.has_key(backgroundName1) and not dataDict.has_key(backgroundName2):
        hasBackground =False
    else:
        print "ERROR: Only 1 dataset has a background file. This is a very very bad idea"
        sys.exit()


    if hasBackground:
        namesList = [name1]    
        extraMap = [name2,dataDict[name2]['background']]
    else:
        namesList = [name1]
        extraMap = [name2]


    return pipeline_dfci.callRose2(dataFile,'',parentFolder,namesList,extraMap,mergedGFFFile,tss=0,stitch=0)


def callMergeSupers(dataFile,superFile1,superFile2,name1,name2,mergeName,genome,parentFolder):

    '''
    this is the main run function for the script
    all of the work should occur here, but no functions should be defined here
    '''
    mergedGFFFile = '%s%s_%s_MERGED_REGIONS_-0_+0.gff' % (parentFolder,string.upper(genome),mergeName)    

    #check to make sure this hasn't been done yet
    roseOutput = "%s%s_ROSE/%s_%s_MERGED_REGIONS_-0_+0_SuperEnhancers_ENHANCER_TO_GENE.txt" % (parentFolder,name1,string.upper(genome),mergeName)

    try:
        foo = utils.parseTable(roseOutput,'\t')
        print "ROSE OUTPUT ALREADY FOUND HERE %s" % (roseOutput)
        return roseOutput
    except IOError:
        
        print "MERGING ENHANCER REGIONS FROM %s and %s" % (superFile1,superFile2)
        mergedGFF = mergeCollections(superFile1,superFile2,name1,name2,mergedGFFFile)

        #call rose on the merged shit    


        roseBashFile = callRoseMerged(dataFile,mergedGFF,name1,name2,parentFolder)
        print('i can has rose bash file %s' % (roseBashFile))

        #run the bash command
        os.system('bash %s' % (roseBashFile))

        #check for and return output
        if utils.checkOutput(roseOutput,1,10):
            return roseOutput
        else:
            #try finding it w/ a different name
            #this will bug out if nothing is there
            roseFolder = "%s%s_ROSE/" % (parentFolder,name1)
            roseFileList = [x for x in os.listdir(roseFolder) if x[0] != '.'] #no hidden files
            if len(roseFileList) == 0:
                print "No files found in %s" % (roseFolder)
                sys.exit()

            enhancerToGeneFile = getFile('_SuperEnhancers_ENHANCER_TO_GENE.txt',roseFileList,roseFolder)
            


def callDeltaRScript(mergedGFFFile,parentFolder,dataFile,name1,name2,allFile1,allFile2,medianScale = False):

    '''
    runs the R script
    '''
    if medianScale:
        median1 = getMedianSignal(allFile1,name1,dataFile)
        median2 = getMedianSignal(allFile2,name2,dataFile)
        print "normalizing signal for %s by median value of %s" % (name1,median1)
        print "normalizing signal for %s by median value of %s" % (name2,median2)

    else:
        median1 =1
        median2 =1

    
    gffName = mergedGFFFile.split('/')[-1].split('.')[0]
    stitchedFile = "%s%s_ROSE/%s_0KB_STITCHED_ENHANCER_REGION_MAP.txt" % (parentFolder,name1,gffName)
    #print(stitchedFile)
    os.chdir(pipelineDir)

    rcmd = "R --no-save %s %s %s %s %s < ./dynamicEnhancer_plot.R" % (stitchedFile,name1,name2,median1,median2)

    return rcmd

def callRankRScript(enhancerRankFile,name1,name2,superFile1,superFile2):

    '''
    runs the R script
    '''

    enhancerCollection1 = makeSECollection(superFile1,name1,False)
    enhancerCollection2 = makeSECollection(superFile2,name2,False)

    nSuper1 = len(enhancerCollection1)
    nSuper2 = len(enhancerCollection2)



    os.chdir(pipelineDir)
    rcmd = "R --no-save %s %s %s %s %s < ./dynamicEnhancer_rank.R" % (enhancerRankFile,name1,name2,nSuper1,nSuper2)

    return rcmd




def callRoseGeneMapper(mergedGFFFile,genome,parentFolder,name1):

    '''
    calls the rose gene mapper w/ 100kb window
    '''
    gffName = mergedGFFFile.split('/')[-1].split('.')[0]
    stitchedFile = "%s%s_ROSE/%s_0KB_STITCHED_ENHANCER_REGION_MAP.txt" % (parentFolder,name1,gffName)
    
    deltaFile = stitchedFile.replace('REGION_MAP','DELTA')
    
    os.chdir(pipelineDir)
    cmd = 'python ROSE2_geneMapper.py -g %s -i %s -w 100000' % (genome,deltaFile)
    os.system(cmd)
    print(cmd)
    


def assignEnhancerRank(enhancerToGeneFile,enhancerFile1,enhancerFile2,name1,name2,rankOutput=''):

    '''
    for all genes in the enhancerToGene Table, assigns the highest overlapping ranked enhancer in the other tables
    '''

    enhancerToGene = utils.parseTable(enhancerToGeneFile,'\t')

    enhancerCollection1 = makeSECollection(enhancerFile1,name1,False)
    enhancerCollection2 = makeSECollection(enhancerFile2,name2,False)

    enhancerDict1 = makeSEDict(enhancerFile1,name1,False)
    enhancerDict2 = makeSEDict(enhancerFile2,name2,False)

    
    #we're going to update the enhancerToGeneTable

    enhancerToGene[0] += ['%s_rank' % name1,'%s_rank' % name2]
    
    for i in range(1,len(enhancerToGene)):

        line = enhancerToGene[i]
        
        locusLine = utils.Locus(line[1],line[2],line[3],'.',line[0])
        
        #if the enhancer doesn't exist, its ranking is dead last on the enhancer list

        enhancer1Overlap = enhancerCollection1.getOverlap(locusLine,'both')
        if len(enhancer1Overlap) == 0:
            enhancer1Rank = len(enhancerCollection1)
        else:
            
            rankList1 = [enhancerDict1[x.ID()]['rank'] for x in enhancer1Overlap]
            enhancer1Rank = min(rankList1)


        enhancer2Overlap = enhancerCollection2.getOverlap(locusLine,'both')
        if len(enhancer2Overlap) == 0:
            enhancer2Rank = len(enhancerCollection2)
        else:
            
            rankList2 = [enhancerDict2[x.ID()]['rank'] for x in enhancer2Overlap]
            enhancer2Rank = min(rankList2)
        enhancerToGene[i]+=[enhancer1Rank,enhancer2Rank]


    if len(rankOutput) == 0:
        return enhancerToGene
    else:
        utils.unParseTable(enhancerToGene,rankOutput,'\t')

#make gain lost gffs

def finishRankOutput(dataFile,rankOutput,genome,mergeFolder,mergeName,name1,name2,cutOff=1.5,window = 100000,superOnly=True,plotBam=True):

    '''
    cleans up the rank output table
    makes a gff of all of the gained/lost supers beyond
    a certain cutoff w/ a window
    makes a list of gained genes and lost genes
    makes a bed of gained loss
    '''
    dataDict = pipeline_dfci.loadDataTable(dataFile)
    #making sure window and cutoff are int/float
    cutOff = float(cutOff)
    window = int(window)
    genome = string.upper(genome)

    #make the output folder
    outputFolder =pipeline_dfci.formatFolder(mergeFolder+'output/',True)
    
    #bring in the old rank table
    rankEnhancerTable = utils.parseTable(rankOutput,'\t')
    
    #make a new formatted table
    header = rankEnhancerTable[0]
    header[-4] = 'DELTA RANK'
    header[-3] = 'IS_SUPER'
    formattedRankTable =[header]

    #the gffs
    gainedGFF = []
    lostGFF = []

    gainedWindowGFF = []
    lostWindowGFF = []

    if superOnly:
        enhancerType = 'SUPERS'
    else:
        enhancerType = 'ENHANCERS'

    #the beds
    if superOnly:
        gainedTrackHeader = 'track name="%s %s only SEs" description="%s super enhancers that are found only in %s vs %s" itemRGB=On color=255,0,0' % (genome,name2,genome,name2,name1)
        gainedBed = [[gainedTrackHeader]]
        conservedTrackHeader = 'track name="%s %s and %s SEs" description="%s super enhancers that are found in both %s vs %s" itemRGB=On color=0,0,0' % (genome,name1,name2,genome,name1,name2)
        conservedBed = [[conservedTrackHeader]]

        lostTrackHeader = 'track name="%s %s only SEs" description="%s super enhancers that are found only in %s vs %s" itemRGB=On color=0,255,0' % (genome,name1,genome,name1,name2)
        lostBed = [[lostTrackHeader]]
    else:
        gainedTrackHeader = 'track name="%s %s only enhancers" description="%s enhancers that are found only in %s vs %s" itemRGB=On color=255,0,0' % (genome,name2,genome,name2,name1)
        gainedBed = [[gainedTrackHeader]]
        conservedTrackHeader = 'track name="%s %s and %s enhancers" description="%s enhancers that are found in both %s vs %s" itemRGB=On color=0,0,0' % (genome,name1,name2,genome,name1,name2)
        conservedBed = [[conservedTrackHeader]]

        lostTrackHeader = 'track name="%s %s only enhancers" description="%s enhancers that are found only in %s vs %s" itemRGB=On color=0,255,0' % (genome,name1,genome,name1,name2)
        lostBed = [[lostTrackHeader]]



    #the genes
    geneTable =[['GENE','ENHANCER_ID','ENHANCER_CHROM','ENHANCER_START','ENHANCER_STOP',header[6],header[7],header[8],'STATUS']]

    for line in rankEnhancerTable[1:]:
        #fixing the enhancer ID
        line[0] = line[0].replace('_lociStitched','')
        formattedRankTable.append(line)

        #getting the genes
        geneList = []
        geneList += line[9].split(',')
        geneList += line[10].split(',')
        geneList += line[11].split(',')
        geneList = [x for x in geneList if len(x) >0]
        geneList = utils.uniquify(geneList)
        geneString = string.join(geneList,',')

        bedLine = [line[1],line[2],line[3],line[0],line[-4]]
        
        #for gained
        if float(line[6]) > cutOff:
            gffLine = [line[1],line[0],'',line[2],line[3],'','.','',geneString]
            gffWindowLine = [line[1],line[0],'',int(line[2])-window,int(line[3])+window,'','.','',geneString]
            gainedGFF.append(gffLine)
            gainedWindowGFF.append(gffWindowLine)
            geneStatus = name2
            gainedBed.append(bedLine)
        #for lost
        elif float(line[6]) < (-1 * cutOff):
            gffLine = [line[1],line[0],'',line[2],line[3],'','.','',geneString]
            gffWindowLine = [line[1],line[0],'',int(line[2])-window,int(line[3])+window,'','.','',geneString]
            lostGFF.append(gffLine)
            lostWindowGFF.append(gffWindowLine)
            geneStatus = name1
            lostBed.append(bedLine)
        #for conserved
        else:
            geneStatus = 'CONSERVED'
            conservedBed.append(bedLine)

        #now fill in the gene Table
        for gene in geneList:
            geneTableLine = [gene,line[0],line[1],line[2],line[3],line[6],line[7],line[8],geneStatus]
            geneTable.append(geneTableLine)

    #concat the bed
    fullBed = gainedBed + conservedBed + lostBed
            
    #start writing the output
    #there's the two gffs, the bed,the formatted table, the gene table
    
    
    #formatted table
    formattedFilename = "%s%s_%s_MERGED_%s_RANK_TABLE.txt" % (outputFolder,genome,mergeName,enhancerType)
    utils.unParseTable(formattedRankTable,formattedFilename,'\t')

    #gffs
    gffFolder = pipeline_dfci.formatFolder(outputFolder+'gff/',True)
    gffFilename_gained = "%s%s_%s_%s_ONLY_%s_-0_+0.gff" % (gffFolder,genome,mergeName,string.upper(name2),enhancerType)
    gffFilenameWindow_gained = "%s%s_%s_%s_ONLY_%s_-%sKB_+%sKB.gff" % (gffFolder,genome,mergeName,string.upper(name2),enhancerType,window/1000,window/1000)

    gffFilename_lost = "%s%s_%s_%s_ONLY_%s_-0_+0.gff" % (gffFolder,genome,mergeName,string.upper(name1),enhancerType)
    gffFilenameWindow_lost = "%s%s_%s_%s_ONLY_%s_-%sKB_+%sKB.gff" % (gffFolder,genome,mergeName,string.upper(name1),enhancerType,window/1000,window/1000)

    utils.unParseTable(gainedGFF,gffFilename_gained,'\t')
    utils.unParseTable(gainedWindowGFF,gffFilenameWindow_gained,'\t')
            
    utils.unParseTable(lostGFF,gffFilename_lost,'\t')
    utils.unParseTable(lostWindowGFF,gffFilenameWindow_lost,'\t')
    
    #bed
    bedFilename = "%s%s_%s_MERGED_%s.bed" % (outputFolder,genome,mergeName,enhancerType)
    utils.unParseTable(fullBed,bedFilename,'\t')

    #geneTable
    geneFilename = "%s%s_%s_MERGED_%s_GENE_TABLE.txt" % (outputFolder,genome,mergeName,enhancerType)
    utils.unParseTable(geneTable,geneFilename,'\t')

    #finally, move all of the plots to the output folder
    cmd = "cp %s%s_ROSE/*.pdf %s%s_%s_MERGED_%s_DELTA.pdf" % (mergeFolder,name1,outputFolder,genome,mergeName,enhancerType)
    os.system(cmd)

    cmd = "cp %s%s_ROSE/*RANK_PLOT.png %s%s_%s_MERGED_%s_RANK_PLOT.png" % (mergeFolder,name1,outputFolder,genome,mergeName,enhancerType)
    os.system(cmd)

    #now execute the bamPlot_turbo.py commands
    if plotBam:
        bam1 = dataDict[name1]['bam']
        bam2 = dataDict[name2]['bam']
        bamString = "%s,%s" % (bam1,bam2)
        nameString = "%s,%s" % (name1,name2)
        colorString = "0,0,0:100,100,100"

        #change dir
        os.chdir(pipelineDir)
    
        if len(gainedGFF) > 0:
            #gained command
            plotTitle = "%s_ONLY_SE" % (name2)
            cmd = 'python bamPlot_turbo.py -g %s -b %s -i %s -o %s -n %s -c %s -t %s -r -y UNIFORM -p MULTIPLE' % (genome,bamString,gffFilename_gained,outputFolder,nameString,colorString,plotTitle)
            os.system(cmd)

            #gained window command
            plotTitle = "%s_ONLY_SE_%sKB_WINDOW" % (name2,window/1000)
            cmd = 'python bamPlot_turbo.py -g %s -b %s -i %s -o %s -n %s -c %s -t %s -r -y UNIFORM -p MULTIPLE' % (genome,bamString,gffFilenameWindow_gained,outputFolder,nameString,colorString,plotTitle)
            os.system(cmd)

        if len(lostGFF) > 0:
            #lost command
            plotTitle = "%s_ONLY_SE" % (name1)
            cmd = 'python bamPlot_turbo.py -g %s -b %s -i %s -o %s -n %s -c %s -t %s -r -y UNIFORM -p MULTIPLE' % (genome,bamString,gffFilename_lost,outputFolder,nameString,colorString,plotTitle)
            os.system(cmd)

            #lost command
            plotTitle = "%s_ONLY_SE_%sKB_WINDOW" % (name1,window/1000)
            cmd = 'python bamPlot_turbo.py -g %s -b %s -i %s -o %s -n %s -c %s -t %s -r -y UNIFORM -p MULTIPLE' % (genome,bamString,gffFilenameWindow_lost,outputFolder,nameString,colorString,plotTitle)
            os.system(cmd)


    return
    

#================================================================================
#===============================MAIN RUN=========================================
#================================================================================


#write the actual script here
def main():



    '''
    main run function
    '''

    from optparse import OptionParser

    usage = "usage: %prog [options] -g [GENOME] -d [DATAFILE] -n [DATA_NAMES] -r [ROSE_FOLDERS] -o [OUTPUT_FOLDER]"
    parser = OptionParser(usage = usage)
    #required flags
    parser.add_option("-g","--genome", dest="genome",nargs = 1, default=None,
                      help = "Enter the genome build (HG18,HG19,MM9,RN4) for the project")
    parser.add_option("-d","--data", dest="data",nargs = 1, default=None,
                      help = "Enter the data file for the project")
    parser.add_option("-r","--rose", dest="rose",nargs = 1, default=None,
                      help = "Enter a comma separated list of rose folder")
    parser.add_option("-o","--output", dest="output",nargs = 1, default=None,
                      help = "Enter the output folder for the project")
    parser.add_option("-n","--names", dest="names",nargs = 1, default=None,
                      help = "Enter a comma separated list of names to go with the datasets")


    #additional options
    parser.add_option("-p","--plot", dest="plot",action = 'store_true', default=False,
                      help = "If flagged, will plot differential regions")
    parser.add_option("-a","--all", dest="all",action = 'store_true', default=False,
                      help = "If flagged, will run analysis for all enhancers and not just supers.")
    parser.add_option("-m","--median", dest="median",action = 'store_true', default=False,
                      help = "If flagged, will use median enhancer scaling")
    parser.add_option("-e","--enhancer-type", dest="enhancer_type",nargs = 1,default='super',
                      help = "specify type of enhancer to analyze: super, stretch, superStretch")

    (options,args) = parser.parse_args()

    print(options)
    print(args)
    
    if options.genome and options.data and options.rose and options.output and options.names:
        genome = string.upper(options.genome)
        dataFile = options.data

        roseFolderString = options.rose
        [roseFolder1,roseFolder2] = roseFolderString.split(',')
        parentFolder = utils.formatFolder(options.output,True)
        

        nameString = options.names
        [name1,name2] =nameString.split(',')

        mergeName = "%s_%s_merged" % (name1,name2)

        #option for median scaling
        medianScale = options.median

        plotBam = options.plot
        if options.all:
            superOnly = False
        else:
            superOnly = True

        if superOnly and plotBam:
            print "Running dynamic enhancer analysis on all super enhancers in %s and %s and plotting output to %s" % (name1,name2,parentFolder)
        if superOnly and not plotBam:
            print "Running dynamic enhancer analysis on all super enhancers in %s and %s and writing output to %s" % (name1,name2,parentFolder)
        if not superOnly and plotBam:
            print "Running dynamic enhancer analysis on all enhancers in %s and %s and plotting output to %s. WARNING: Plotting all differential enhancers could take a while" % (name1,name2,parentFolder)
        if not superOnly and not plotBam:
            print "Running dynamic enhancer analysis on all enhancers in %s and %s and writing output to %s." % (name1,name2,parentFolder)

        #part 1
        print "PART1: analyzing ROSE output from %s and %s" % (name1,name2)
        #start with the all enhancer tables from the initial rose calls

        roseFolder1 = pipeline_dfci.formatFolder(roseFolder1,False)
        roseFolder2 = pipeline_dfci.formatFolder(roseFolder2,False)

        roseDict1 = makeRoseDict(roseFolder1)
        roseDict2 = makeRoseDict(roseFolder2)


        #choosing the type of enhancer to analyze
        enhancerCallType = string.lower(options.enhancer_type)
        if superOnly:
            print("ANALYZING ENHANCER TYPE: %s" % (string.upper(enhancerCallType)))
        superFile1 = roseDict1[enhancerCallType]
        superFile2 = roseDict2[enhancerCallType]

        allFile1 = roseDict1['AllEnhancer']
        allFile2 = roseDict2['AllEnhancer']

        print('\tMERGING ENHANCERS AND CALLING ROSE')
        if superOnly:
            if len(superFile1) ==0:
                print "ERROR: UNABLE TO FIND %s FILES IN %s" % (enhancerCallType,roseFolder1)
                sys.exit()
            if len(superFile2) == 0:
                print "ERROR: UNABLE TO FIND %s FILES IN %s" % (enhancerCallType,roseFolder2)
                sys.exit()
            roseOutput = callMergeSupers(dataFile,superFile1,superFile2,name1,name2,mergeName,genome,parentFolder)

        else:

            roseOutput = callMergeSupers(dataFile,allFile1,allFile2,name1,name2,mergeName,genome,parentFolder)



        print('\tCALCULATING ENHANCER DELTA AND MAKING PLOTS')

        #part2 is the R script
        mergedGFFFile = '%s%s_%s_MERGED_REGIONS_-0_+0.gff' % (parentFolder,string.upper(genome),mergeName)    
        rcmd = callDeltaRScript(mergedGFFFile,parentFolder,dataFile,name1,name2,allFile1,allFile2,medianScale)
        print(rcmd) 
        os.system(rcmd)

        time.sleep(30)
        callRoseGeneMapper(mergedGFFFile,genome,parentFolder,name1)

        #rank the genes


        #part 3
        #rank the delta
        print "PART 3: assinging ranks to differential enhancers"
        print('\tASSIGNING SUPER RANK TO MERGED ENHANCERS')

        gffName = '%s_%s_MERGED_REGIONS_-0_+0' % (string.upper(genome),mergeName)
        enhancerToGeneFile = "%s%s_ROSE/%s_0KB_STITCHED_ENHANCER_DELTA_ENHANCER_TO_GENE_100KB.txt" % (parentFolder,name1,gffName)
        if utils.checkOutput(enhancerToGeneFile):
            rankOutput = "%s%s_ROSE/%s_0KB_STITCHED_ENHANCER_DELTA_ENHANCER_TO_GENE_100KB_RANK.txt" % (parentFolder,name1,gffName)
            assignEnhancerRank(enhancerToGeneFile,allFile1,allFile2,name1,name2,rankOutput)
        else:
            print('ERROR: DELTA SCRIPT OR ROSE GENE MAPPER FAILED TO RUN')
            sys.exit()

        #make the rank plot
        print('MAKING RANK PLOTS')
        if utils.checkOutput(rankOutput):
            rcmd = callRankRScript(rankOutput,name1,name2,superFile1,superFile2)
            print(rcmd)
            os.system(rcmd)
        else:
            print('ERROR: RANK PLOT SCRIPT FAILED TO RUN')
            sys.exit()

        time.sleep(30)

        print('FINISHING OUTPUT')
        finishRankOutput(dataFile,rankOutput,genome,parentFolder,mergeName,name1,name2,1,100000,superOnly,plotBam)
    else:
        parser.print_help()
        sys.exit()

main()
