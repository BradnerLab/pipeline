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


#================================================================================
#============================GLOBAL PARAMETERS===================================
#================================================================================

#add locations of files and global parameters in this section

pipelineDir = '/ark/home/cl512/pipeline/'
roseDir = '/ark/home/cl512/rose/'

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

#write your specific functions here


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
                
                signal = float(line[6]) - float(line[7])
                rank = int(line[-2])
                enhancerID = name+'_'+line[0]
                seDict[enhancerID] = {'rank':rank,'signal':signal}

        else:

            signal = float(line[6]) - float(line[7])
            rank = int(line[-2])
            enhancerID = name+'_'+line[0]
            seDict[enhancerID] = {'rank':rank,'signal':signal}

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

    namesList = [name1]    
    extraMap = [name2,dataDict[name2]['background']]


    return pipeline_dfci.callRose(dataFile,'',parentFolder,namesList,extraMap,mergedGFFFile,tss=0,stitch=0)


def callMergeSupers(dataFile,superFile1,superFile2,name1,name2,mergedGFFFile,parentFolder):

    '''
    this is the main run function for the script
    all of the work should occur here, but no functions should be defined here
    '''
    
    mergedGFF = mergeCollections(superFile1,superFile2,name1,name2,mergedGFFFile)

    #call rose on the merged shit    


    roseBashFile = callRoseMerged(dataFile,mergedGFF,name1,name2,parentFolder)
    print('i can has rose bash file %s' % (roseBashFile))

    #run the bash command
    os.system('bash %s' % (roseBashFile))

def callDeltaRScript(mergedGFFFile,parentFolder,name1,name2):

    '''
    runs the R script
    '''

    gffName = mergedGFFFile.split('/')[-1].split('.')[0]
    stitchedFile = "%s%s_ROSE/%s_0KB_STITCHED_ENHANCER_REGION_MAP.txt" % (parentFolder,name1,gffName)
    #print(stitchedFile)
    os.chdir(pipelineDir)

    rcmd = "R --no-save %s %s %s < ./dynamicEnhancer_plot.R" % (stitchedFile,name1,name2)

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
    
    os.chdir(roseDir)
    cmd = 'python ROSE_geneMapper.py -g %s -i %s -w 100000' % (genome,deltaFile)
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

    usage = "usage: %prog [options] -g [GENOME] -d [DATAFILE] -r [ROSE_FOLDERS] -o [OUTPUT_FOLDER]"
    parser = OptionParser(usage = usage)
    #required flags
    parser.add_option("-g","--genome", dest="genome",nargs = 1, default=None,
                      help = "Enter the genome build (HG18,HG19,MM9) for the project")
    parser.add_option("-d","--data", dest="data",nargs = 1, default=None,
                      help = "Enter the data file for the project")
    parser.add_option("-r","--rose", dest="rose",nargs = 1, default=None,
                      help = "Enter a comma separated list of rose folder")
    parser.add_option("-o","--output", dest="output",nargs = 1, default=None,
                      help = "Enter the output folder for the project")

    #additional options
    parser.add_option("-n","--names", dest="names",nargs = 1, default=None,
                      help = "Enter a comma separated list of names to go with the datasets")
    parser.add_option("-p","--plot", dest="plot",action = 'store_true', default=False,
                      help = "If flagged, will plot differential regions")
    parser.add_option("-a","--all", dest="all",action = 'store_true', default=False,
                      help = "If flagged, will run analysis for all enhancers and not just supers.")

    (options,args) = parser.parse_args()

    print(options)
    print(args)
    
    if options.genome and options.data and options.rose and options.output:
        genome = string.upper(options.genome)
        dataFile = options.data

        roseFolderString = options.rose
        [roseFolder1,roseFolder2] = roseFolderString.split(',')
        parentFolder = utils.formatFolder(options.output,True)
        
        if options.names:
            nameString = options.names
            [name1,name2] =nameString.split(',')
        else:
            name1 = roseFolder1.split('/')[-1]
            name1 = string.replace(name1,'_ROSE','')

            name2 = roseFolder2.split('/')[-1]
            name2 = string.replace(name2,'_ROSE','')

        mergeName = "%s_%s_merged" % (name1,name2)

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
        superFile1 = '%s%s_peaks_SuperEnhancers.table.txt' % (roseFolder1,name1)
        superFile2 = '%s%s_peaks_SuperEnhancers.table.txt' % (roseFolder2,name2)

        allFile1 = '%s/%s_peaks_AllEnhancers.table.txt' % (roseFolder1,name1)
        allFile2 = '%s/%s_peaks_AllEnhancers.table.txt' % (roseFolder2,name2)

        print('\tMERGING ENHANCERS AND CALLING ROSE')
        if superOnly:
            mergedGFFFile = '%s%s_%s_MERGED_SUPERS_-0_+0.gff' % (parentFolder,string.upper(genome),mergeName)
            callMergeSupers(dataFile,superFile1,superFile2,name1,name2,mergedGFFFile,parentFolder)

        else:
            mergedGFFFile = '%s%s_%s_MERGED_ENHANCERS_-0_+0.gff' % (parentFolder,string.upper(genome),mergeName)
            callMergeSupers(dataFile,allFile1,allFile2,name1,name2,mergedGFFFile,parentFolder)


        if superOnly:
            superOutput = "%s%s_ROSE/%s_%s_MERGED_SUPERS_-0_+0_SuperEnhancers_ENHANCER_TO_GENE.txt" % (parentFolder,name1,string.upper(genome),mergeName)
        else:
            superOutput = "%s%s_ROSE/%s_%s_MERGED_ENHANCERS_-0_+0_SuperEnhancers_ENHANCER_TO_GENE.txt" % (parentFolder,name1,string.upper(genome),mergeName)

        print('\tCALCULATING ENHANCER DELTA AND MAKING PLOTS')
        if utils.checkOutput(superOutput):
            #part2 is the R script
            rcmd = callDeltaRScript(mergedGFFFile,parentFolder,name1,name2)
            print(rcmd) 
            os.system(rcmd)
            time.sleep(30)
            callRoseGeneMapper(mergedGFFFile,genome,parentFolder,name1)
        else:
            print('ERROR: ROSE CALL FAILED')
            sys.exit()

        #rank the genes


        #part 3
        #rank the delta
        print "PART 3: assinging ranks to differential enhancers"
        print('\tASSIGNING SUPER RANK TO MERGED ENHANCERS')
        if superOnly:
            gffName = '%s_%s_MERGED_SUPERS_-0_+0' % (string.upper(genome),mergeName)
        else:
            gffName = '%s_%s_MERGED_ENHANCERS_-0_+0' % (string.upper(genome),mergeName)
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

# #=====================================================
# #=======================MAC_SUPERS====================
# #=====================================================
# #MAC H4K12AC supers
# mergeName = 'MAC_MERGED_SUPERS'
# genome ='mm9'
# dataFile = '/ark/home/cl512/projects/mouse_macrophage/MAC_TABLE.txt'
# mergeFolder = '/ark/home/cl512/projects/mouse_macrophage/mergeTest/%s/' % (mergeName)
# mergeFolder = pipeline_dfci.formatFolder(mergeFolder,True)
# roseFolder1 = '/ark/home/cl512/projects/mouse_macrophage/rose/MAC_H4K12AC_0H_MINUS_ROSE'
# roseFolder2 = '/ark/home/cl512/projects/mouse_macrophage/rose/MAC_H4K12AC_1H_MINUS_ROSE'

# name1 = 'MAC_H4K12AC_0H_MINUS'
# name2 = 'MAC_H4K12AC_1H_MINUS'


# main(dataFile,genome,mergeFolder,roseFolder1,roseFolder2,name1,name2,mergeName,True,True)



#=====================================================
#===================MYC NB LEVELS=====================
#=====================================================
#MYC vs MYCN in BE(2)C
mergeName = 'BE2C_MYC_SUPERS'
genome ='hg18'
dataFile = '/ark/home/cl512/projects/neuroblastoma/NEURO_TABLE.txt'
mergeFolder = '/ark/home/cl512/projects/neuroblastoma/MYC_analysis/dynamicEnhancer/%s/' % (mergeName)
mergeFolder = pipeline_dfci.formatFolder(mergeFolder,True)
roseFolder1 = '/ark/home/cl512/projects/neuroblastoma/MYC_analysis/be2c_mycn_rose/'
roseFolder2 = '/ark/home/cl512/projects/neuroblastoma/MYC_analysis/be2c_myc_rose/'

name1 = 'BE2C_MYCN'
name2 = 'BE2C_MYC'


#main(dataFile,genome,mergeFolder,roseFolder1,roseFolder2,name1,name2,mergeName,True,True)




# #=====================================================
# #=======================786-O_SUPERS==================
# #=====================================================
# #786O
# mergeName = 'RCC_786O_MERGED_SUPERS'
# genome ='hg18'
# dataFile = '/ark/home/cl512/projects/renal/RENAL_TABLE.txt'
# mergeFolder = '/ark/home/cl512/projects/renal/dynamicEnhancer/%s/' % (mergeName)
# mergeFolder = pipeline_dfci.formatFolder(mergeFolder,True)
# roseFolder1 = '/ark/home/cl512/projects/renal/rose/786-O_H3K27AC_ROSE'
# roseFolder2 = '/ark/home/cl512/projects/renal/rose/28-6_H3K27AC_ROSE'

# name1 = '786-O_H3K27AC'
# name2 = '28-6_H3K27AC'


# main(dataFile,genome,mergeFolder,roseFolder1,roseFolder2,name1,name2,mergeName,True,True)



# #=====================================================
# #=======================786-O_TYPICALS================
# #=====================================================
# #786O
# mergeName = 'RCC_786O_MERGED_ENHANCERS'
# genome ='hg18'
# dataFile = '/ark/home/cl512/projects/renal/RENAL_TABLE.txt'
# mergeFolder = '/ark/home/cl512/projects/renal/dynamicEnhancer/%s/' % (mergeName)
# mergeFolder = pipeline_dfci.formatFolder(mergeFolder,True)
# roseFolder1 = '/ark/home/cl512/projects/renal/rose/786-O_H3K27AC_ROSE'
# roseFolder2 = '/ark/home/cl512/projects/renal/rose/28-6_H3K27AC_ROSE'

# name1 = '786-O_H3K27AC'
# name2 = '28-6_H3K27AC'


# main(dataFile,genome,mergeFolder,roseFolder1,roseFolder2,name1,name2,mergeName,False,False)


# #=====================================================
# #===================EC_ALL_ENHANCERS==================
# #=====================================================
# #EC_BRD4 ALL enhancers in con and tnf
# mergeName = 'EC_MERGED_ENHANCERS'
# mergeFolder = '/ark/home/cl512/projects/athero/mergeTest/%s/' % (mergeName)
# mergeFolder = pipeline_dfci.formatFolder(mergeFolder,True)
# roseFolder1 = '/ark/home/cl512/projects/athero/rose/EC_BRD4_CON_ROSE/'
# roseFolder2 = '/ark/home/cl512/projects/athero/rose/EC_BRD4_TNF_ROSE/'

# name1 = 'EC_BRD4_CON'
# name2 = 'EC_BRD4_TNF'


# main(dataFile,genome,mergeFolder,roseFolder1,roseFolder2,name1,name2,mergeName,False,False)



# #=====================================================
# #===================EC_SUPER_ENHANCERS==================
# #=====================================================
# #EC_BRD4 ALL enhancers in con and tnf
# mergeName = 'EC_MERGED_SUPERS'
# mergeFolder = '/ark/home/cl512/projects/athero/mergeTest/%s/' % (mergeName)
# mergeFolder = pipeline_dfci.formatFolder(mergeFolder,True)
# roseFolder1 = '/ark/home/cl512/projects/athero/rose/EC_BRD4_CON_ROSE/'
# roseFolder2 = '/ark/home/cl512/projects/athero/rose/EC_BRD4_TNF_ROSE/'

# name1 = 'EC_BRD4_CON'
# name2 = 'EC_BRD4_TNF'


# main(dataFile,genome,mergeFolder,roseFolder1,roseFolder2,name1,name2,mergeName,True,True)




# #=====================================================
# #====================786O_H3K27AC=====================
# #=====================================================
# #786O_H3K27AC
# mergeName = '786O_H3K27AC'
# mergeFolder = '/home/clin/projects/131106_seComp/mergeAnalysis/%s/' % (mergeName)
# mergeFolder = pipeline_dfci.formatFolder(mergeFolder,True)
# roseFolder1 = '/home/clin/projects/131106_seComp/rose/786O/786O_H3K27AC_VHL_WT_ROSE/'
# roseFolder2 = '/home/clin/projects/131106_seComp/rose/786O/786O_H3K27AC_VHL_NULL_ROSE/'
# name1 = '786O_H3K27AC_VHL_WT'
# name2 = '786O_H3K27AC_VHL_NULL'


# main(dataFile,genome,mergeFolder,roseFolder1,roseFolder2,name1,name2,mergeName)

# # #=====================================================
# # #====================UMRC2_H3K27AC=====================
# # #=====================================================
# # #UMRC2_H3K27AC
# # mergeName = 'UMRC2_H3K27AC'
# # mergeFolder = '/home/clin/projects/131106_seComp/mergeAnalysis/%s/' % (mergeName)
# # mergeFolder = pipeline_dfci.formatFolder(mergeFolder,True)
# # roseFolder1 = '/home/clin/projects/131106_seComp/rose/UMRC2/UMRC2_H3K27AC_VHL_WT_ROSE/'
# # roseFolder2 = '/home/clin/projects/131106_seComp/rose/UMRC2/UMRC2_H3K27AC_VHL_NULL_ROSE/'
# # name1 = 'UMRC2_H3K27AC_VHL_WT'
# # name2 = 'UMRC2_H3K27AC_VHL_NULL'


# # main(dataFile,genome,mergeFolder,roseFolder1,roseFolder2,name1,name2,mergeName)


# #=====================================================
# #====================RCC4_H3K27AC=====================
# #=====================================================
# #RCC4_H3K27AC
# mergeName = 'RCC4_H3K27AC'
# mergeFolder = '/home/clin/projects/131106_seComp/mergeAnalysis/%s/' % (mergeName)
# mergeFolder = pipeline_dfci.formatFolder(mergeFolder,True)
# roseFolder1 = '/home/clin/projects/131106_seComp/rose/RCC4/RCC4_H3K27AC_VHL_WT_ROSE/'
# roseFolder2 = '/home/clin/projects/131106_seComp/rose/RCC4/RCC4_H3K27AC_VHL_NULL_ROSE/'
# name1 = 'RCC4_H3K27AC_VHL_WT'
# name2 = 'RCC4_H3K27AC_VHL_NULL'


# main(dataFile,genome,mergeFolder,roseFolder1,roseFolder2,name1,name2,mergeName)


# #=====================================================
# #====================RCC4_MED1=====================
# #=====================================================
# #RCC4_MED1
# mergeName = 'RCC4_MED1'
# mergeFolder = '/home/clin/projects/131106_seComp/mergeAnalysis/%s/' % (mergeName)
# mergeFolder = pipeline_dfci.formatFolder(mergeFolder,True)
# roseFolder1 = '/home/clin/projects/131106_seComp/rose/RCC4/RCC4_MED1_VHL_WT_ROSE/'
# roseFolder2 = '/home/clin/projects/131106_seComp/rose/RCC4/RCC4_MED1_VHL_NULL_ROSE/'
# name1 = 'RCC4_MED1_VHL_WT'
# name2 = 'RCC4_MED1_VHL_NULL'


# main(dataFile,genome,mergeFolder,roseFolder1,roseFolder2,name1,name2,mergeName)




# #=====================================================
# #====================A704_H3K27AC=====================
# #=====================================================
# #A704_H3K27AC
# mergeName = 'A704_H3K27AC'
# mergeFolder = '/home/clin/projects/131106_seComp/mergeAnalysis/%s/' % (mergeName)
# mergeFolder = pipeline_dfci.formatFolder(mergeFolder,True)
# roseFolder1 = '/home/clin/projects/131106_seComp/rose/A704/A704_H3K27AC_BAF180_WT_ROSE/'
# roseFolder2 = '/home/clin/projects/131106_seComp/rose/A704/A704_H3K27AC_BAF180_NULL_ROSE/'
# name1 = 'A704_H3K27AC_BAF180_WT'
# name2 = 'A704_H3K27AC_BAF180_NULL'


# main(dataFile,genome,mergeFolder,roseFolder1,roseFolder2,name1,name2,mergeName)


# #=====================================================
# #====================CAKI2_H3K27AC=====================
# #=====================================================
# #CAKI2_H3K27AC
# mergeName = 'CAKI2_H3K27AC'
# mergeFolder = '/home/clin/projects/131106_seComp/mergeAnalysis/%s/' % (mergeName)
# mergeFolder = pipeline_dfci.formatFolder(mergeFolder,True)
# roseFolder1 = '/home/clin/projects/131106_seComp/rose/CAKI2/CAKI2_H3K27AC_BAF180_WT_ROSE/'
# roseFolder2 = '/home/clin/projects/131106_seComp/rose/CAKI2/CAKI2_H3K27AC_BAF180_NULL_ROSE/'
# name1 = 'CAKI2_H3K27AC_BAF180_WT'
# name2 = 'CAKI2_H3K27AC_BAF180_NULL'


# main(dataFile,genome,mergeFolder,roseFolder1,roseFolder2,name1,name2,mergeName)


# #=====================================================
# #====================CAKI2_MED1=====================
# #=====================================================
# #CAKI2_MED1
# mergeName = 'CAKI2_MED1'
# mergeFolder = '/home/clin/projects/131106_seComp/mergeAnalysis/%s/' % (mergeName)
# mergeFolder = pipeline_dfci.formatFolder(mergeFolder,True)
# roseFolder1 = '/home/clin/projects/131106_seComp/rose/CAKI2/CAKI2_MED1_BAF180_WT_ROSE/'
# roseFolder2 = '/home/clin/projects/131106_seComp/rose/CAKI2/CAKI2_MED1_BAF180_NULL_ROSE/'
# name1 = 'CAKI2_MED1_BAF180_WT'
# name2 = 'CAKI2_MED1_BAF180_NULL'


# main(dataFile,genome,mergeFolder,roseFolder1,roseFolder2,name1,name2,mergeName)
