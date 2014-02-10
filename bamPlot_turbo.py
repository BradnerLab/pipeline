#bamPlot.py

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

import sys


#MUST SET THIS STRING TO THE PATH OF THE bamliquidator program
bamliquidatorString = 'bamliquidator'

#as of now the number of bins to sample the space is hard wired
nBins = 200

from utils import *
import pipeline_dfci
import time
import subprocess
import os
import string


#script that takes in a list of bams, makes an intermediate table file, and then calls R to make the plot

#updated for batching with bamliquidator magic


#==========================================================================
#======================HELPER FUNCTIONS====================================
#==========================================================================




def loadAnnotFile(genome):

    '''
    load in the annotation and create a geneDict and transcription collection
    '''
    genomeDict = {
        'HG18':'./annotation/hg18_refseq.ucsc',
        'MM9': './annotation/mm9_refseq.ucsc',
        'hg18':'./annotation/hg18_refseq.ucsc',
        'mm9': './annotation/mm9_refseq.ucsc',
        'HG19':'./annotation/hg19_refseq.ucsc',
        'hg19':'./annotation/hg19_refseq.ucsc'
        }

    annotFile = genomeDict[genome]
    #geneList =['NM_002460','NM_020185']
    geneList = []
    geneDict = makeGenes(annotFile,geneList,True)
    txCollection =makeTranscriptCollection(annotFile,0,0,500,geneList) 

    return geneDict,txCollection

def tasteTheRainbow(n):
    
    '''
    samples rainbow color space
    '''
    from colorsys import hsv_to_rgb
    
    colorList = []
    nRange = [x/float(n) for x in range(0,n)]
    for i in nRange:
        color = [int(255 * x) for x in list(hsv_to_rgb(i,.9,.9))]
    
        colorList.append(color)
    return colorList

def mapGFFLineToAnnot(gffLine,outFolder,nBins,geneDict,txCollection,sense='both',header=''):

    '''
    for every line produces a file with all of the rectangles to draw
    '''

    if len(header) == 0:
        gffString = '%s_%s_%s_%s' % (gffLine[0],gffLine[6],gffLine[3],gffLine[4])
    else:
        gffString = header
    diagramTable = [[0,0,0,0]]
    nameTable = [['',0,0]]
    gffLocus = Locus(gffLine[0],int(gffLine[3]),int(gffLine[4]),gffLine[6],gffLine[1])    

    scaleFactor = float(nBins)/gffLocus.len()
    #plotting buffer for diagrams
    plotBuffer = int(gffLocus.len()/float(nBins)*20)

    overlapLoci = txCollection.getOverlap(gffLocus,sense='both')
    geneList = [locus.ID() for locus in overlapLoci]
    
    if gffLine[6] == '-':
        refPoint = int(gffLine[4])
    else:
        refPoint = int(gffLine[3])
    offsetCollection=LocusCollection([],500)
    for geneID in geneList:

        gene = geneDict[geneID]

        print(gene.commonName())
        if len(gene.commonName())>1:
            name = gene.commonName()
        else:
            name = geneID
        offset = 4*len(offsetCollection.getOverlap(gene.txLocus()))
        offsetCollection.append(makeSearchLocus(gene.txLocus(),plotBuffer,plotBuffer))
        #write the name of the gene down
        if gene.sense()=='+':
            geneStart = gene.txLocus().start()
        else:
            geneStart = gene.txLocus().end()
        geneStart = abs(geneStart-refPoint)*scaleFactor
        nameTable.append([name,geneStart,-2-offset])
        #draw a line across the entire txLocus

        [start,stop] = [abs(x-refPoint)*scaleFactor for x in gene.txLocus().coords()]
        diagramTable.append([start,-0.01-offset,stop,0.01-offset])

        
        #now draw thin boxes for all txExons
        if len(gene.txExons()) > 0:
            for txExon in gene.txExons():

                [start,stop] = [abs(x-refPoint)*scaleFactor for x in txExon.coords()]

                diagramTable.append([start,-0.5-offset,stop,0.5-offset])

        #now draw fatty boxes for the coding exons if any
        if len(gene.cdExons()) > 0:
            for cdExon in gene.cdExons():
        
                [start,stop] = [abs(x-refPoint)*scaleFactor for x in cdExon.coords()]

                diagramTable.append([start,-1-offset,stop,1-offset])

    unParseTable(diagramTable,outFolder+gffString+'_diagramTemp.txt','\t')
    unParseTable(nameTable,outFolder+gffString+'_nameTemp.txt','\t')    

def mapBamToGFFLine(bamFile,MMR,name,gffLine,color,nBins,sense = 'both',extension = 200):
    '''maps reads from a bam to a gff'''


    
    print('using a MMR value of %s' % (MMR))

    
    line = gffLine[0:9]
    gffLocus = Locus(line[0],int(line[3]),int(line[4]),line[6],line[1])

    #setting up the output clusterline
    colorLine = color
    bamName = bamFile.split('/')[-1]
    clusterLine = [bamName,gffLocus.ID(),name,gffLocus.__str__()]+colorLine


    binSize = gffLocus.len()/nBins
    #some regions will be too short to get info on
    #we just kick these back and abandon them
    if binSize == 0:
        clusterLine+=['NA']*int(nBins)
        return clusterLine


    #flippy flip if sense is negative
    senseTrans = maketrans('-+.','+-+')
    if sense == '-':
        bamSense = string.translate(gffLocus.sense(),senseTrans)
    elif sense == '+':
        bamSense = gffLocus.sense()
    else:
        bamSense = '.'
    #using the bamLiquidator to get the readstring            
    #print('using nBin of %s' % nBin)

    bamCommand = "%s %s %s %s %s %s %s %s" % (bamliquidatorString,bamFile,gffLocus.chr(),gffLocus.start(),gffLocus.end(),bamSense,nBins,extension)
    #print(bamCommand)
    getReads = subprocess.Popen(bamCommand,stdin = subprocess.PIPE,stderr = subprocess.PIPE,stdout = subprocess.PIPE,shell = True)
    readString = getReads.communicate()
    denList = readString[0].split('\n')[:-1]

    #flip the denList if the actual gff region is -
    if gffLocus.sense() == '-':
        denList = denList[::-1]

    #converting from units of total bp of read sequence per bin to rpm/bp
    denList = [round(float(x)/binSize/MMR,4) for x in denList]

    clusterLine += denList
        
    return clusterLine


def callRPlot(summaryFile,outFile,yScale,plotStyle):

    '''
    calls the R plotting thingy
    '''

    cmd = 'R --no-save %s %s %s %s < ./bamPlot_turbo.R' % (summaryFile,outFile,yScale,plotStyle)
    print('calling command %s' % (cmd))
    return cmd

    




def makeBamPlotTables(gff,genome,bamFileList,colorList,nBins,sense,extension,rpm,outFolder,names,title):

    '''
    makes a plot table for each line of the gff mapped against all the bams in the bamList
    '''
    
    #load in the gff
    if type(gff) == str:
        gff = parseTable(gff,'\t')

    #load in the annotation
    print('loading in annotation for %s' % (genome))
    geneDict,txCollection = loadAnnotFile(genome)

    #make an MMR dict so MMRs are only computed once
    print('Getting information about read depth in bams')
    mmrDict = {}
    for bamFile in bamFileList:
        #millionMappedReads
        bam = Bam(bamFile)
        if rpm:    
            MMR= round(float(bam.getTotalReads('mapped'))/1000000,4)
        else:
            MMR = 1
        mmrDict[bamFile] = MMR
        #mmrDict[bamFile] = 21.5377
        
    ticker = 1
    #go line by line in the gff
    summaryTable = [['DIAGRAM_TABLE','NAME_TABLE','PLOT_TABLE','CHROM','ID','SENSE','START','END']]
    for gffLine in gff:
        gffString = 'line_%s_%s_%s_%s_%s_%s' % (ticker,gffLine[0],gffLine[1],gffLine[6],gffLine[3],gffLine[4])
        ticker+=1
        print('writing the gene diagram table for region %s' % (gffLine[1]))
        mapGFFLineToAnnot(gffLine,outFolder,nBins,geneDict,txCollection,sense='both',header=gffString)

        outTable = []

        outTable.append(['BAM','GENE_ID','NAME','LOCUSLINE','COLOR1','COLOR2','COLOR3'] + ['bin_'+str(n) for n in range(1,int(nBins)+1,1)])
        
        for i in range(0,len(bamFileList),1):
            bamFile = bamFileList[i]
            name = names[i]
            
            color = colorList[i]
            print('getting data for location %s in dataset %s' % (gffLine[1],bamFile))
            mmr = mmrDict[bamFile]
            newLine = mapBamToGFFLine(bamFile,mmr,name,gffLine,color,nBins,sense,extension)


            outTable.append(newLine)

        #get the gene name
        if geneDict.has_key(gffLine[1]):
            geneName = geneDict[gffLine[1]].commonName()
        else:
            geneName = gffLine[1]
        unParseTable(outTable,outFolder+gffString+'_plotTemp.txt','\t')
        diagramTable = outFolder+gffString+'_diagramTemp.txt'
        plotTable = outFolder+gffString+'_plotTemp.txt'
        nameTable = outFolder +gffString+'_nameTemp.txt'
        summaryTable.append([diagramTable,nameTable,plotTable,gffLine[0],geneName,gffLine[6],gffLine[3],gffLine[4]])
    summaryTableFileName = "%s%s_summary.txt" % (outFolder,title)
    unParseTable(summaryTable,summaryTableFileName,'\t')
    return summaryTableFileName






def main():
    
    '''
    main run function
    '''

    from optparse import OptionParser

    usage = "usage: %prog [options] -g [GENOME] -b [SORTED BAMFILE(S)] -i [INPUTFILE] -o [OUTPUTFOLDER]"
    parser = OptionParser(usage = usage)
    #required flags
    parser.add_option("-b","--bam", dest="bam",nargs = 1, default=None,
                      help = "Enter a comma separated list of .bam files to be processed.")
    parser.add_option("-i","--input", dest="input",nargs = 1, default=None,
                      help = "Enter .gff or genomic region e.g. chr1:+:1-1000.")
    parser.add_option("-g","--genome",dest="genome",nargs =1, default = None,
                      help = "specify a genome, options are hg18 or mm9 right now")


    #output flag
    parser.add_option("-o","--output", dest="output",nargs = 1, default=None,
                      help = "Enter the output folder.")
    #additional options
    parser.add_option("-c","--color", dest="color",nargs = 1, default=None,
                      help = "Enter a colon separated list of colors e.g. 255,0,0:255,125,0, default samples the rainbow")
    parser.add_option("-s","--sense", dest="sense",nargs = 1, default='both',
                      help = "Map to '+','-' or 'both' strands. Default maps to both.")
    parser.add_option("-e","--extension", dest="extension",nargs = 1, default=200,
                      help = "Extends reads by n bp. Default value is 200bp")
    parser.add_option("-r","--rpm", dest="rpm",action = 'store_true', default=False,
                      help = "Normalizes density to reads per million (rpm) Default is True")
    parser.add_option("-y","--yScale",dest="yScale",nargs =1, default = "relative",
                      help = "Choose either relative or uniform y axis scaling. options = 'relative,uniform' Default is relative scaling")
    parser.add_option("-n","--names",dest="names",nargs =1, default = None,
                      help = "Enter a comma separated list of names for your bams")
    parser.add_option("-p","--plot",dest="plot",nargs =1, default = "multiple",
                      help = "Choose either all lines on a single plot or multiple plots. options = 'single,multiple'")
    parser.add_option("-t","--title",dest ="title",nargs=1,default = '',
                      help = "Specify a title for the output plot(s), default will be the coordinate region")
                  


    (options,args) = parser.parse_args()

    print(options)
    print(args)
    
    if options.bam and options.input and options.genome and options.output:

        #bring in the bams
        bamFileList = options.bam.split(',')
        
        #bring in the gff
        try:
            gff = parseTable(options.input,'\t')
            gffName = options.input.split('/')[-1].split('.')[0]
        except IOError:
            #means a coordinate line has been given e.g. chr1:+:1-100

            chromLine = options.input.split(':')
            
            chrom = chromLine[0]
            sense = chromLine[1]
            [start,end] = chromLine[2].split('-')
            if chrom[0:3] != 'chr':
                print('ERROR: UNRECOGNIZED GFF OR CHROMOSOME LINE INPUT')
                exit()
            gffLine = [chrom,'',options.input,start,end,'',sense,'','']
            gffName = "%s_%s_%s_%s" % (chrom,sense,start,end)
            gff = [gffLine]

        #bring in the genome
        genome = upper(options.genome)
        if ['HG18','HG19','MM9','RN5'].count(genome) == 0:
            print('ERROR: UNSUPPORTED GENOME TYPE %s. USE HG19,HG18, RN5, OR MM9' % (genome))
            parser.print_help()
            exit()

        #bring in the rest of the options
        
        #output
        rootFolder = options.output
        if rootFolder[-1] != '/':
            rootFolder+='/'
        try:
            foo = os.listdir(rootFolder)
        except OSError:
            print('ERROR: UNABLE TO FIND OUTPUT DIRECTORY %S' % (rootFolder))
            exit()

        #Get analysis title
        if len(options.title) == 0:
            title = gffName
        else:
            title = options.title

        #make a temp folder
        tempFolder = rootFolder + title + '/'
        print("CREATING TEMP FOLDER %s" % (tempFolder))
        pipeline_dfci.formatFolder(tempFolder,create=True)
                         
        #colors
        if options.color:
            colorList = options.color.split(':')
            colorList = [x.split(',') for x in colorList]
            if len(colorList) < len(bamFileList):
                print('WARNING: FEWER COLORS THAN BAMS SPECIFIED. COLORS WILL BE RECYCLED')
                #recycling the color list
                colorList += colorList*(len(bamFileList)/len(colorList))
                colorList = colorList[0:len(bamFileList)]

        else:
            #cycles through the colors of the rainbow
            colorList = tasteTheRainbow(len(bamFileList))

        #sense
        sense = options.sense
        
        extension = int(options.extension)

        rpm = options.rpm

        yScale = upper(options.yScale)
                                     
        #names
        if options.names:
            names = options.names.split(',')
        
            if len(names) != len(bamFileList):
                print('ERROR: NUMBER OF NAMES AND NUMBER OF BAMS DO NOT CORRESPOND')
                parser.print_help()
                exit()
        else:
            names = [x.split('/')[-1] for x in bamFileList]

        #plot style
        plotStyle = upper(options.plot)
        if ['SINGLE','MULTIPLE'].count(plotStyle) == 0:
            print('ERROR: PLOT STYLE %s NOT AN OPTION' % (plotStyle))
            parser.print_help()
            exit()


        #now run!
        summaryTableFileName = makeBamPlotTables(gff,genome,bamFileList,colorList,nBins,sense,extension,rpm,tempFolder,names,title)
        print ("%s is the summary table" % (summaryTableFileName))


        outFile = "%s%s_plots.pdf" % (rootFolder,title)
        rCmd = callRPlot(summaryTableFileName,outFile,yScale,plotStyle)

        #open a bash file to get shit done
        bashFileName = "%s%s_Rcmd.sh" % (tempFolder,title)
        bashFile = open(bashFileName,'w')
        bashFile.write(rCmd)
        bashFile.close()
        print("Wrote R command to %s" % (bashFileName))
        os.system("bash %s" % (bashFileName))



        
    else:
        parser.print_help()
        exit()

if __name__ == "__main__":
    main()

