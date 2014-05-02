#!/usr/bin/python
#bamToGFF.py

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


#script to grab reads from a bam that align to a .gff file

from utils import *

from collections import defaultdict

import os
import string

def parseSamHeader(samFile):
    '''parses any sam type file with a 3 column tab del header'''
    samDict = {}
    sam = open(samFile,'r')

    for line in sam:
        if line[0] == '@':
            headerLine = line[:-1].split('\t')
            samDict[headerLine[1]] = headerLine[2]
        else:
            break
    sam.close()

    return samDict



#THIS FUNCTION USED TO TRY TO LOOK UP THE UNIQUE NUMBER OF READS OR COMPUTE IT BUT WAS LATER EDITED
#BRIAN BEING THAT GUY, YOU CAN FIX IT UP AND MAKE IT FIND YOUR UNIQUE MMR
def getUniquelyMappingReads(bamFile):

    '''
    scripts designed to return the total number of uniquely mapping sequence tags from a bam file,
    first by looking for a corresponding stats file in the folder, and second by manually computing the number.
    '''
    
    #a uniquely mapping sequence tag is defined by collapsing all tags that map to the exact same location with the exact same sequence

    #first try to extract the number from a stats file
    fullPath = os.path.abspath(bamFile)
    bamName = fullPath.split('/')[-1].split('.')[0]
    pathFolder = join(fullPath.split('/')[0:-1],'/')
    bamFiles = os.listdir(pathFolder)
    statsFile = filter(lambda x: x.count(bamName) ==1 and x.count('stats.') ==1,bamFiles)
    if len(statsFile) == 1:
        print('USING STATS FILE %s' % (pathFolder+'/'+statsFile[0]))
        samDict = parseSamHeader(pathFolder+'/'+statsFile[0])
        return int(samDict['UniquelyMappingSequenceTags'])
    elif len(statsFile) > 1:
        statsFile = filter(lambda x: x.count(bamName) ==1 and x.count('stats.concise') ==1,bamFiles)
        print('USING STATS FILE %s' % (pathFolder+'/'+statsFile[0]))
        samDict = parseSamHeader(pathFolder+'/'+statsFile[0])
        return int(samDict['UniquelyMappingSequenceTags'])

    else:
        print('no precomputed stats file found for %s.' % (bamFile))
        return None
        
        



def mapBamToGFF(bamFile,gff,sense = 'both',unique = 0,extension = 200,floor = 0,density = False,rpm = False,binSize = 25,clusterGram = None,matrix = None,raw = False,includeJxnReads = False):
    '''maps reads from a bam to a gff'''
    floor = int(floor)
    bam = Bam(bamFile)
    #if cluster is specified, override certain flags
    if clusterGram:
        density =True
        binSize = int(binSize)
    else:
        binSize = 25

    if matrix:
        density = True
    #new GFF to write to
    newGFF = []
    #millionMappedReads


    #SECTION CHANGED FOR BRIAN BEING THAT GUY
    if float(unique) == 0.0:
        unique = False
        if rpm:    
            MMR= round(float(bam.getTotalReads('mapped'))/1000000,4)
        else:
            MMR = 1
    else:

        if rpm:
            MMR = float(unique)
        else:
            MMR = 1
        unique = True



    print('using a MMR value of %s' % (MMR))
    
    senseTrans = maketrans('-+.','+-+')

    
    if type(gff) == str:
        gff = parseTable(gff,'\t')

    #setting up a clustergram table
    if clusterGram:
        #first grab a header line
        line = gff[0]
        gffLocus = Locus(line[0],int(line[3]),int(line[4]),line[6],line[1])
        
        nBins = gffLocus.len()/binSize
        binSizeList = [nBins]

        #now go through each line of the gff and make sure they're all the same length
        for i in range(0,len(gff),1):
            line = gff[i]
            gffLocus = Locus(line[0],int(line[3]),int(line[4]),line[6],line[1])
            binSizeList.append(gffLocus.len()/binSize)
        binSizeList = uniquify(binSizeList)
        if len(binSizeList) > 1: 
            print('WARNING: lines in gff are of different length. Output clustergram will have variable row length')
        newGFF.append(['GENE_ID','locusLine'] + [str(x*binSize)+'_'+bamFile.split('/')[-1] for x in range(1,max(binSizeList)+1,1)])        
        
    #setting up a maxtrix table
    if matrix:
        newGFF.append(['GENE_ID','locusLine'] + ['bin_'+str(n)+'_'+bamFile.split('/')[-1] for n in range(1,int(matrix)+1,1)])        

    #getting and processing reads for gff lines
    ticker = 0
    print('Number lines processed')
    for line in gff:
        line = line[0:9]
        if ticker%100 == 0:
            print ticker
        ticker+=1
        gffLocus = Locus(line[0],int(line[3]),int(line[4]),line[6],line[1])
        searchLocus = makeSearchLocus(gffLocus,int(extension),int(extension))
        
        reads = bam.getReadsLocus(searchLocus,'both',unique,'none',includeJxnReads)
        #now extend the reads and make a list of extended reads
        extendedReads = []
        for locus in reads:
            if locus.sense() == '+' or locus.sense() == '.':
                locus = Locus(locus.chr(),locus.start(),locus.end()+extension,locus.sense(), locus.ID())
            if locus.sense() == '-':
                locus = Locus(locus.chr(),locus.start()-extension,locus.end(),locus.sense(),locus.ID())
            extendedReads.append(locus)
        if gffLocus.sense() == '+' or gffLocus.sense == '.':
            senseReads = filter(lambda x:x.sense() == '+' or x.sense() == '.',extendedReads)
            antiReads = filter(lambda x:x.sense() == '-',extendedReads)
        else:
            senseReads = filter(lambda x:x.sense() == '-' or x.sense() == '.',extendedReads)
            antiReads = filter(lambda x:x.sense() == '+',extendedReads)
    
        #at this point can output starts onto the GFF unless density is called

        if density:

            senseHash = defaultdict(int)
            antiHash = defaultdict(int)

            #filling in the readHashes             
            if sense == '+' or sense == 'both' or sense =='.':
                for read in senseReads:
                    for x in range(read.start(),read.end()+1,1):
                        senseHash[x]+=1
            if sense == '-' or sense == 'both' or sense == '.':
                #print('foo')
                for read in antiReads:
                    for x in range(read.start(),read.end()+1,1):
                        antiHash[x]+=1

            #now apply flooring and filtering for coordinates
            keys = uniquify(senseHash.keys()+antiHash.keys())
            if floor > 0:
                    
                keys = filter(lambda x: (senseHash[x]+antiHash[x]) > floor,keys)
            #coordinate filtering
            keys = filter(lambda x: gffLocus.start() < x < gffLocus.end(),keys)

            if clusterGram or matrix:
                clusterLine = [gffLocus.ID(),gffLocus.__str__()]
                if matrix:
                    binSize = (gffLocus.len()-1)/int(matrix)
                    nBins = int(matrix)
                if clusterGram:
                    nBins = gffLocus.len()/binSize
                if binSize == 0:
                    clusterLine+=['NA']*int(matrix)
                    newGFF.append(clusterLine)
                    continue
                n=0
                if gffLocus.sense() == '+' or gffLocus.sense() =='.' or gffLocus.sense() == 'both':
                    i = gffLocus.start()
                    
                    while n <nBins:
                        n+=1
                        binKeys = filter(lambda x: i < x < i+binSize,keys)
                        binDen = float(sum([senseHash[x]+antiHash[x] for x in binKeys]))/binSize
                        clusterLine+=[round(binDen/MMR,4)]
                        i = i+binSize
                else:
                    i = gffLocus.end()
                    while n < nBins:
                        n+=1
                        binKeys = filter(lambda x: i-binSize < x < i,keys)
                        binDen = float(sum([senseHash[x]+antiHash[x] for x in binKeys]))/binSize
                        clusterLine+=[round(binDen/MMR,4)]
                        i = i-binSize
                newGFF.append(clusterLine)
        
            #for regular old density calculation
            else:
                senseTotalDen = float(sum([senseHash[x] for x in keys]))/gffLocus.len()
                antiTotalDen = float(sum([antiHash[x] for x in keys]))/gffLocus.len()
                if rpm:
                    senseTotalDen = senseTotalDen/MMR
                    antiTotalDen = antiTotalDen/MMR
                if sense == 'both' or sense == '.':
                    if gffLocus.sense() == '+' or gffLocus.sense() == '.':
                        readLine = '+'+':%s' % (round(senseTotalDen,4)) + ';' +'-' + ':%s' % (round(antiTotalDen,4))
                    else:
                        readLine = '+'+':%s' % (round(antiTotalDen,4)) + ';' +'-' + ':%s' % (round(senseTotalDen,4))
                elif sense == '+':
                    readLine = '+'+':%s' % (round(senseTotalDen,4))
                elif sense == '-':
                    readLine = '-' + ':%s' % (round(antiTotalDen,4))
                newGFF.append(line + [readLine])             
        #if not cluster or density simply return reads 
        elif raw:
            if sense == 'both' or sense == '.':
                if gffLocus.sense() == '+' or gffLocus.sense() == '.':
                    readLine = '+'+':'+ join([str(locus.start()) for locus in senseReads],',') +';' + '-'+':'+ join([str(locus.start()) for locus in antiReads],',')
                else:
                    readLine = '+'+':'+ join([str(locus.start()) for locus in antiReads],',')+';'+'-'+':'+ join([str(locus.start()) for locus in senseReads],',')
            elif sense == '+':
                readLine = gffLocus.sense()+':'+ join([str(locus.start()) for locus in senseReads],',')
            elif sense == '-':
                readLine = string.translate(gffLocus.sense(),senseTrans)+':'+ join([str(locus.start()) for locus in antiReads],',')
            newGFF.append(line+[readLine])
        #if not raw and not density gives total
        else:


            if sense == 'both' or sense == '.':
                readLine = str((len(antiReads) + len(senseReads))/MMR)   
            elif sense == '+':
                readLine = str(len(senseReads)/MMR)
            elif sense == '-':
                readLine = str(len(antiReads)/MMR)
            newGFF.append(line+[readLine])
            
    return newGFF
        
                
                
            


            
                
    
def convertEnrichedRegionsToGFF(enrichedRegionFile):
    '''converts a young lab enriched regions file into a gff'''
    newGFF = []
    enrichedRegions = open(enrichedRegionFile,'r')
    header = enrichedRegions.readline()
    i = 0
    for line in enrichedRegions:
        line = line[:-1].split('\t')
        newLine = ['chr'+line[0],'row_'+str(i),line[4],line[1],line[2],'','.','','row_'+str(i),'']
        newGFF.append(newLine)
        i+=1
    return newGFF

        
#python bamToGFF.py --density --floor 0 -b test.sam.sorted.bam -g pol2_sample.gff -o pol2_sample_mapped.gff

def main():
    from optparse import OptionParser
    usage = "usage: %prog [options] -b [SORTED BAMFILE] -i [INPUTFILE] -o [OUTPUTFILE]"
    parser = OptionParser(usage = usage)
    #required flags
    parser.add_option("-b","--bam", dest="bam",nargs = 1, default=None,
                      help = "Enter .bam file to be processed.")
    parser.add_option("-i","--input", dest="input",nargs = 1, default=None,
                      help = "Enter .gff or ENRICHED REGION file to be processed.")
    #output flag
    parser.add_option("-o","--output", dest="output",nargs = 1, default=None,
                      help = "Enter the output filename.")
    #additional options
    parser.add_option("-s","--sense", dest="sense",nargs = 1, default='both',
                      help = "Map to '+','-' or 'both' strands. Default maps to both.")
    parser.add_option("-u","--unique", dest="unique",nargs = 1, default=0,
                      help = "Takes only unique sequence tags (avoids pcr repeats). must provide number of million uniquely mapping reads") #BJA changed default from None to 0
    parser.add_option("-d","--density", dest="density",action='store_true', default=False,
                      help = "Calculates a read density for each region, returns a single value per region")
    parser.add_option("-f","--floor", dest="floor",nargs =1, default=0,
                      help = "Sets a read floor threshold necessary to count towards density")    
    parser.add_option("-e","--extension", dest="extension",nargs = 1, default=200,
                      help = "Extends reads by n bp. Default value is 200bp")
    parser.add_option("-r","--rpm", dest="rpm",action = 'store_true', default=False,
                      help = "Normalizes density to reads per million (rpm)")
    parser.add_option("-t","--total", dest="total",action = 'store_true', default=False,
                      help = "Gives the total read count in the region. Overrides density, floor, and rpm flags")
    parser.add_option("-c","--cluster", dest="cluster",nargs = 1, default=None,
                      help = "Outputs a fixed bin size clustergram. user must specify bin size.")
    parser.add_option("-m","--matrix", dest="matrix",nargs = 1, default=None,
                      help = "Outputs a variable bin sized matrix. User must specify number of bins.")
    parser.add_option("-j","--jxn", dest="jxn",action = 'store_true', default=False,
                      help = "if flagged, includes jxn reads")
    (options,args) = parser.parse_args()

    print(options)
    print(args)

    if options.bam:
        bamFile = options.bam
        fullPath = os.path.abspath(bamFile)
        bamName = fullPath.split('/')[-1].split('.')[0]
        pathFolder = join(fullPath.split('/')[0:-1],'/')
        fileList = os.listdir(pathFolder)
        hasBai = False
        for fileName in fileList:
            if fileName.count(bamName) == 1 and fileName.count('.bai') == 1:
                hasBai = True

        if not hasBai:
            print('ERROR: no associated .bai file found with bam. Must use a sorted bam with accompanying index file')
            parser.print_help()
            exit()
   
    if options.sense:
        if ['+','-','.','both'].count(options.sense) == 0:
            print('ERROR: sense flag must be followed by +,-,.,both')
            parser.print_help()
            exit()

    if options.cluster and options.matrix:
        print('ERROR: Cannot specify both matrix and clustergram flags.')
        parser.print_help()
        exit()

    if options.matrix:
        try:
            int(options.matrix)
        except:
            print('ERROR: User must specify an integer bin number for matrix (try 50)')
            parser.print_help()
            exit()
            
    if options.cluster:
        try:
            int(options.cluster)
        except:
            print('ERROR: User must specify an integer bin size for clustergram (try 25)')
            parser.print_help()
            exit()

    
    
    if options.input and options.bam:
        inputFile = options.input
        if inputFile.split('.')[-1] != 'gff':
            print('converting file to a .gff')
            gffFile = convertEnrichedRegionsToGFF(inputFile)
        else:
            gffFile = inputFile

        bamFile = options.bam
        
        if options.output == None:
            output = os.getcwd() + inputFile.split('/')[-1]+'.mapped'
        else:
            output = options.output
        if options.cluster:
            print('mapping to GFF and making clustergram with fixed bin width')
            newGFF = mapBamToGFF(bamFile,gffFile,options.sense,options.unique,int(options.extension),options.floor,options.density,options.rpm,options.cluster,True,None,False,options.jxn)
        elif options.matrix:
            print('mapping to GFF and making a matrix with fixed bin number')
            newGFF = mapBamToGFF(bamFile,gffFile,options.sense,options.unique,int(options.extension),options.floor,options.density,options.rpm,25,None,options.matrix,False,options.jxn)
            
        else:
            print('mapping to GFF and returning reads')
            if options.total:

                newGFF = mapBamToGFF(bamFile,gffFile,options.sense,options.unique,int(options.extension),options.floor,options.density,options.rpm,25,None,None,False,options.jxn)
            else:
                newGFF = mapBamToGFF(bamFile,gffFile,options.sense,options.unique,int(options.extension),options.floor,options.density,options.rpm,25,None,None,True)
        unParseTable(newGFF,output,'\t')
    else:
        parser.print_help()
        

        
 
if __name__ == "__main__":
    main()
