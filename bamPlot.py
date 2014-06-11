#bamPlot.py

from utils import *
import string
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

#script that takes in a list of bams, makes an intermediate table file, and then calls R to make the plot

#borrows liberally from bamToGFF.py

#as of now the number of bins to sample the space is hard wired
nBins = 200

def loadAnnotFile(genome):

    '''
    load in the annotation and create a geneDict and transcription collection
    '''
    genomeDict = {
        'HG18':'/ark/home/cl512/pipeline/annotation/hg18_refseq.ucsc',
        'MM9': '/ark/home/cl512/pipeline/annotation/mm9_refseq.ucsc',
        'hg18':'/ark/home/cl512/pipeline./annotation/hg18_refseq.ucsc',
        'mm9': '/ark/home/cl512/pipeline/annotation/mm9_refseq.ucsc',
        'hg19':'/ark/home/cl512/pipeline/annotation/hg19_refseq.ucsc',
        'HG19':'/ark/home/cl512/pipeline/annotation/hg19_refseq.ucsc',
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

def mapGFFLineToAnnot(gffLine,outFolder,nBins,geneDict,txCollection,sense='both'):

    '''
    for every line produces a file with all of the rectangles to draw
    '''

    gffString = '%s_%s_%s_%s' % (gffLine[0],gffLine[6],gffLine[3],gffLine[4])
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

def mapBamToGFFLine(bamFile,name,gffLine,color,nBins,sense = 'both',unique = False,extension = 200,floor = 0,rpm = True,includeJxnReads = False):
    '''maps reads from a bam to a gff'''

    floor = int(floor)
    bam = Bam(bamFile)

    #=============
    #millionMappedReads
    if rpm:    
        MMR= round(float(bam.getTotalReads('mapped'))/1000000,4)
    else:
        MMR = 1
    
    print('using a MMR value of %s' % (MMR))
    
    senseTrans = string.maketrans('-+.','+-+')
    
    line = gffLine[0:9]
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

    #this is where we add Shit About the region like color, name etc...
    print(color)
    colorLine = color
    bamName = bamFile.split('/')[-1]
    clusterLine = [bamName,gffLocus.ID(),name,gffLocus.__str__()]+colorLine

    binSize = (gffLocus.len()-1)/int(nBins)

    if binSize == 0:
        clusterLine+=['NA']*int(nBins)
        return clusterLine

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

            
    return clusterLine


def callRPlot(nameTable,diagramTable,plotTable,yScale,plotStyle,fileName):
    '''
    calls the R plotting thingy
    '''

    cmd = 'R --no-save %s %s %s %s %s %s < /ark/home/cl512/pipeline/bamPlot.R' % (nameTable,diagramTable,plotTable,yScale,plotStyle,fileName)
    print('calling command %s' % (cmd))
    os.system(cmd)

    




def makeBamPlotTables(gff,genome,bamFileList,colorList,nBins,sense,unique,extension,floor,rpm,outFolder,yScale,names,plotStyle,fileName):

    '''
    makes a plot table for each line of the gff mapped against all the bams in the bamList
    '''
    
    #load in the gff
    if type(gff) == str:
        gff = parseTable(gff,'\t')

    #load in the annotation

    print('loading in annotation for %s' % (genome))

    geneDict,txCollection = loadAnnotFile(genome)

    #go line by line in the gff
    for gffLine in gff:
        gffString = '%s_%s_%s_%s' % (gffLine[0],gffLine[6],gffLine[3],gffLine[4])
        print('writing the gene diagram table for %s' % (gffLine[1]))
        mapGFFLineToAnnot(gffLine,outFolder,nBins,geneDict,txCollection,sense='both')

        outTable = []

        outTable.append(['BAM','GENE_ID','NAME','LOCUSLINE','COLOR1','COLOR2','COLOR3'] + ['bin_'+str(n) for n in range(1,int(nBins)+1,1)])
        
        for i in range(0,len(bamFileList),1):
            bamFile = bamFileList[i]
            name = names[i]
            
            color = colorList[i]
            print('getting data for location %s in dataset %s' % (gffLine[1],bamFile))
            newLine = mapBamToGFFLine(bamFile,name,gffLine,color,nBins,sense,unique,extension,floor,rpm,includeJxnReads = False)


            outTable.append(newLine)


        unParseTable(outTable,outFolder+gffString+'_plotTemp.txt','\t')
        diagramTable = outFolder+gffString+'_diagramTemp.txt'
        plotTable = outFolder+gffString+'_plotTemp.txt'
        nameTable = outFolder +gffString+'_nameTemp.txt'
        callRPlot(nameTable,diagramTable,plotTable,yScale,plotStyle,fileName)




def main():
    
    '''
    main run function
    '''

    from optparse import OptionParser

    usage = "usage: %prog [options] -b [SORTED BAMFILE(S)] -i [INPUTFILE] -o [OUTPUTFOLDER]"
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
    parser.add_option("-f","--floor", dest="floor",nargs =1, default=0,
                      help = "Sets a read floor threshold necessary to count towards density")    
    parser.add_option("-e","--extension", dest="extension",nargs = 1, default=200,
                      help = "Extends reads by n bp. Default value is 200bp")
    parser.add_option("-r","--rpm", dest="rpm",action = 'store_true', default=True,
                      help = "Normalizes density to reads per million (rpm) Default is True")
    parser.add_option("-u","--unique", dest="unique",action = 'store_false', default=True,
                      help = "Uses only unique reads")
    parser.add_option("-y","--yScale",dest="yScale",nargs =1, default = "relative",
                      help = "Choose either relative or uniform y axis scaling. options = 'relative,uniform' Default is relative scaling")
    parser.add_option("-n","--names",dest="names",nargs =1, default = None,
                      help = "Enter a comma separated list of names for your bams")
    parser.add_option("-p","--plot",dest="plot",nargs =1, default = "multiple",
                      help = "Choose either all lines on a single plot or multiple plots. options = 'single,multiple'")
    parser.add_option("-t","--title",dest ="title",nargs=1,default = '',
                      help = "Specify a title for the output plot, default will be the coordinate region")
                  


    (options,args) = parser.parse_args()

    print(options)
    print(args)
    
    if options.bam and options.input and options.genome:

        #bring in the bams
        bamFileList = options.bam.split(',')
        
        #bring in the gff
        try:
            gff = parseTable(options.input,'\t')
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
            gff = [gffLine]

        #bring in the genome
        genome = string.upper(options.genome)
        if ['HG19','HG18','MM9','RN5'].count(genome) == 0:
            print('ERROR: UNSUPPORTED GENOME TYPE %s. USE HG18, HG19, RN5, OR MM9' % (genome))
            parser.print_help()
            exit()

        #bring in the rest of the options
        
        #output
        outFolder = options.output
        if outFolder[-1] != '/':
            outFolder+='/'
        try:
            foo = os.listdir(outFolder)
        except OSError:
            print('ERROR: UNABLE TO FIND OUTPUT DIRECTORY %s' % (outFolder))
            exit()

                 
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
        
        #floor
        floor = int(options.floor)
        
        extension = int(options.extension)

        rpm = options.rpm

        unique = options.unique
        yScale = string.upper(options.yScale)
                                     
        fileName = options.title 
        #names
        if options.names:
            names = options.names.split(',')
        
            if len(names) != len(bamFileList):
                print('ERROR: NUMBER OF NAMES AND NUMBER OF BAMS DO NOT CORRESPOND')
                parser.print_help()
                exit()
        else:
            names = bamFileList

        #plot style
        plotStyle = string.upper(options.plot)
        if ['SINGLE','MULTIPLE'].count(plotStyle) == 0:
            print('ERROR: PLOT STYLE %s NOT AN OPTION' % (plotStyle))
            parser.print_help()
            exit()


        #now run!

        makeBamPlotTables(gff,genome,bamFileList,colorList,nBins,sense,unique,extension,floor,rpm,outFolder,yScale,names,plotStyle,fileName)

    else:
        parser.print_help()
        exit()

if __name__ == "__main__":
    main()

