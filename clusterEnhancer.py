#clusterEnhancer.py

#add X11 here


'''
The MIT License (MIT)

Copyright (c) 2014 Charles Lin

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


'''
program to perform 2D clustering by enhancer signal
can perform initial enhancer mapping or draw from a set of finished rose outputs
'''

#generic pipeline template for human data


#==========================================================================
#=============================DEPENDENCIES=================================
#==========================================================================


import sys
sys.path.append('/ark/home/cl512/src/pipeline/')

import pipeline_dfci
import utils
import string
import numpy
import os


#==========================================================================
#============================PARAMETERS====================================
#==========================================================================

###Change cl512 to jr246

#==========================================================================
#==============================FUNCTIONS===================================
#==========================================================================

def makeNameDict(dataFile,roseFolder,namesList=[]):

    '''
    for each name, check for the presence of an enriched file or  allEnhancer table
    these are the files required for enhancer clustering
    '''

    dataDict = pipeline_dfci.loadDataTable(dataFile)
    
    #draw the parent folder from the dataFile
    parentFolder = utils.getParentFolder(dataFile)

    #check to see if a rose folder exists already
    if utils.formatFolder(roseFolder,False):
        roseExists = True
        roseFolder = utils.formatFolder(roseFolder,False)
    else:
        roseExists = False
        roseFolder = utils.formatFolder(roseFolder,True)

    #check namesList to see if datasets exist
    if len(namesList) == 0:
        namesList = [name for name in dataDict.keys() if dataDict[name]['background'] != 'NONE']
        #this filters out control WCE datatsets

    #now check that all of the datasets at a minimum have a rose output OR enriched region file

    nameDict = {}
    for name in namesList:
        
        nameDict[name] = {}
        #assumes standard folder structure for enriched file
        enrichedFile = "%smacsEnriched/%s" % (parentFolder,dataDict[name]['enrichedMacs'])
        print enrichedFile
        try:
            foo = open(enrichedFile,'r')
            foo.close()
            nameDict[name]['enrichedFile'] = enrichedFile
        except IOError:
            nameDict[name]['enrichedFile'] = ''

        #roseOutput looks for standard format rose output
        #need an allEnhancers table and a region table to proceed
        #if the rose folder doesn't exist, don't bother
        if roseExists:

            roseOutputFiles = os.listdir("%s%s_ROSE" % (roseFolder,name))
            allEnhancerFileList = [x for x in roseOutputFiles if x.count("AllEnhancers.table.txt") == 1 and x[0] != '.' ] #no weird hidden or temp files
            if len(allEnhancerFileList) > 0:
                nameDict[name]['enhancerFile'] = "%s%s_ROSE/%s" % (roseFolder,name,allEnhancerFileList[0])
            else:
                nameDict[name]['enhancerFile'] = ''
        
        if nameDict[name]['enhancerFile'] == '' and nameDict[name]['enrichedFile'] =='':
            print "INSUFFICIENT DATA TO RUN ENAHNCER ANALYSIS ON %s. PLEASE MAKE SURE ROSE OUTPUT OR MACS ENRICHED REGION PEAKS FILE EXISTS" % (name)
            sys.exit()
    return nameDict



def launchEnhancerMapping(dataFile,nameDict,outputFolder):

    '''
    launches enhancer mapping if needed from enriched region files
    '''

    namesList = nameDict.keys()

    #check to see if everything is good, if so return True and call it a day
    if len([name for name in namesList if len(name) > 0]) == len(namesList):
        print "ENHANCER OUTPUT FOUND FOR ALL DATASETS"
        return nameDict

    #if not, have to call rose
    roseOutputFolder = utils.formatFolder(outputFolder + 'rose/',True)
    
    queueList =[]
    for name in namesList:

        #check to see if we need to call rose
        if nameDict[name]['enhancerFile'] == '':
     
            #get the enriched file
            enrichedFile = nameDict[name]['enrichedFile']
            #call rose
            print "CALLING ROSE FOR %s" % (name)
            bashFileName = pipeline_dfci.callRose(dataFile,'',roseOutputFolder,[name],[],enrichedFile)
            #os.system('bash %s &' % (bashFileName))
            #add name to queue list
            queueList,append(name)

    #now check for completion of datasets
            
    for name in queueList:

        #check for the AllEnhancers table
        enhancerFile = "%s%s_ROSE/%s_peaks_AllEnhancers.table.txt" % (roseOutputFolder,name,name)
        print "CHECKING FOR %s ROSE OUTPUT IN %s" % (name,enhancerFile)
        if utils.checkOutput(enhancerFile,5,60):
            
            print "FOUND ENHANCER OUTPUT FOR %s" % (name)
            nameDict[name]['enhancerFile'] = enhancerFile
        else:
            print "UNABLE TO FIND ENHANCER OUTPUT FOR %s. QUITTING NOW" % (name)
            sys.exit()

    return nameDict
        


def makeMedianDict(nameDict):

    '''
    for each dataset returns the median background subtracted enhancer signal
    '''

    medianDict = {}

    for name in nameDict:

        #open up the allenhancerTable
        enhancerTable = utils.parseTable(nameDict[name]['enhancerFile'],'\t')

        #assume header ends after line 5
        enhancerVector = [float(line[6]) - float(line[7]) for line in enhancerTable[6:]]

        medianDict[name] = numpy.median(enhancerVector)

    return medianDict
    

def makeSECollection(enhancerFile,name,superOnly = True):
    '''
    returns a locus collection from a super table
    top gives the number of rows
    '''
    enhancerTable = utils.parseTable(enhancerFile,'\t')
    enhancerLoci = []


    for line in enhancerTable:
        if line[0][0] == '#' or line[0][0] == 'R':
            continue
        else:

            if superOnly and int(line[-1]) == 0:
                break
            enhancerLoci.append(utils.Locus(line[1],line[2],line[3],'.',name+'_'+line[0]))

    return utils.LocusCollection(enhancerLoci,50)



def mergeCollections(nameDict,analysisName,output='',superOnly=True):

    '''
    merges them collections
    '''

    allLoci = []
    namesList = nameDict.keys()
    for name in namesList:

        allLoci += makeSECollection(nameDict[name]['enhancerFile'],name).getLoci()

    mergedCollection = utils.LocusCollection(allLoci,50)

    #stitch the collection together
    stitchedCollection = mergedCollection.stitchCollection()

    stitchedLoci = stitchedCollection.getLoci()

    #sort by size and provide a unique ID

    sizeList = [locus.len() for locus in stitchedLoci]

    sizeOrder = utils.order(sizeList,decreasing=True)
    
    orderedLoci = [stitchedLoci[i] for i in sizeOrder]

    for i in range(len(orderedLoci)):
        orderedLoci[i]._ID = 'merged_%s_%s' % (analysisName,str(i+1))

    mergedGFF = []
    for locus in orderedLoci:
        newLine = [locus.chr(),locus.ID(),'',locus.start(),locus.end(),'',locus.sense(),'',locus.ID()]
        mergedGFF.append(newLine)


    if len(output) == 0:
        return mergedGFF
    else:
        print "writing merged gff to %s" % (output)
        utils.unParseTable(mergedGFF,output,'\t')
        return output



def mapMergedGFF(dataFile,nameDict,mergedGFFFile,analysisName,outputFolder):

    '''
    calls rose on the mergedGFFFile for all datasets
    '''
    dataDict= pipeline_dfci.loadDataTable(dataFile)
    roseParentFolder = "%srose/" % (outputFolder)
    gffName = mergedGFFFile.split('/')[-1].split('.')[0]
    bashFileName = "%srose/%s_roseCall.sh" % (outputFolder,analysisName)
    #namesList is just the first dataset
    #extrmap will have to have all other datasets + their backgrounds




    namesList = nameDict.keys()
    extraMap = []
    for name in namesList[1:]:
        
        backgroundName = dataDict[name]['background']
        extraMap+=[name,backgroundName]


    #first check to see if this has already been done
    mergedRegionMap = "%srose/%s_ROSE/%s_0KB_STITCHED_ENHANCER_REGION_MAP.txt" % (outputFolder,namesList[0],gffName)
    if utils.checkOutput(mergedRegionMap,1,1):
        return mergedRegionMap



    bashFileName = pipeline_dfci.callRose(dataFile,'',roseParentFolder,[namesList[0]],extraMap,mergedGFFFile,0,0,bashFileName) 
    
    bashCommand = "bash %s" % (bashFileName)
    os.system(bashCommand)
    print "Running enhancer mapping command:\n%s" % (bashCommand)


    if utils.checkOutput(mergedRegionMap,5,60):
        return mergedRegionMap
    else:
        print "UNABLE TO CALL ROSE ENHANCER MAPPING ON CONSENSUS ENHANCER FILE %s.\nEXITING NOW" % (mergedGFFFile)
        sys.exit()
    

def makeEnhancerSignalTable(mergedRegionMap,medianDict,analysisName,genome,outputFolder):

    '''
    makes a table where each row is an enhancer and each column is the log2 
    background corrected signal vs. median
    '''

    #load in the region map
    regionMap = utils.parseTable(mergedRegionMap,'\t')
    namesList = medianDict.keys()
    signalTable = [['REGION_ID','CHROM','START','STOP','NUM_LOCI','CONSTITUENT_SIZE'] + namesList]
    for line in regionMap[1:]:

        newLine = line[0:6]
        for i in range(len(namesList)):
            enhancerIndex = (i*2) + 6
            controlIndex = (i*2) + 7
            enhancerSignal = float(line[enhancerIndex]) - float(line[controlIndex])
            if enhancerSignal < 0:
                enhancerSignal = 0
            enhancerSignal = enhancerSignal/medianDict[namesList[i]]
            newLine.append(enhancerSignal)

        signalTable.append(newLine)

    outputFile = "%s%s_%s_signalTable.txt" % (outputFolder,genome,analysisName)
    print "WRITING MEDIAN NORMALIZED SIGNAL TABLE TO %s" % (outputFile)
    utils.unParseTable(signalTable,outputFile,'\t')
    return outputFile

    

def callRScript(genome,outputFolder,analysisName,signalTableFile):

    '''
    calls the R script to do clustering and heatmap
    '''
            
            
    clusterTable = "%s%s_%s_clusterTable.txt" % (outputFolder,genome,analysisName)

    rCmd = 'R --no-save %s %s %s %s < /ark/home/cl512/pipeline/clusterEnhancer.R' % (genome,outputFolder,analysisName,signalTableFile)
    print("Calling command %s" % rCmd)

    os.system(rCmd)

    print "Checking for cluster table output at %s" % (clusterTable)
    if utils.checkOutput(clusterTable,1,5):

        return clusterTable
    else:
        print "ERROR: CLUSTERING TABLE FAILED TO GENERATE"
        sys.exit()

#==========================================================================
#====================4. MAP H3K27AC TO CONSENSUS SUPERS====================
#==========================================================================


# ### We can use the callRose function to map all of the datasets to the merged supers

# parentFolder = '%smergedRose/' % (projectFolder)
# namesList = ['LY1_H3K27AC'] ### we will just rank wby LY1 for simplicity
# extraMap = dataDict.keys() ### now we map all data to these regions
# inputFile = gffFolder + 'HG18_DLBCL_SUPERS_-0_+0.gff'  ### this is the merged supers gff

# pipeline_dfci.callRose(dataFile,macsEnrichedFolder,parentFolder,namesList,extraMap,inputFile)



#==========================================================================
#==============5. CREATE THE H3K27AC VECTOR FOR COMBINED SUPERS============
#==========================================================================

### Use the *_ENHANCER_REGION_MAP.txt file in the cominbed supers rose output folder
### This has signal for h3k27ac and wce for each dataset at each super-enhancer
###simplify this table by subtracting the background from each k27ac and make a new table
### each row here is an enhancer ID and the columns should have background subtracted k27ac values
### then bring back your median signal table and divide each value by the correct median
### to generate median fold values
### use numpy.log2 function to get the log2 median fold values

###x = 0.2
###y = numpy.log2(x) <- returns the log transformed value




#==========================================================================
#==============6. IDENTIFY SUPER ASSOCIATED GENES==========================
#==========================================================================

###run the following commands on the terminal

### cd /ark/home/cl512/rose/ <- change into the rose source code folder

### python ROSE_geneMapper.py -g HG18 -i  XXX -w 500000 <- XXX is the full path to the AllEnhancers.table.txt in the combined supers rose output.  this will create a file that's AllEnhancers_ENHANCER_TO_GENE.txt  you can use this to identify genes associated with each super
### the -w tells the algorithm to search +/- 500kb from the center of the super for nearby genes





#==========================================================================
#=================7. CREATE THE EXPRESSION TABLE===========================
#==========================================================================

### the expression folder contains a table with all expression values JQ1.hsentrezg.rma.gct
### there are a lot of columns here that are uneccessary
### you want a single expression value for each gene in each sample

### sampleKey.txt tells you which column IDs correspond to each dataset
### use the relevant columns to get the mean expression for each gene
### numpy.mean(l) gives the mean of any list l



#==========================================================================
#=================8. CORRELATE H3K27AC TO EXPRESSION=======================
#==========================================================================

###for each enhancer, we want to correlate its h3k27ac vector to the expression vector of its associated genes


### you should by now have all the files you need
### use dictionaries to help you assign various expression or h3k27ac vectors to genes or enhancerIDs

### use dictionaries to assign genes to each enhancerID

### use the following function to correlate the k27ac vector to the gene expression vector

###rom scipy.stats.stats import pearsonr

###p = pearsonr(x,y)  where x and y are lists of equal length
### p is a tuple where p[0] is the correlation statistic and p[1] is the significance of correlation

### for each enhancerID and gene pair, write a table like this

'''
[['ENHANCER_ID','ASSOCIATED_GENE','PEARSON_COR','PEARSON_SIG'],
['ENH_1','GENE_A',.79,0.05]]
'''

###this is your final output table




#==========================================================================
#=============================MAIN METHOD==================================
#==========================================================================


def main():

    from optparse import OptionParser

    usage = "usage: %prog [options] -d [DATA_FILE] -n [NAMES_LIST] -r [ROSE_FOLDER] -o [OUTPUTFOLDER]"
    parser = OptionParser(usage = usage)
    #required flags
    parser.add_option("-d","--data", dest="data",nargs = 1, default=None,
                      help = "Enter a data file for datasets to be processed")
    parser.add_option("-o","--output",dest="output",nargs =1, default = None,
                      help = "specify an output folder to write results to")

    #additional options
    parser.add_option("-i","--input", dest="input",nargs = 1, default=None,
                      help = "Enter a comma separated list of names to analyze. Default will be all datasets")

    parser.add_option("-n","--name", dest="name",nargs=1,default=None,
                      help = "Enter a name for the analysis")

    parser.add_option("-r","--rose", dest="rose",nargs = 1, default=None,
                      help = "Enter a folder to detect or write rose output")

    parser.add_option("-a","--all", dest="all",action = 'store_true', default=False,
                      help = "flag to run analysis on ALL enhancers (this is much slower)")

    (options,args) = parser.parse_args()

    print(options)
    print(args)
    
    if options.data and options.output:

        #check to see if minimum arguments are met

        #pull in arguments
        
        #pull in the datafile and create a datadict
        dataFile = options.data




        #now the output folder
        outputFolder = utils.formatFolder(options.output,True) #check and create the output folder
        #now the rose folder
        if options.rose:
            roseFolder = options.rose
        else:
            roseFolder = "%srose/" % (outputFolder)

        if options.input:
            namesList = options.input.split(',')
        else:
            namesList = []

        #get the genome
        dataDict = pipeline_dfci.loadDataTable(dataFile)
        genome = dataDict[dataDict.keys()[0]]['genome']

        #check if using only supers
        if options.all:
            superOnly = False
        else:
            superOnly = True

        #get the anlysis name
        if options.name:
            analysisName = options.name
        else:
            analysisName = "enhancers"

        #=====================================================
        #=================SUMMARIZE INPUTS====================
        #=====================================================
        
        print "WORKING IN GENOME %s" % (genome)
        print "DRAWING DATA FROM %s AND ROSE FOLDER %s" % (dataFile,roseFolder)
        print "USING %s AS THE OUTPUT FOLDER" % (outputFolder)
        print "STARTING ANALYSIS ON THE FOLLOWING DATASETS:"
        print namesList

        #=====================================================
        #==============ESTABLISH ALL WORKING FILES============
        #=====================================================

        print "\n\n\nESTABLISHING WORKING FILES"
        nameDict = makeNameDict(dataFile,roseFolder,namesList)

            
        print nameDict
        
        #=====================================================
        #==============LAUNCH ENHANCER MAPPING================
        #=====================================================
        
        print "\n\n\nLAUNCHING ENHANCER MAPPING (IF NECESSARY)"
        nameDict = launchEnhancerMapping(dataFile,nameDict,outputFolder)
        print nameDict


        #=====================================================
        #====================GET MEDIAN SIGNAL================
        #=====================================================
        
        print "\n\n\nGETTING MEDIAN ENHANCER SIGNAL FROM EACH SAMPLE"
        medianDict = makeMedianDict(nameDict)

        print medianDict
        
        #=====================================================
        #====================MERGING ENHANCERS================
        #=====================================================
        
        print "\n\n\nIDENTIFYING CONSENSUS ENHANCER REGIONS"
        mergedGFFFile = "%s%s_%s_-0_+0.gff" % (outputFolder,genome,analysisName)
        mergeCollections(nameDict,analysisName,mergedGFFFile,superOnly)


        #=====================================================
        #===============MAP TO MERGED REGIONS=================
        #=====================================================

        print "\n\n\nMAPPING DATA TO CONSENSUS ENHANCER REGIONS"
        mergedRegionMap = mapMergedGFF(dataFile,nameDict,mergedGFFFile,analysisName,outputFolder)
        
        #=====================================================
        #==============CORRECT FOR MEDIAN SIGNAL==============
        #=====================================================

        print "\n\n\nCREATING ENHANCER SIGNAL TABLE"
        signalTableFile = makeEnhancerSignalTable(mergedRegionMap,medianDict,analysisName,genome,outputFolder)
        #=====================================================
        #===============CALL CLUSTERING R SCRIPT==============
        #=====================================================

        print "\n\n\nGENERATING CLUSTERING OUTPUT"
        clusterTableFile = callRScript(genome,outputFolder,analysisName,signalTableFile)
        #output should be
        #png of cluster gram with rows as genes
        #png of cluster gram of samples w/ tree
        #ordered table w/ cluster assignment
        #similarity matrix for samples

        #=====================================================
        #=============GENE MAPPING BY CLUSTER=================
        #=====================================================

        os.chdir('/ark/home/cl512/rose/')
        cmd = 'python /ark/home/cl512/rose/ROSE_geneMapper.py -g %s -i %s' % (genome,clusterTableFile)
        os.system(cmd)

        print "FINISHED"


    else:
        parser.print_help()
        sys.exit()




if __name__ == "__main__":
    main()
