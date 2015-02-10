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

#pythonTemplate.py <- change to title of your script
#130801 <- date
#Name 


#Description:

#This is a generic python template that has functions from utils.py imported and can be used on CFCE1



#================================================================================
#=============================DEPENDENCIES=======================================
#================================================================================

import sys

print "Using python version %s" % sys.version


#importing utils package
sys.path.append('/ark/home/cl512/src/pipeline/')
import utils
import string
import os
import time
import numpy
from collections import defaultdict
#================================================================================
#============================GLOBAL PARAMETERS===================================
#================================================================================

#add locations of files and global parameters in this section

#set for the samtools string
samtoolsString = 'samtools'

#================================================================================
#===================================CLASSES======================================
#================================================================================

#user defined classes here

#================================================================================
#=================================FUNCTIONS======================================
#================================================================================

#write your specific functions here


def makeAnnotDict(annotFile):

    '''
    makes a dictionary keyed by guideID
    '''

    guideDict = defaultdict(str)
    geneDict = defaultdict(list)

    geckoAnnot = utils.parseTable(annotFile,'\t')
    
    for line in geckoAnnot[1:]:
        guideDict[line[1]] = line[0]
        geneDict[line[0]].append(line[1])

    return guideDict,geneDict



def makeFoldTable(annotFile,analysisName,testName,controlName,testMMR,controlMMR,testIdxFile,controlIdxFile,outputFolder,epsilon = 1):

    '''
    makes the fold table and writes to disk
    fold table is ranked by fold change
    first column is guideID, second column is gene name, third is fold change
    '''

    guideDict,geneDict = makeAnnotDict(annotFile)

    testIdx = utils.parseTable(testIdxFile,'\t')
    controlIdx = utils.parseTable(controlIdxFile,'\t')

    #for each guide, divide the count by the MMR then add 1 then take the log2 ratio

    outTable = [['GUIDE_ID','GENE','LOG2_RATIO',testName,controlName]]
    for i in range(len(testIdx)):

        guideID = testIdx[i][0]
        gene = guideDict[guideID]
        
        testCount = float(testIdx[i][2])/testMMR + epsilon
        controlCount = float(controlIdx[i][2])/controlMMR + epsilon

        log2Ratio = numpy.log2(testCount/controlCount)

        newLine = [guideID,gene,log2Ratio,round(testCount,4),round(controlCount,4)]

        outTable.append(newLine)

    outputFile = '%s%s_log2Ratio.txt' % (outputFolder,analysisName)
    utils.unParseTable(outTable,outputFile,'\t')
    return outputFile



#go from the fold table to a riger acceptable table

def makeRigerTable(foldTableFile,output=''):

    '''
    blah
    '''

    #need a table of this format
    rigerTable = [['Construct','GeneSymbol','NormalizedScore','Construct Rank','HairpinWeight']]
    #set weight to 1 for now

    foldTable = utils.parseTable(foldTableFile,'\t')

    constructOrder = utils.order([float(line[2]) for line in foldTable[1:]],decreasing=True)

    #make geneCountDict
    print("making gene count dictionary")
    geneCountDict= defaultdict(int)
    for line in foldTable[1:]:
        geneCountDict[line[1]] +=1

    print("iterating through constructs")
    constructRank = 1
    for i in constructOrder:
        rowIndex = i+1 # accounts for the header
        geneName = foldTable[rowIndex][1]
        if geneCountDict[geneName] == 1:
            print("Gene %s only has one guide RNA. Excluding from FRIGER analysis" % (geneName))
            continue

        newLine = foldTable[rowIndex][0:3] + [constructRank,1]
        rigerTable.append(newLine)
        constructRank += 1

    if len(output) == 0:
        output = string.replace(foldTableFile,'_log2Ratio.txt','_friger.txt')
    
    utils.unParseTable(rigerTable,output,'\t')

    return output


def callRiger(rigerTableFile,scoring='KSbyScore',output='',callRiger = True):

    '''
    calls riger using the KS scoring metric (default)
    '''
    rigerDirectory = '/raider/temp/riger/'

    rigerTableAbsFile = os.path.abspath(rigerTableFile)
    outputFolder = utils.getParentFolder(rigerTableAbsFile)
    if len(output) == 0:
        output = string.replace(rigerTableAbsFile,'_friger.txt','_friger_%s_out.txt' % (scoring))
    rigerBashFileName = string.replace(rigerTableAbsFile,'_friger.txt','_callRiger.sh')

    rigerBashFile = open(rigerBashFileName,'w')

    rigerBashFile.write('#!/usr/bin/bash\n')
    rigerBashFile.write('cd %s\n\n' % (rigerDirectory))
    rigerCmd = 'java -cp commons-cli-1.2.jar:rigerj-1.6.2.jar org.broadinstitute.rnai.rigerj.RigerJMain -scoringMethod %s -inputFile %s -outputFile %s' % (scoring,rigerTableAbsFile,output)
    rigerBashFile.write(rigerCmd)
    rigerBashFile.write('\n')

    rigerBashFile.close()
    print("WROTE RIGER CMD TO %s" % (rigerBashFileName))
    if callRiger == True:
        print("Calling RIGER with %s scoring method" % (scoring))
        print("RIGER CMD: %s" % (rigerCmd))
        os.system(rigerBashFileName)
    return rigerBashFileName

    
#================================================================================
#===============================MAIN RUN=========================================
#================================================================================


#write the actual script here


def main():


    '''
    main run function
    '''

    from optparse import OptionParser

    usage = "usage: %prog [options] -t [TEST_BAM] -c [CONTROL_BAM] -g [GENOME]"
    parser = OptionParser(usage = usage)
    #required flags
    parser.add_option("-t","--test", dest="test",nargs = 1, default=None,
                      help = "Enter the full path of the test bam")
    parser.add_option("-c","--control", dest="control",nargs = 1, default=None,
                      help = "Enter the full path of the control bam")
    parser.add_option("-g","--genome", dest="genome",nargs = 1, default=None,
                      help = "Enter the build for the GeCKO library (currently only supports geckov2)")


    #optional arguments
    parser.add_option("-n","--name",dest="name",nargs =1, default = 0,
                      help = "Comma separated test,control name")
    parser.add_option("-s","--scoring",dest="scoring",nargs =1, default = 'WtSum',
                      help = "Scoring method (KSbyScore,WtSum,SecondBestRank) defulat: WtSum")
    parser.add_option("-o","--output", dest="output",nargs = 1, default=None,
                      help = "Enter the full path of the output folder. Default is the current working directory")


    (options,args) = parser.parse_args()

    #three required parameters to get started
    if options.test and options.control and options.genome:

        #get the names of the datasets
        if options.name:
            if len(options.name.split(',')) == 2:
                [testName,controlName] = options.name.split(',')
            else:
                print("ERROR: Must provide a comma separated test,control name if using -n flag")
                parser.print_help()
                sys.exit()
        else:
            #try to extract names from file
            #strip extension from filename
            testName = options.test.split('/')[-1].split('.')[0]
            controlName = options.control.split('/')[-1].split('.')[0]

        #names
        print("using %s as name for test dataset" % (testName))
        print("using %s as name for control dataset" % (controlName))

        #get the analysis name
        analysisName = '%s_%s' % (testName,controlName)
        print("using %s as analysis name" % (analysisName))
        
        #get the scoring method
        scoringMethod = options.scoring
        if ['KSbyScore','WtSum','SecondBestRank'].count(scoringMethod)==0:
            print("ERROR: please specify one of the following scoring methods:('KSbyScore','WtSum','SecondBestRank') or leave blank (default WtSum)")
            parser.print_help()
            sys.exit()
                  
        
        #set up output folder
        if options.output:
            outputFolder = utils.formatFolder(options.output,True)
        else:
            outputFolder = utils.formatFolder('./%s/' % (analysisName),True)

        print("using %s as an output folder" % (outputFolder))

        #get the right annotation
        genomeDict = {'geckov2':'/grail/genomes/gecko/GeCKOv2/Annotation/Human_GeCKOv2_Library.txt',
                      }

        #load the annotation dictionary
        annotFile = genomeDict[string.lower(options.genome)]
        print("using %s as the annotation file" % (annotFile))
        
        #guideDict,geneDict = makeAnnotDict(annotFile)
        
        #now set up each bam
        testBam = utils.Bam(options.test)
        controlBam = utils.Bam(options.control)

        #get the MMR for each
        testMMR = round(float(testBam.getTotalReads())/1000000,4)
        controlMMR = round(float(controlBam.getTotalReads())/1000000,4)

        print("Test dataset: %s has an MMR of %s" % (testName,testMMR))
        print("Control dataset: %s has an MMR of %s" % (controlName,controlMMR))

        #now get the idxstats output
        testIdxFile = '%s%s_idxstats.txt' % (outputFolder,testName)
        testIdxCmd = '%s idxstats %s > %s' % (samtoolsString,options.test,testIdxFile)
        print("Test idxstats command:")
        print(testIdxCmd)
        os.system(testIdxCmd)

        controlIdxFile = '%s%s_idxstats.txt' % (outputFolder,controlName)
        controlIdxCmd = '%s idxstats %s > %s' % (samtoolsString,options.control,controlIdxFile)
        print("Control idxstats command:")
        print(controlIdxCmd)
        os.system(controlIdxCmd)

        print("Checking for output")
        if not utils.checkOutput(testIdxFile,0.1,5):
            print("ERROR: UNABLE TO GENERATE IDX OUTPUT FOR %s" % (options.test))
        print("Found test IdxStats file")
        if not utils.checkOutput(controlIdxFile,0.1,5):
            print("ERROR: UNABLE TO GENERATE IDX OUTPUT FOR %s" % (options.control))
        print("Found control IdxStats file")

        #now make the fold table

        foldTableFile =makeFoldTable(annotFile,analysisName,testName,controlName,testMMR,controlMMR,testIdxFile,controlIdxFile,outputFolder,epsilon = 1)
        
        print('writing output to %s' % (foldTableFile))
        
        print("MAING FRIGER TABLE")
        rigerTableFile = makeRigerTable(foldTableFile,output='')
        print('writing FRIGER table to %s' % (rigerTableFile))

        rigerBashFileName = callRiger(rigerTableFile,scoring=scoringMethod,output='',callRiger=True)

        

    else:
        parser.print_help()
        sys.exit()




if __name__ == "__main__":
    main()



