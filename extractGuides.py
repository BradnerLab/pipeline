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

'''
Extracts putative guide RNAs out of a fastq for processing in the GECKO pipeline
'''



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

from collections import defaultdict


#================================================================================
#============================GLOBAL PARAMETERS===================================
#================================================================================

#GUIDE EXTRACTION PARAMETERS
cutSeq = 'CACCG'
cutOffset = 5
guideLength = 20





#================================================================================
#===================================CLASSES======================================
#================================================================================

#user defined classes here

#================================================================================
#=================================FUNCTIONS======================================
#================================================================================


def extractGuideFastq(fastqFile,outputFolder='',gzip=False):

    '''
    takes a fastq and extracts candidate guide RNAs
    '''
    
    #get the full absolute path for the fastq File
    fastqFile = os.path.abspath(fastqFile)
    fastq = utils.open(fastqFile,'r')
    
    #get the fastq name and root
    if len(outputFolder) == 0:
        outputFolder = utils.getParentFolder(fastqFile)

    #makes sure the output folder exists
    utils.formatFolder(outputFolder,True)

    #grab the name info from the fastq
    fastqName = fastqFile.split('/')[-1]
    fastqRoot = string.replace(fastqName,'.fastq','')
    fastqRoot = string.replace(fastqRoot,'.gz','')
    
    #guideFastqFile output
    guideFastqFile = '%s%s.gecko.fastq' % (outputFolder,fastqRoot)
    guideFastq = utils.open(guideFastqFile,'w')

    print('processing %s' % (fastqName))
    print('million reads processed:')
    ticker = 0
    found = 0
    while True:
        
        if ticker%1000000 == 0:
            print(ticker/1000000)

        fastqLines = []

        #now load the fastq lines
        try:
            for i in range(4):
                fastqLines.append(fastq.next())
        except StopIteration:
            break

    
        #see if you can find a cut site
        seq = fastqLines[1].rstrip()
        try:
            cutPosition = seq.index(cutSeq)
            found+=1
        except ValueError:
            ticker+=1
            continue

        guideStart = cutPosition + cutOffset
        guideStop = guideStart + guideLength
        
        #pulling out the guide seq in the fastqLines
        fastqLines[1] = fastqLines[1][guideStart:guideStop] + '\n'
        fastqLines[3] = fastqLines[3][guideStart:guideStop] + '\n'

        for line in fastqLines:
            guideFastq.write(line)

        ticker+=1
        # if ticker == 100000:
        #     print(ticker)
        #     print(found)
        #     print(float(found)/float(ticker))
        #     break

    print('SUMMARY STATISTICS')
    print(ticker)
    print(found)
    print(float(found)/float(ticker))

    #close the fastq
    guideFastq.close()

    #gzip the fastq
    if gzip:
        os.system('gzip %s &' % (guideFastqFile))
        guideFastqFile += '.gz'

    return guideFastqFile








#================================================================================
#===============================MAIN RUN=========================================
#================================================================================


#write the actual script here


def main():


    '''
    main run function
    '''

    from optparse import OptionParser

    usage = "usage: %prog [options] -f [FASTQFILE]"
    parser = OptionParser(usage = usage)
    #required flags
    parser.add_option("-f","--fastq", dest="fastq",nargs = 1, default=None,
                      help = "Enter a comma separated list of fastq files to be extracted")


    #optional arguments
    parser.add_option("-o","--output",dest="output",nargs =1, default = 0,
                      help = "specify an output folder. default is the same directory as the input fastq")
    parser.add_option("-z","--gzip",dest="gzip",action = 'store_true', default = False,
                      help = "If flagged gzips the guide extracted fastq. Default is False")



    (options,args) = parser.parse_args()

    if not options.fastq:
        parser.print_help()
        exit()

    fastqFileList = options.fastq.split(',')

    if options.output:
        outputFolder = options.output
    else:
        outputFolder = ''

    if options.gzip:
        gzip = True
    else:
        gzip = False

    for fastqFile in fastqFileList:
        print("EXTRACTING GUIDES FROM %s" % (fastqFile))
        guideFastqFile = extractGuideFastq(fastqFile,outputFolder,gzip)
        print("GUIDES EXTRACTED TO %s" % (guideFastqFile))
        sys.stdout.write(guideFastqFile)
if __name__ == "__main__":
    main()



