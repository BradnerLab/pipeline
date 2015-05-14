#!/usr/bin/python
#pipeline_template.py

'''
The MIT License (MIT)

Copyright (c) 2015 Charles Lin

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

#generic pipeline template for human data


#==========================================================================
#=============================DEPENDENCIES=================================
#==========================================================================


import sys
sys.path.append('/ark/home/cl512/src/pipeline/')

import pipeline_dfci

#==========================================================================
#============================PARAMETERS====================================
#==========================================================================



projectName = 'projectName'
dataFile = '/grail/projects/%s/DATA_TABLE.txt' % (projectName)
genome ='hg19'
annotFile = '/ark/home/cl512/src/pipeline/annotation/%s_refseq.ucsc' % (genome)

#project folders
projectFolder = '/grail/projects/%s/' % (projectName) #PATH TO YOUR PROJECT FOLDER

#standard folder names
gffFolder ='%sgff/' % (projectFolder)
macsFolder = '%smacsFolder/' % (projectFolder)
macsEnrichedFolder = '%smacsEnriched/' % (projectFolder)
mappedEnrichedFolder = '%smappedEnriched/' % (projectFolder)
mappedFolder = '%smappedFolder/' % (projectFolder)
wiggleFolder = '%swiggles/' % (projectFolder)
metaFolder = '%smeta/' % (projectFolder)

#making folders
folderList = [gffFolder,macsFolder,macsEnrichedFolder,mappedEnrichedFolder,mappedFolder,wiggleFolder,metaFolder]

for folder in folderList:
    pipeline_dfci.formatFolder(folder,True)



#==========================================================================
#========================FORMATTING SAMPLE TABLE===========================
#==========================================================================

##THIS SECTION CREATES A DATA TABLE FROM A WHITEHEAD ANNOTATION SPREADSHEET

##give full path
##sampleTableFile = 'YOUR_WIGTC_ANNOTATION.xls' #<- the .xls file in the seq data folder provided by WI

#dirpath = ''  <- provide full path of folder containing raw seq files
##e.g. /ark/home/jr246/raw/130925_..../QualityScore/

##bamPath <- where we store our bams.  Must have write access if you want to call bowtie
##e.g. /ark/home/jr246/bam/
#bamPath = '/ark/home/jr246/bam/'

#pipeline_dfci.makePipelineTable(sampleTableFile,dirPath,bamPath,dataFile)

#dataDict = pipeline_dfci.loadDataTable(dataFile)

#namesList = dataDict.keys()

#print(namesList)

#==========================================================================
#=======================LOADING DATA ANNOTATION============================
#==========================================================================

##THIS SECTION LOADS A DATA TABLE.  MUST BE UNCOMMENTED FOR REST OF CODE TO WORK


#LOADING THE DATA TABLE
dataDict = pipeline_dfci.loadDataTable(dataFile)
print(dataDict.keys())

pipeline_dfci.summary(dataFile)

#==========================================================================
#==========================CALLING BOWTIE==================================
#==========================================================================

##THIS SECTION CALLS BOWTIE ON RAW READ FILES TO GENERATE SORTED AND INDEXED BAMS IN THE BAM FOLDER


#namesList = []  <- fill this in if you want to only map a subset of the data. otherwise leave blank

##SET LAUNCH TO False to debug
#pipeline_dfci.makeBowtieBashJobs(dataFile,namesList,launch=True)

#==========================================================================
#=============================CALL MACS====================================
#==========================================================================

##THIS SECTION CALLS THE MACS ERROR MODEL


#namesList = dataDict.keys()

#print(namesList)
#pipeline_dfci.callMacs(dataFile,macsFolder,namesList,overwrite=False,pvalue='1e-9')


#==========================================================================
#=======================FORMAT MACS OUTPUT=================================
#==========================================================================

##THIS SECTION FORMATS THE OUTPUT FROM MACS, CREATES THE MACSENRICHED FOLDER AND MOVES WIGGLES TO THE DESTINATION

#pipeline_dfci.formatMacsOutput(dataFile,macsFolder,macsEnrichedFolder,wiggleFolder,wigLink='/ark/wiggles/')


#==========================================================================
#====================ADDITIONAL PIPELINE ANALYSIS==========================
#==========================================================================




