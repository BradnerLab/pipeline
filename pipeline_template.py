#pipeline_template.py

#generic pipeline template for human data


#==========================================================================
#=============================DEPENDENCIES=================================
#==========================================================================


import sys
sys.path.append('/ark/home/cl512/src/pipeline/')

import pipeline_dfci
import datetime
import random
import string

#==========================================================================
#============================PARAMETERS====================================
#==========================================================================

dataFile = '/ark/home/USERNAME/ressrv19/projects/MYPROJECT/DATA_TABLE.txt' #PATH TO YOUR DATA TABLE
genome ='hg18'
annotFile = '/ark/home/cl512/src/pipeline/annotation/hg18_refseq.ucsc'

#project folders
projectFolder = '/ark/home/USERNAME/ressrv19/projects/MYPROJECT/' #PATH TO YOUR PROJECT FOLDER

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




