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
gffFolder ='%sgff/'
macsFolder = '%smacsFolder/' % (projectFolder)
macsEnrichedFolder = '%smacsEnriched/' % (projectFolder)
mappedEnrichedFolder = '%smappedEnriched/' % (projectFolder)
mappedFolder = '%smappedFolder/' % (projectFolder)
wiggleFolder = '%swiggles/' % (projectFolder)
metaFolder = '%smeta/' % (projectFolder)



#LOADING THE DATA TABLE
dataDict = pipeline_dfci.loadDataTable(dataFile)
print(dataDict.keys())

#==========================================================================
#========================FORMATTING ANNOTATION=============================
#==========================================================================

#sampleTableFile = 'YOUR_WIGTC_ANNOTATION.xls'

#pipeline_dfci.makePipelineTable(sampleTableFile,dirPath,bamPath,dataFile)

#dataDict = pipeline_dfci.loadDataTable(dataFile)

#namesList = dataDict.keys()

print(namesList)


#==========================================================================
#==========================CALLING BOWTIE==================================
#==========================================================================

#namesList = []

#pipeline_dfci.makeBowtieBashJobs(dataFile,namesList,launch=True)

#==========================================================================
#=============================CALL MACS====================================
#==========================================================================

#namesList = dataDict.keys()

#print(namesList)


#pipeline_dfci.callMacs(dataFile,macsFolder,namesList,overwrite=False,pvalue='1e-9')


#pipeline_dfci.formatMacsOutput(dataFile,macsFolder,macsEnrichedFolder,wiggleFolder,wigLink='/ark/wiggles/')



