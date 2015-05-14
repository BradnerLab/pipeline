#!/usr/bin/python
#given a fastq location, create a .sh script to run mapping through generation of sorted bam

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


#===================================================================
#========================MODULES AND DEPENDENCIES===================
#===================================================================

import string

import random
import utils

#===================================================================
#==========================GLOBAL PARAMETERS========================
#===================================================================

#command arguments
bowtieString = 'bowtie2'
samtoolsString = 'samtools'
tempParentFolder = '/grail/BOWTIE_TEMP/'
fastqcString = '/usr/local/FastQC/fastqc'
fastqDelimiter = '::'

#tempParentFolder = '/mnt/d0-0/share/bradnerlab/projects/anna/BOWTIE_TEMP/'

#===================================================================
#=============================FUNCTIONS=============================
#===================================================================


def stripExtension(fileName):

    '''
    tries to strip the extension of a filename
    can strip .tar.gz, .txt, .fastq,.gz,.zip
    '''
    extensionList = ['.tar.gz','.tar', '.txt','.fastq','.fasta','.gz','.zip']

    for extension in extensionList:
        fileName = fileName.replace(extension,'')

    return fileName
    

#print stripExtension('CGATGT-s_8_1_sequence.txt')


def makeFileNameDict(fastqFile,genome,tempString,tempParentFolder,finalFolder,linkFolder,uniqueID='',pairedEnd = False):

    '''
    creates a dictionary w/ all filenames
    '''


    fileNameDict = {}
    if pairedEnd:
        fastqFileList = fastqFile.split(fastqDelimiter)
        fileNameDict['fastqFile_1'] = fastqFileList[0]
        fileNameDict['fastqFile_2'] = fastqFileList[1]
    else:
        fileNameDict['fastqFile'] = fastqFile

    if uniqueID == '':
        fastqName = fastqFile.split('/')[-1]
        fastqName = stripExtension(fastqName)

    else:
        fastqName = uniqueID

    fileNameDict['fastqName'] = fastqName

    #make the temp Folder
    tempFolder = tempParentFolder + 'bwt2_' + fastqName + tempString + '/'    
    fileNameDict['tempFolder'] = tempFolder

    if pairedEnd:
        tempFastqFile1 = tempFolder + fastqName + '_1.rawFastq'
        fileNameDict['tempFastqFile_1'] = tempFastqFile1

        tempFastqFile2 = tempFolder + fastqName + '_2.rawFastq'
        fileNameDict['tempFastqFile_2'] = tempFastqFile2

    else:
        tempFastqFile = tempFolder + fastqName + '.rawFastq'
        fileNameDict['tempFastqFile'] = tempFastqFile

    tempSamFile = tempFolder + fastqName + '.sam'
    fileNameDict['tempSamFile'] = tempSamFile

    tempBamFile = tempFolder + fastqName + '.bam'
    fileNameDict['tempBamFile'] = tempBamFile

    tempSortedBamFile = tempFolder + fastqName + '.%s.bwt2.sorted' % (genome)
    fileNameDict['tempSortedBamFile'] = tempSortedBamFile

    sortedSamFile = fastqName + '.%s.bwt2.sorted.sam' % (genome)
    fileNameDict['sortedSamFile'] = sortedSamFile

    groupHeader = tempFolder + fastqName + '.%s.bwt2' % (genome)
    fileNameDict['groupHeader'] = groupHeader

    fileNameDict['finalFolder'] = finalFolder

    fileNameDict['linkFolder'] = linkFolder

    return fileNameDict




def extractFastqCmd(fileNameDict,pairedEnd = False):

    '''
    creates a command to extract/copy the fastq to a temp location
    '''
    if pairedEnd:
        fastqList = []
        fastqFile1 = fileNameDict['fastqFile_1']
        tempFastqFile1 = fileNameDict['tempFastqFile_1']
        fastqList.append([fastqFile1,tempFastqFile1])

        fastqFile2 = fileNameDict['fastqFile_2']
        tempFastqFile2 = fileNameDict['tempFastqFile_2']
        fastqList.append([fastqFile2,tempFastqFile2])
        
    else:
        fastqFile = fileNameDict['fastqFile']
        tempFastqFile = fileNameDict['tempFastqFile']
        fastqList = [[fastqFile,tempFastqFile]]

    cmdList = []
    for [fastqFile,tempFastqFile] in fastqList:
    #there are 3 possibilities, a gzipped, tarballed, or naked fastq
        if string.lower(fastqFile).count('tar.gz') == 1:
            cmd = "tar --strip-components 5 --to-stdout -xzvf %s > %s" % (fastqFile,tempFastqFile)
        elif string.lower(fastqFile).count('tar') == 1:
            cmd = "tar -xzvf %s > %s" % (fastqFile,tempFastqFile)
        elif string.lower(fastqFile.split('.')[-1]) == 'gz':
            cmd = 'cp %s %s.gz\n' % (fastqFile,tempFastqFile)
            cmd+= 'gunzip %s.gz' % (tempFastqFile)
        else:
            cmd = 'cp %s %s' % (fastqFile,tempFastqFile)

        cmdList.append(cmd)
    fullCmd = string.join(cmdList,'\n')

    return fullCmd

def runFastQC(fastqcString,fileNameDict,pairedEnd = False):

    '''
    cmd to run fastqc
    '''
    if pairedEnd:
        fastqName = fileNameDict['fastqName']
        tempFastqFile1 = fileNameDict['tempFastqFile_1']
        tempFastqFile2 = fileNameDict['tempFastqFile_2']
        finalFolder = fileNameDict['finalFolder']

        if finalFolder[-1] != '/':
            finalFolder+='/'
        finalFolder1 = finalFolder + '%s_1_fastqc' % (fastqName)
        finalFolder2 = finalFolder + '%s_2_fastqc' % (fastqName)
        cmd = 'mkdir %s\n' % (finalFolder1)
        cmd += 'mkdir %s\n' % (finalFolder2)
        cmd += '%s -o %s %s\n' % (fastqcString,finalFolder1,tempFastqFile1)
        cmd += '%s -o %s %s' % (fastqcString,finalFolder2,tempFastqFile2)

    else:
        fastqName = fileNameDict['fastqName']
        tempFastqFile = fileNameDict['tempFastqFile']
        finalFolder = fileNameDict['finalFolder']

        if finalFolder[-1] != '/':
            finalFolder+='/'
        finalFolder += '%s_fastqc' % (fastqName)
        cmd = 'mkdir %s\n' % (finalFolder)
        cmd += '%s -o %s %s' % (fastqcString,finalFolder,tempFastqFile)

    return cmd

def bowtieCmd(bowtieString,paramString,bowtieIndex,fileNameDict,pairedEnd=False):

    '''
    creates the bowtie command call
    '''

    
    #calling bowtie
    if pairedEnd:

        tempFastqFile1 = fileNameDict['tempFastqFile_1']
        tempFastqFile2 = fileNameDict['tempFastqFile_2']
        tempSamFile = fileNameDict['tempSamFile']
        cmd = "%s %s -x %s -1 %s -2 %s -S %s" % (bowtieString,paramString,bowtieIndex,tempFastqFile1,tempFastqFile2,tempSamFile)

    else:
        tempFastqFile = fileNameDict['tempFastqFile']
        tempSamFile = fileNameDict['tempSamFile']

        cmd = "%s %s -x %s -U %s -S %s" % (bowtieString,paramString,bowtieIndex,tempFastqFile,tempSamFile)
    return cmd


def removeTempFastqCmd(fileNameDict,pairedEnd = False):

    '''
    removes the temp fastq
    '''
    if pairedEnd:

        tempFastqFile1 = fileNameDict['tempFastqFile_1']    
        tempFastqFile2 = fileNameDict['tempFastqFile_2']    
        cmd = '/bin/rm -f %s\n' % (tempFastqFile1)
        cmd += '/bin/rm -f %s' % (tempFastqFile2)

    else:
        tempFastqFile = fileNameDict['tempFastqFile']    
        cmd = '/bin/rm -f %s' % (tempFastqFile)
    return cmd

#generate a bam file

def generateTempBamCmd(samtoolsString,fileNameDict):

    '''
    uses samtools to convert the sam to a bam
    '''
    tempSamFile = fileNameDict['tempSamFile']
    tempBamFile = fileNameDict['tempBamFile']

    cmd = "%s view -bS '%s' > '%s'" % (samtoolsString,tempSamFile,tempBamFile)
    
    return cmd

#change into temp directory

def changeTempDir(fileNameDict):

    '''
    changes into the temp directory
    '''
    tempFolder = fileNameDict['tempFolder']
    cmd = "cd %s" % (tempFolder)
    return cmd

#sort
def sortBamCmd(samtoolsString,fileNameDict):

    '''
    uses smatools to sort the bam
    '''
    tempBamFile = fileNameDict['tempBamFile']
    tempSortedBamFile = fileNameDict['tempSortedBamFile']

    cmd = "%s sort '%s' '%s'" % (samtoolsString,tempBamFile,tempSortedBamFile)
    return cmd


#index
def indexBamCmd(samtoolsString,fileNameDict):

    '''
    uses samtools to index the bam
    '''

    tempSortedBamFile=fileNameDict['tempSortedBamFile']
    cmd = "%s index '%s.bam'" % (samtoolsString,tempSortedBamFile)

    return cmd

def rmSamCmd(fileNameDict):

    '''
    remove the sam
    '''
    tempSamFile = fileNameDict['tempSamFile']
    cmd = "/bin/rm -f '%s'" % (tempSamFile)
    return cmd

def mvSamCmd(fileNameDict):

    '''
    rename and move the sam
    '''
    tempSamFile = fileNameDict['tempSamFile']
    finalFolder = fileNameDict['finalFolder']
    sortedSamFile=fileNameDict['sortedSamFile']
    cmd = "mv %s %s%s" % (tempSamFile,finalFolder,sortedSamFile)
    return cmd

def mvBamCmd(fileNameDict):

    '''
    moves and renames the bam w/o the temp string
    '''
    groupHeader = fileNameDict['groupHeader']
    finalFolder = fileNameDict['finalFolder']

    cmd = "mv %s* %s" % (groupHeader,finalFolder)
    
    return cmd

def linkBamCmd(fileNameDict):

    '''
    moves and renames the bam w/o the temp string
    '''
    groupHeader = fileNameDict['groupHeader']
    finalFolder = fileNameDict['finalFolder']
    linkFolder = fileNameDict['linkFolder']
    groupName = groupHeader.split('/')[-1]

    cmd = "ln %s%s* %s" % (finalFolder,groupName,linkFolder)
    
    return cmd


def rmTempFiles(fileNameDict):

    '''
    removes everything left in the temp folder
    '''
    groupHeader = fileNameDict['groupHeader']
    cmd = "/bin/rm -f '%s*'" % (groupHeader)
    return cmd

#===================================================================
#=============================MAIN==================================
#===================================================================

def main():

    '''
    main run function
    '''

    from optparse import OptionParser

    usage = "usage: %prog [options] -f [FASTQFILE] -g [GENOME] -u [UNIQUEID] -o [OUTPUTFOLDER]"
    parser = OptionParser(usage = usage)
    #required flags
    parser.add_option("-f","--fastq", dest="fastq",nargs = 1, default=None,
                      help = "Enter the full path of a fastq file to be mapped")
    parser.add_option("-g","--genome",dest="genome",nargs =1, default = None,
                      help = "specify a genome, options are hg19,hg18, mm9 or geckov2 right now")
    parser.add_option("-u","--unique",dest="unique",nargs =1, default = None,
                      help = "specify a uniqueID")
    parser.add_option("-o","--output",dest="output",nargs =1, default = None,
                      help = "Specify an output folder")


    #optional arguments
    parser.add_option("--param",dest="paramString",nargs =1, default = '',
                      help = "A string of bowtie parameters")
    parser.add_option("--link-folder",dest="linkFolder",nargs =1, default = None,
                      help = "Specify a folder to symlink the bam")
    parser.add_option("-p","--paired",dest="paired",action='store_true',default = False,
                      help = "Flag for paired end data")
    parser.add_option("-S","--sam",dest="sam",action='store_true',default = False,
                      help = "Flag to save sam")
    parser.add_option("-q","--qc",dest="qc",action='store_true',default = False,
                      help = "Flag to run fastqc")



    (options,args) = parser.parse_args()

    if not options.fastq or not options.genome or not options.unique or not options.output:
        parser.print_help()
        exit()


    #retrive the arguments
    fastqFile = options.fastq
    genome = string.lower(options.genome)
    uniqueID = options.unique
    outputFolder = options.output
    
    #make the output folder
    outputFolder = utils.formatFolder(outputFolder,True)

    #retrieve optional arguments
    paramString = options.paramString
    if options.linkFolder:

        linkFolder = options.linkFolder
    else:
        linkFolder =''
    pairedEnd = options.paired

    #get the bowtie index
    bowtieDict = {
        'mm9':'/raider/index/mm9/Bowtie2Index/genome',
        'hg19':'/raider/index/hg19/Bowtie2Index/genome',
        'hg18':'/grail/genomes/Homo_sapiens/human_gp_mar_06_no_random/bowtie/hg18',
        'geckov2':'/grail/genomes/gecko/GeCKOv2/Sequence/Bowtie2Index/gecko',
        'ribo':'/raider/temp/rDNA/hg19_45S_index/genome',
        'hg19_ribo':'/grail/genomes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index_ribo/genome',
        }

    bowtieIndex = bowtieDict[string.lower(genome)]

    #get the temp string
    tempString = '_%s' % str(random.randint(1,10000))
    
    fileNameDict = makeFileNameDict(fastqFile,genome,tempString,tempParentFolder,outputFolder,linkFolder,uniqueID,pairedEnd)

    #open the bashfile to write to
    bashFileName = "%s%s_bwt2.sh" % (outputFolder,uniqueID)
    bashFile = open(bashFileName,'w')

    #shebang
    bashFile.write('#!/usr/bin/bash\n')

    #make temp directory
    cmd = 'mkdir %s' % (fileNameDict['tempFolder'])
    bashFile.write(cmd+'\n')

    #extract fastq
    cmd = extractFastqCmd(fileNameDict,pairedEnd)
    bashFile.write(cmd+'\n')

    #call fastqc
    if options.qc:
        cmd =runFastQC(fastqcString,fileNameDict,pairedEnd)
        bashFile.write(cmd+'\n')

    #call bowtie
    cmd = bowtieCmd(bowtieString,paramString,bowtieIndex,fileNameDict,pairedEnd)
    bashFile.write(cmd+'\n')

    #remove temp fastq
    cmd = removeTempFastqCmd(fileNameDict,pairedEnd)
    bashFile.write(cmd+'\n')

    #generate a bam
    cmd = generateTempBamCmd(samtoolsString,fileNameDict)
    bashFile.write(cmd+'\n')

    #change into the temp directory
    cmd = changeTempDir(fileNameDict)
    bashFile.write(cmd+'\n')

    #sort the bam
    cmd = sortBamCmd(samtoolsString,fileNameDict)
    bashFile.write(cmd+'\n')

    #index
    cmd = indexBamCmd(samtoolsString,fileNameDict)
    bashFile.write(cmd+'\n')

    #remove sam
    if not options.sam:
        cmd = rmSamCmd(fileNameDict)
        bashFile.write(cmd+'\n')
    
    #or move the sam
    if options.sam:
        cmd = mvSamCmd(fileNameDict)
        bashFile.write(cmd+'\n')
    #mv bams
    cmd = mvBamCmd(fileNameDict)
    bashFile.write(cmd+'\n')

    #link bams
    if options.linkFolder:
        cmd = linkBamCmd(fileNameDict)
        bashFile.write(cmd+'\n')

    #cleanup
    cmd = rmTempFiles(fileNameDict)
    bashFile.write(cmd+'\n')


    bashFile.close()

    print "Wrote mapping command to %s" % (bashFileName)
if __name__ == "__main__":
    main()

