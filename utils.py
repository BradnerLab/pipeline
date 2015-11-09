#!/usr/bin/python
#SET OF GENERAL UTILITY FUNCTIONS FOR SEQ DATA
#last modified 141217

#please edit this to the location of the samtools program
samtoolsString ='samtools'

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


#Locus, LocusCollection, and Gene classes were generously provided by Graham Ruby and Charles Lin
#Additional functions are found from online sources and are noted in the comments



#==================================================================
#===========================DEPENDENCIES===========================
#==================================================================

import os
import gzip
import time
import re
import sys
import math

# Very pretty error reporting, where available
try:
    from IPython.core import ultratb
    sys.excepthook = ultratb.FormattedTB(mode='Context', color_scheme='Linux')
except ImportError:
    pass

from string import join, maketrans

import subprocess
import datetime

from collections import defaultdict

#==================================================================
#======================TABLE OF CONTENTS===========================
#==================================================================

#1. Input/Output and file handling functions

#def open(file,mode='r'):  <- replaces open with a version that can handle gzipped files
#def parseTable(fn, sep, header = False,excel = False): <- opens standard delimited files
#def unParseTable(table, output, sep): <- writes standard delimited files, opposite of parseTable
#def formatBed(bed,output=''):
#def bedToGFF(bed,output=''):
#def gffToBed(gff,output= ''): <- converts standard UCSC gff format files to UCSC bed format files
#def formatFolder(folderName,create=False): <- checks for the presence of any folder and makes it if create =True

#2. Gene annotation functions
#def makeStartDict(annotFile,geneList = []): <- takes a standard UCSC refseq table and creates a dictionary keyed by refseq ID with info about each transcript
#def getTSSs(geneList,refseqTable,refseqDict): <- returns the TSS location of any gene
#def importRefseq(refseqFile, returnMultiples = False): <- imports a standard UCSC refseq annotation file into a dictionary
#def makeGenes(annotFile,geneList=[],asDict = False): <- takes a UCSC refseq annotation file and a gene list and makes a list or dictionary of Gene class objects
#def makeTranscriptCollection(annotFile,upSearch,downSearch,window = 500,geneList = []): <- takes a UCSC refseq annotation file and makes a LocusCollection where each locus is a full transcript
#def importBoundRegion(boundRegionFile,name): <- imports a bound region file (a standard bed or macs output bed)

#3. Locus class
#class Locus(chr,start,end,sense,ID) <- standard locus class for tracking genomic loci
#class LocusCollection(lociList,windowSize=500) <- a collection of locus objects used for querying large sets of loci

#4. Gene class
#class Gene(name,chr,sense,txCoords,cdCoords,exStarts,exEnds,commonName=''): <- gene class object that contains all annotation information about a given transcript

#5. Locus functions
#def locusCollectionToGFF(locusCollection): <- turns a locus collection into a gff
#def gffToLocusCollection(gff,window =500): <- turns a gff into a locus collection (reverse of gff)
#def makeTSSLocus(gene,startDict,upstream,downstream): <- from a start dict makes a locus surrounding the tss
#def makeSearchLocus(locus,upSearch,downSearch): <- takes an existing locus and makes a larger flanking locus
#def makeSECollection(enhancerFile,name,top=0):


#6. Bam class
#class Bam(bamFile) <- a class for handling and manipulating bam objects.  requires samtools

#7. Misc. functions
#def uniquify(seq, idfun=None):  <- makes a list unique
#def order(x, NoneIsLast = True, decreasing = False): <- returns the ascending or descending order of a list


#8 Sequence functions
#def fetchSeq(directory,chrom,start,end,UCSC=False,lineBreaks=True,header = True): <- grabs sequence from a region
#def gffToFasta(species,directory,gff,UCSC = True): <- converts a gff to a fasta
#def revComp(seq,rev = True, RNA=False): <- is awesome


#==================================================================
#==========================I/O FUNCTIONS===========================
#==================================================================

# TODO: Overriding internal functions is evil!
bopen=open
def open(fileName,mode='r'):
    if fileName.split('.')[-1] == 'gz':
        return gzip.open(fileName, mode + 'b')
    else:
        return bopen(fileName, mode)

#parseTable 4/14/08
#takes in a table where columns are separated by a given symbol and outputs
#a nested list such that list[row][col]
#example call:
#table = parseTable('file.txt','\t')
def parseTable(fn, sep, header = False,excel = False):
    fh = open(fn)
    if header == True:
        header = fh.readline() #disposes of the header

    table = []
    for line in fh:
        line = line.rstrip().split(sep)
        table.append(line)

    fh.close()

    return table



#unParseTable 4/14/08
#takes in a table generated by parseTable and writes it to an output file
#takes as parameters (table, output, sep), where sep is how the file is delimited
#example call unParseTable(table, 'table.txt', '\t') for a tab del file

def unParseTable(table, output, sep):
    fh_out = open(output, 'w')
    if len(sep) == 0:
        for i in table:
            fh_out.write(str(i) + '\n')
    else:
        for line in table:
            fh_out.write(sep.join([str(x) for x in line]) + '\n')

    fh_out.close()




def formatBed(bed,output=''):

    '''
    formats a bed file from UCSC or MACS into a WUSTL gateway compatible bed
    '''


    newBed = []

    if type(bed) == str:
        bed = parseTable(bed,'\t')

    indexTicker = 1
    for line in bed:

        newLine = line[0:4]
        try:
            strand = line[5]
        except IndexError:
            strand = '.'
        newLine+= [indexTicker,strand]
        indexTicker +=1
        newBed.append(newLine)

    if len(output) > 0:
        unParseTable(newBed,output,'\t')
    else:
        return newBed

def bedToGFF(bed, output=''):

    '''
    turns a bed into a gff file
    '''
    if isinstance(bed, str):
        bed = parseTable(bed, '\t')

    bed = formatBed(bed)

    gff = []

    for line in bed:

        gffLine = [line[0],line[3],'',line[1],line[2],line[4],line[5],'',line[3]]
        gff.append(gffLine)


    if len(output) > 0:
        unParseTable(gff,output,'\t')
    else:
        return gff



#100912
#gffToBed

def gffToBed(gff,output= ''):
    '''
    turns a gff to a bed file
    '''
    bed = []
    for line in gff:
        newLine = [line[0],line[3],line[4],line[1],0,line[6]]
        bed.append(newLine)
    if len(output) == 0:
        return bed
    else:
        unParseTable(bed,output,'\t')

def formatFolder(folderName, create=False):

    '''
    makes sure a folder exists and if not makes it
    returns a bool for folder
    '''

    if folderName[-1] != '/':
        folderName +='/'

    try:
        os.listdir(folderName)
        return folderName
    except OSError:
        print('folder %s does not exist' % (folderName))
        if create:
            os.mkdir(folderName)
            return folderName
    return False


def checkOutput(fileName, waitTime = 1, timeOut = 30):

    '''
    checks for the presence of a file every N minutes
    if it exists, returns True
    default is 1 minute with a max timeOut of 30 minutes
    '''
    waitTime = int(waitTime*60)

    timeOut = int(timeOut*60)

    maxTicker = timeOut/waitTime
    ticker = 0

    fileExists = False
    while not fileExists:
        try:
            size1 = os.stat(fileName).st_size
            time.sleep(.5)
            size2 = os.stat(fileName).st_size
            if size1 == size2:
                fileExists = True
            else:
                time.sleep(waitTime)
                ticker+=1
        except OSError:
            time.sleep(waitTime)
            ticker+=1
        if ticker == maxTicker:
            break


    time.sleep(.1)
    if fileExists:
        return True
    else:
        print('OPERATION TIMED OUT. FILE %s NOT FOUND' % (fileName))
        return False

def getParentFolder(inputFile):

    '''
    returns the parent folder for any file
    '''

    parentFolder = join(inputFile.split('/')[:-1],'/') +'/'
    if parentFolder =='':
        return './'
    else:
        return parentFolder


#==================================================================
#===================ANNOTATION FUNCTIONS===========================
#==================================================================


def makeStartDict(annotFile, geneList=[]):
    '''
    makes a dictionary keyed by refseq ID that contains information about 
    chrom/start/stop/strand/common name
    '''

    if type(geneList) == str:
        geneList = parseTable(geneList, '\t')
        geneList = [line[0] for line in geneList]

    if annotFile.upper().count('REFSEQ') == 1:
        refseqTable,refseqDict = importRefseq(annotFile)
        if len(geneList) == 0:
            geneList = refseqDict.keys()
        startDict = {}
        for gene in geneList:
            if refseqDict.has_key(gene) == False:
                continue
            startDict[gene]={}
            startDict[gene]['sense'] = refseqTable[refseqDict[gene][0]][3]
            startDict[gene]['chr'] = refseqTable[refseqDict[gene][0]][2]
            startDict[gene]['start'] = getTSSs([gene],refseqTable,refseqDict)
            if startDict[gene]['sense'] == '+':
                startDict[gene]['end'] =[int(refseqTable[refseqDict[gene][0]][5])]
            else:
                startDict[gene]['end'] = [int(refseqTable[refseqDict[gene][0]][4])]
            startDict[gene]['name'] = refseqTable[refseqDict[gene][0]][12]
    return startDict


#generic function to get the TSS of any gene
def getTSSs(geneList,refseqTable,refseqDict):
    #refseqTable,refseqDict = importRefseq(refseqFile)
    if len(geneList) == 0:
        refseq = refseqTable
    else:
        refseq = refseqFromKey(geneList,refseqDict,refseqTable)
    TSS = []
    for line in refseq:
        if line[3] == '+':
            TSS.append(line[4])
        if line[3] == '-':
            TSS.append(line[5])
    TSS = map(int,TSS)

    return TSS




#10/13/08
#importRefseq
#takes in a refseq table and makes a refseq table and a refseq dictionary for keying the table

def importRefseq(refseqFile, returnMultiples = False):

    '''
    opens up a refseq file downloaded by UCSC
    '''
    refseqTable = parseTable(refseqFile,'\t')
    refseqDict = {}
    ticker = 1
    for line in refseqTable[1:]:
        if refseqDict.has_key(line[1]):
            refseqDict[line[1]].append(ticker)
        else:
            refseqDict[line[1]] = [ticker]
        ticker = ticker + 1

    multiples = []
    for i in refseqDict:
        if len(refseqDict[i]) > 1:
            multiples.append(i)

    if returnMultiples == True:
        return refseqTable,refseqDict,multiples
    else:
        return refseqTable,refseqDict


#12/29/08
#refseqFromKey(refseqKeyList,refseqDict,refseqTable)
#function that grabs refseq lines from refseq IDs
def refseqFromKey(refseqKeyList,refseqDict,refseqTable):
    typeRefseq = []
    for name in refseqKeyList:
        if refseqDict.has_key(name):
            typeRefseq.append(refseqTable[refseqDict[name][0]])
    return typeRefseq




#10/13/08
#make genes
#turns a refseq ID into a Gene object from utility module
def makeGenes(annotFile,geneList=[],asDict = False):
    '''
    takes in a refseq or ensembl annotation file and enters all identifiers in the geneList into a list as gene objects
    '''
    if asDict:
        genes = {}
    else:
        genes = []
    if type(geneList) == str:
        print('importing gene list from %s' % (geneList))
        geneList = parseTable(geneList,'\t')
        geneList = [line[0] for line in geneList]


    if annotFile.upper().count('REFSEQ') == 1:
        refTable,refDict = importRefseq(annotFile)

        if len(geneList) == 0:
            geneList = refDict.keys()

        for refseqID in geneList:
            if refseqID not in refDict:
                #print('no such gene ' + str(refseqID))
                continue

            geneIndex = refDict[refseqID][0]
            geneLine = refTable[int(geneIndex)]
            exonStarts = map(int,geneLine[9].split(',')[:-1])
            exonEnds = map(int,geneLine[10].split(',')[:-1])

            gene = Gene(refseqID,geneLine[2],geneLine[3],[int(geneLine[4]),int(geneLine[5])],[int(geneLine[6]),int(geneLine[7])],exonStarts,exonEnds,geneLine[12])

            if asDict:
                genes[refseqID] = gene
            else:
                genes.append(gene)

    return genes


#04/07/09
#makes a LocusCollection w/ each transcript as a locus
#bob = makeTranscriptCollection('/Users/chazlin/genomes/mm8/mm8refseq.txt')
def makeTranscriptCollection(annotFile,upSearch,downSearch,window = 500,geneList = []):
    '''
    makes a LocusCollection w/ each transcript as a locus
    takes in either a refseqfile or an ensemblGFF
    '''

    if annotFile.upper().count('REFSEQ') == 1:
        refseqTable,refseqDict = importRefseq(annotFile)
        locusList = []
        ticker = 0
        if len(geneList) == 0:
            geneList = refseqDict.keys()

        for line in refseqTable[1:]:
            if line[1] in geneList:
                if line[3] == '-':
                    locus = Locus(line[2],int(line[4])-downSearch,int(line[5])+upSearch,line[3],line[1])
                else:
                    locus = Locus(line[2],int(line[4])-upSearch,int(line[5])+downSearch,line[3],line[1])
                locusList.append(locus)
                ticker = ticker + 1
                if ticker%1000 == 0:
                    print(ticker)

    transCollection = LocusCollection(locusList, window)
    return transCollection


#140213

def nameToRefseq(geneNamesList,annotFile,unique=True):

    '''
    takes a list of names and gets you the refseqID
    '''
    startDict=  makeStartDict(annotFile)
    nameDict = defaultdict(list)
    for refID in startDict.keys():

        name = startDict[refID]['name']
        nameDict[name].append(refID)

    newTable = []
    for name in geneNamesList:
        refIDList = nameDict[name]

        #unique preserves the initial number of genes in geneNamesList
        #by taking only 1 refID per geneName
        #will take the first in the refList, which should usually be the lower
        #refID and thus the more relevant, but no guarantees
        if unique:
            newTable.append([name,refIDList[0]])
        else:
            for refID in refIDList:
                newTable.append([name,refID])

    return newTable


#06/11/09                                                                                                                              
#import bound region                                                                                                                   
#imports a bound region file and turns it into a locus collection                                                                      
#bound region files are output by my pipeline as Name_boundFile.txt files                                                              
def importBoundRegion(boundRegionFile,name):
    '''                                                                                                                                
    imports bound regions in either bed format or in error model format                                                                
    '''

    bound = parseTable(boundRegionFile,'\t')
    lociList = []
    ticker = 1
    #                                                                                                                                  
    if boundRegionFile.split('.')[-1] == 'bed':
        bed = True
    else:
        bed = False
    if bed:
        for line in bound:
            if ticker%1000 == 0:
                print(ticker)
            lociList.append(Locus(line[0],int(line[1]),int(line[2]),'.',ID = name + '_' + str(ticker)))
            ticker = ticker + 1
    else:
        for line in bound:
            if ticker%1000 == 0:
                print(ticker)

            lociList.append(Locus('chr'+line[0],int(line[1]),int(line[2]),'.',ID = name + '_' + str(ticker)))
            ticker = ticker + 1
    return LocusCollection(lociList,500)



#==================================================================
#========================LOCUS INSTANCE============================
#==================================================================

#Locus and LocusCollection instances courtesy of Graham Ruby


class Locus:
    # this may save some space by reducing the number of chromosome strings
    # that are associated with Locus instances (see __init__).
    __chrDict = dict()
    __senseDict = {'+':'+', '-':'-', '.':'.'}
    # chr = chromosome name (string)
    # sense = '+' or '-' (or '.' for an ambidexterous locus)
    # start,end = ints of the start and end coords of the locus;
    #      end coord is the coord of the last nucleotide.
    def __init__(self,chr,start,end,sense,ID='',score=0):
        coords = sorted([int(start), int(end)])
        # this method for assigning chromosome should help avoid storage of
        # redundant strings.
        if not(self.__chrDict.has_key(chr)): self.__chrDict[chr] = chr
        self._chr = self.__chrDict[chr]
        self._sense = self.__senseDict[sense]
        self._start = int(coords[0])
        self._end = int(coords[1])
        self._ID = ID
        self._score = score
    def ID(self): return self._ID
    def chr(self): return self._chr
    def start(self): return self._start  ## returns the smallest coordinate
    def end(self): return self._end   ## returns the biggest coordinate
    def len(self): return self._end - self._start + 1
    def score(self): return self._score
    def getAntisenseLocus(self):
        if self._sense=='.': return self
        else:
            switch = {'+':'-', '-':'+'}
            return Locus(self._chr,self._start,self._end,switch[self._sense])
    def coords(self): return [self._start,self._end]  ## returns a sorted list of the coordinates
    def sense(self): return self._sense
    # returns boolean; True if two loci share any coordinates in common
    def overlaps(self,otherLocus):
        if self.chr()!=otherLocus.chr(): return False
        elif not(self._sense=='.' or \
                             otherLocus.sense()=='.' or \
                             self.sense()==otherLocus.sense()): return False
        elif self.start() > otherLocus.end() or otherLocus.start() > self.end(): return False
        else: return True

    # returns boolean; True if all the nucleotides of the given locus overlap
    #      with the self locus
    def contains(self,otherLocus):
        if self.chr()!=otherLocus.chr(): return False
        elif not(self._sense=='.' or \
                             otherLocus.sense()=='.' or \
                             self.sense()==otherLocus.sense()): return False
        elif self.start() > otherLocus.start() or otherLocus.end() > self.end(): return False
        else: return True

    # same as overlaps, but considers the opposite strand
    def overlapsAntisense(self,otherLocus):
        return self.getAntisenseLocus().overlaps(otherLocus)
    # same as contains, but considers the opposite strand
    def containsAntisense(self,otherLocus):
        return self.getAntisenseLocus().contains(otherLocus)
    def __hash__(self): return self._start + self._end
    def __eq__(self,other):
        if self.__class__ != other.__class__: return False
        if self.chr()!=other.chr(): return False
        if self.start()!=other.start(): return False
        if self.end()!=other.end(): return False
        if self.sense()!=other.sense(): return False
        return True
    def __ne__(self,other): return not(self.__eq__(other))
    def __str__(self): return self.chr()+'('+self.sense()+'):'+'-'.join(map(str,self.coords()))
    def plotStr(self): return self.chr() + ':' + self.sense() + ':' + '-'.join(map(str,self.coords()))
    def checkRep(self):
        pass
    def gffLine(self): return [self.chr(),self.ID(),'',self.start(),self.end(),'',self.sense(),'',self.ID()]
    def getConservation(self,phastConFolder):
        '''
        uses tabix to get a per base conservation score from an indexed conservation bedgraph
        '''
        tabixString = 'tabix' #set the path/location of tabix

        #figure out which file is the correct one
        chrom = self.chr()

        phastFile = [phastConFolder + x for x in os.listdir(phastConFolder) if x.count('bg.gz') == 1 and x.count('tbi') == 0 and x.split('.')[0] == chrom][0]

        tabixCmd = '%s %s %s:%s-%s' % (tabixString,phastFile,self.chr(),self.start(),self.end())

        phast = subprocess.Popen(tabixCmd,stdin = subprocess.PIPE,stderr = subprocess.PIPE,stdout = subprocess.PIPE,shell = True)

        phastLines = phast.stdout.readlines()
        phast.stdout.close()

        phastTable = [x.rstrip().split('\t') for x in phastLines]

        #print phastTable
        #set up a conservation sum
        phastSum = 0.0
        #and number of bases w/ conservation data
        phastBases = 0

        for line in phastTable:
            if int(line[1]) < self.start():
                lineLen = min(int(line[2]),self.end()) - self.start() +1
                phastSum += lineLen*float(line[3])
                phastBases += lineLen
            elif int(line[2]) > self.end():
                lineLen = self.end() - max(int(line[1]),self.start())
                phastSum += lineLen*float(line[3])
                phastBases += lineLen
            else:
                lineLen = int(line[2]) - int(line[1])
                phastSum += lineLen*float(line[3])
                phastBases += lineLen

        if phastBases > self.len():
            print "this locus is sad %s. please debug me" % (self.__str__())
            print "locus length is %s" % (self.len())
            print "phastBases are %s" % (phastBases)


        return phastSum/self.len()


class LocusCollection:
    def __init__(self,loci,windowSize=50):
        ### top-level keys are chr, then strand, no space
        self.__chrToCoordToLoci = dict()
        self.__loci = dict()
        self.__winSize = windowSize
        for lcs in loci: self.__addLocus(lcs)

    def __addLocus(self,lcs):
        if not(self.__loci.has_key(lcs)):
            self.__loci[lcs] = None
            if lcs.sense()=='.': chrKeyList = [lcs.chr()+'+', lcs.chr()+'-']
            else: chrKeyList = [lcs.chr()+lcs.sense()]
            for chrKey in chrKeyList:
                if not(self.__chrToCoordToLoci.has_key(chrKey)): self.__chrToCoordToLoci[chrKey] = dict()
                for n in self.__getKeyRange(lcs):
                    if not(self.__chrToCoordToLoci[chrKey].has_key(n)): self.__chrToCoordToLoci[chrKey][n] = []
                    self.__chrToCoordToLoci[chrKey][n].append(lcs)

    def __getKeyRange(self,locus):
        start = locus.start() / self.__winSize
        end = locus.end() / self.__winSize + 1 ## add 1 because of the range
        return range(start,end)

    def __len__(self): return len(self.__loci)

    def append(self,new): self.__addLocus(new)
    def extend(self,newList):
        for lcs in newList: self.__addLocus(lcs)
    def hasLocus(self,locus):
        return self.__loci.has_key(locus)
    def remove(self,old):
        if not(self.__loci.has_key(old)): raise ValueError("requested locus isn't in collection")
        del self.__loci[old]
        if old.sense()=='.': senseList = ['+','-']
        else: senseList = [old.sense()]
        for k in self.__getKeyRange(old):
            for sense in senseList:
                self.__chrToCoordToLoci[old.chr()+sense][k].remove(old)

    def getWindowSize(self): return self.__winSize
    def getLoci(self): return self.__loci.keys()

    def getSize(self):
           size = 0
           for locus in self.__loci:
                  newsize = int(locus.end())-int(locus.start())
                  size += newsize
           return size

    def getChrList(self):
        # i need to remove the strand info from the chromosome keys and make
        # them non-redundant.
        tempKeys = dict()
        for k in self.__chrToCoordToLoci.keys(): tempKeys[k[:-1]] = None
        return tempKeys.keys()

    def __subsetHelper(self,locus,sense):
        sense = sense.lower()
        if ['sense','antisense','both'].count(sense)!=1:
            raise ValueError("sense command invalid: '"+sense+"'.")
        matches = dict()
        senses = ['+','-']
        if locus.sense()=='.' or sense=='both': lamb = lambda s: True
        elif sense=='sense': lamb = lambda s: s==locus.sense()
        elif sense=='antisense': lamb = lambda s: s!=locus.sense()
        else: raise ValueError("sense value was inappropriate: '"+sense+"'.")
        for s in filter(lamb, senses):
            chrKey = locus.chr()+s
            if self.__chrToCoordToLoci.has_key(chrKey):
                for n in self.__getKeyRange(locus):
                    if self.__chrToCoordToLoci[chrKey].has_key(n):
                        for lcs in self.__chrToCoordToLoci[chrKey][n]:
                            matches[lcs] = None
        return matches.keys()

    # sense can be 'sense' (default), 'antisense', or 'both'
    # returns all members of the collection that overlap the locus
    def getOverlap(self,locus,sense='sense'):
        matches = self.__subsetHelper(locus,sense)
        ### now, get rid of the ones that don't really overlap
        realMatches = dict()
        if sense=='sense' or sense=='both':
            for i in filter(lambda lcs: lcs.overlaps(locus), matches):
                realMatches[i] = None
        if sense=='antisense' or sense=='both':
            for i in filter(lambda lcs: lcs.overlapsAntisense(locus), matches):
                realMatches[i] = None
        return realMatches.keys()

    # sense can be 'sense' (default), 'antisense', or 'both'
    # returns all members of the collection that are contained by the locus
    def getContained(self,locus,sense='sense'):
        matches = self.__subsetHelper(locus,sense)
        ### now, get rid of the ones that don't really overlap
        realMatches = dict()
        if sense=='sense' or sense=='both':
            for i in filter(lambda lcs: locus.contains(lcs), matches):
                realMatches[i] = None
        if sense=='antisense' or sense=='both':
            for i in filter(lambda lcs: locus.containsAntisense(lcs), matches):
                realMatches[i] = None
        return realMatches.keys()

    # sense can be 'sense' (default), 'antisense', or 'both'
    # returns all members of the collection that contain the locus
    def getContainers(self,locus,sense='sense'):
        matches = self.__subsetHelper(locus,sense)
        ### now, get rid of the ones that don't really overlap
        realMatches = dict()
        if sense=='sense' or sense=='both':
            for i in filter(lambda lcs: lcs.contains(locus), matches):
                realMatches[i] = None
        if sense=='antisense' or sense=='both':
            for i in filter(lambda lcs: lcs.containsAntisense(locus), matches):
                realMatches[i] = None
        return realMatches.keys()

    def stitchCollection(self,stitchWindow=1,sense='both'):

        '''
        reduces the collection by stitching together overlapping loci
        returns a new collection
        '''

        #initializing stitchWindow to 1 
        #this helps collect directly adjacent loci



        locusList = self.getLoci()
        oldCollection = LocusCollection(locusList,500)

        stitchedCollection = LocusCollection([],500)

        for locus in locusList:
            #print(locus.coords())
            if oldCollection.hasLocus(locus):
                oldCollection.remove(locus)
                overlappingLoci = oldCollection.getOverlap(Locus(locus.chr(),locus.start()-stitchWindow,locus.end()+stitchWindow,locus.sense(),locus.ID()),sense)

                stitchTicker = 1
                while len(overlappingLoci) > 0:
                    stitchTicker+=len(overlappingLoci)
                    overlapCoords = locus.coords()


                    for overlappingLocus in overlappingLoci:
                        overlapCoords+=overlappingLocus.coords()
                        oldCollection.remove(overlappingLocus)
                    overlapCoords = [int(x) for x in overlapCoords]
                    if sense == 'both':
                        locus = Locus(locus.chr(),min(overlapCoords),max(overlapCoords),'.',locus.ID())
                    else:
                        locus = Locus(locus.chr(),min(overlapCoords),max(overlapCoords),locus.sense(),locus.ID())
                    overlappingLoci = oldCollection.getOverlap(Locus(locus.chr(),locus.start()-stitchWindow,locus.end()+stitchWindow,locus.sense()),sense)
                locus._ID = '%s_%s_lociStitched' % (stitchTicker,locus.ID())

                stitchedCollection.append(locus)

            else:
                continue
        return stitchedCollection



#==================================================================
#========================GENE INSTANCE============================
#==================================================================
# this is a gene object.  unlike the previous gene_object, this actually represents
# a gene, as opposed to a whole set of genes, which was a poor design in the first place.
class Gene:
    # name = name of the gene (string)
    # txCoords = list of coords defining the boundaries of the transcipt
    # cdCoords = list of coords defining the beginning and end of the coding region
    # exStarts = list of coords marking the beginning of each exon
    # exEnds = list of coords marking the end of each exon
    # IF THIS IS A NON-CODING GENE, cdCoords => [0,0]
    #    def __init__(self,name,chr,sense,txCoords,cdCoords,exStarts,exEnds):
    #        self._name = name
    #        self._txLocus = Locus(chr,min(txCoords),max(txCoords),sense)
    #        self._cdLocus = Locus(chr,min(cdCoords),max(cdCoords),sense)
    #
    #        exStarts = map(lambda i: i, exStarts)
    #        exEnds = map(lambda i: i, exEnds)
    #        exStarts.sort()
    #        exEnds.sort()
    #
    #        self._txExons = []
    #        self._cdExons = []
    #        self._introns = []
    #
    #        for n in range(len(exStarts)):
    #            if n==0:
    #                self._txExons.append(Locus(chr,txCoords[0],exEnds[n]-1,sense))
    #                self._cdExons.append(Locus(chr,cdCoords[0],exEnds[n]-1,sense))
    #            elif n==len(exStarts)-1:
    #                self._txExons.append(Locus(chr,txCoords[0],txCoords[1],sense))
    #                self._cdExons.append(Locus(chr,cdCoords[0],cdCoords[1],sense))
    #            else:
    #                newExon = Locus(chr,exStarts[n],exEnds[n]-1,sense)
    #                self._txExons.append(newExon)
    #                self._cdExons.append(newExon)
    #            if n < len(exStarts)-1: self._introns.append(Locus(chr,exEnds[n],exStarts[n+1]-1,sense))
    #
    #        if sense=='+':
    #            self._fpUtr = Locus(chr,txCoords[0],cdCoords[0]-1,sense)
    #            self._tpUtr = Locus(chr,cdCoords[1]+1,txCoords[1],sense)
    #        elif sense=='-':
    #            self._fpUtr = Locus(chr,cdCoords[1]+1,txCoords[1],sense)
    #            self._tpUtr = Locus(chr,txCoords[0],cdCoords[0]-1,sense)
    def __init__(self,name,chr,sense,txCoords,cdCoords,exStarts,exEnds,commonName=''):
        self._name = name
        self._commonName = commonName
        self._txLocus = Locus(chr,min(txCoords),max(txCoords),sense,self._name)
        if cdCoords == None:
            self._cdLocus = None
        else:
            self._cdLocus = Locus(chr,min(cdCoords),max(cdCoords),sense)

        exStarts = sorted(map(lambda i: i, exStarts))
        exEnds = sorted(map(lambda i: i, exEnds))

        self._txExons = []
        self._cdExons = []
        self._introns = []

        cd_exon_count = 0

        for n in range(len(exStarts)):
            first_locus = Locus(chr,exStarts[n],exStarts[n],sense)
            second_locus = Locus(chr,exEnds[n],exEnds[n],sense)

            # Add the transcription unit exon
            tx_exon = Locus(chr,exStarts[n],exEnds[n],sense)

            self._txExons.append(tx_exon)

            # Add Coding Exons
            # Need to make sure that the current exon is actually in the coding region of the gene first
            if self.isCoding() and tx_exon.overlaps(self._cdLocus):
                if not first_locus.overlaps(self._cdLocus):
                    first_coord = min(cdCoords)
                else:
                    first_coord = exStarts[n]

                if not second_locus.overlaps(self._cdLocus):
                    second_coord = max(cdCoords)
                else:
                    second_coord = exEnds[n]

                new_cd_exon = Locus(chr,first_coord,second_coord,sense)
                self._cdExons.append(new_cd_exon)

            # Add Introns
            if n < len(exStarts)-1:
                self._introns.append(Locus(chr,exEnds[n]+1,exStarts[n +1]-1,sense))

        if self.isCoding():
            if sense=='+':
                self._fpUTR = Locus(chr,min(txCoords),min(cdCoords)-1,sense)
                self._tpUTR = Locus(chr,max(cdCoords)+1,max(txCoords),sense)
            elif sense=='-':
                self._fpUTR = Locus(chr,max(cdCoords)+1,max(txCoords),sense)
                self._tpUTR = Locus(chr,min(txCoords),min(cdCoords)-1,sense)
        else:
            self._fpUTR = None
            self._tpUTR = None

    def commonName(self): return self._commonName
    def name(self): return self._name
    def chr(self): return self._txLocus.chr()
    def sense(self): return self._txLocus.sense()
    def txLocus(self): return self._txLocus   ## locus of full transcript
    def cdLocus(self): return self._cdLocus   ## locus from start codon to end codon
    def txExons(self): return map(lambda i: i, self._txExons)  ## list of loci
    def cdExons(self): return map(lambda i: i, self._cdExons)  ## list of loci
    def introns(self): return map(lambda i: i, self._introns)  ## list of loci
    def fpUtr(self): return self._fpUTR  ## locus
    def tpUtr(self): return self._tpUTR  ## locus
    def isCoding(self): return not(self._cdLocus.start()==0 and self._cdLocus.end()==0)  # boolean; is this gene protein-coding?
    def tss(self,upstream = 0,downstream = 0):
        if self._txLocus.sense() == '-':
            return Locus(self._txLocus.chr(),self._txLocus.end()-downstream,self._txLocus.end()+upstream,self._txLocus.sense(),self._name)
        else:
            return Locus(self._txLocus.chr(),self._txLocus.start()-upstream,self._txLocus.start()+downstream,self._txLocus.sense(),self._name)
    def __hash__(self): return self._txLocus.__hash__()



#==================================================================
#========================LOCUS FUNCTIONS===========================
#==================================================================
#06/11/09
#turns a locusCollection into a gff
#does not write to disk though


def locusCollectionToBed(locusCollection):

    lociList = locusCollection.getLoci()
    bed = []
    for locus in lociList:
        coords = locus.coords()
        coords.sort()
        newLine = [locus.chr(),coords[0],coords[1],locus.sense(),locus.ID()]
        bed.append(newLine)
    return bed



def locusCollectionToGFF(locusCollection):
    lociList = locusCollection.getLoci()
    gff = []
    for locus in lociList:
        newLine = [locus.chr(),locus.ID(),'',locus.coords()[0],locus.coords()[1],'',locus.sense(),'',locus.ID()]
        gff.append(newLine)
    return gff



def bedToLocusCollection(bedfile):

    table = parseTable(bedfile, '\t')
    loci = [Locus(x[0], x[1], x[2], '.', ID=x[3]) for x in table]
    collection = LocusCollection(loci, 50)

    return collection

def gffToLocusCollection(gff,window =500):

    '''
    opens up a gff file and turns it into a LocusCollection instance
    '''

    lociList = []
    if type(gff) == str:
        gff = parseTable(gff,'\t')

    for line in gff:
        #USE line[1] as the locus ID.  If that is empty use line[8]
        if len(line[1]) > 0:
            name = line[1]
        elif len(line[8]) >0:
            name = line[8]
        else:
            name = '%s:%s:%s-%s' % (line[0],line[6],line[3],line[4])

        lociList.append(Locus(line[0],line[3],line[4],line[6],name))
    return LocusCollection(lociList,window)



def makeTSSLocus(gene,startDict,upstream,downstream):
    '''
    given a startDict, make a locus for any gene's TSS w/ upstream and downstream windows
    '''

    start = startDict[gene]['start'][0]
    if startDict[gene]['sense'] == '-':
        return Locus(startDict[gene]['chr'],start-downstream,start+upstream,'-',gene)
    else:
        return Locus(startDict[gene]['chr'],start-upstream,start+downstream,'+',gene)


#06/11/09
#takes a locus and expands it by a fixed upstream/downstream amount. spits out the new larger locus
def makeSearchLocus(locus,upSearch,downSearch):
    if locus.sense() == '-':
        searchLocus = Locus(locus.chr(),locus.start()-downSearch,locus.end()+upSearch,locus.sense(),locus.ID())
    else:
        searchLocus = Locus(locus.chr(),locus.start()-upSearch,locus.end()+downSearch,locus.sense(),locus.ID())
    return searchLocus


def makeSECollection(enhancerFile,name,top=0):
    '''
    returns a locus collection from a super table
    top gives the number of rows
    '''
    enhancerTable = parseTable(enhancerFile,'\t')
    superLoci = []

    ticker = 0
    for line in enhancerTable:
        if line[0][0] == '#' or line[0][0] == 'R':
            continue
        else:
            ticker+=1

            superLoci.append(Locus(line[1],line[2],line[3],'.',name+'_'+line[0]))

            if ticker == top:
                break
    return LocusCollection(superLoci,50)


#==================================================================
#==========================BAM CLASS===============================
#==================================================================

#11/11/10
#makes a new class Bam for dealing with bam files and integrating them into the SolexaRun class

def convertBitwiseFlag(flag):
    if int(flag) & 16:
        return "-";
    else:
        return "+";

class Bam:
    '''A class for a sorted and indexed bam file that allows easy analysis of reads'''
    def __init__(self,bamFile):
        self._bam = bamFile

    def getTotalReads(self,readType = 'mapped'):
        command = '%s flagstat %s' % (samtoolsString,self._bam)
        stats = subprocess.Popen(command,stdin = subprocess.PIPE,stderr = subprocess.PIPE,stdout = subprocess.PIPE,shell = True)
        statLines = stats.stdout.readlines()
        stats.stdout.close()
        if readType == 'mapped':
            for line in statLines:
                if line.count('mapped (') == 1:

                    return int(line.split(' ')[0])
        if readType == 'total':
            return int(statLines[0].split(' ')[0])

    def convertBitwiseFlag(self,flag):
        if flag & 16:
            return "-";
        else:
            return "+";

    def getRawReads(self,locus,sense,unique = False,includeJxnReads = False,printCommand = False):
        '''
        gets raw reads from the bam using samtools view.
        can enforce uniqueness and strandedness
        '''
        locusLine = locus.chr()+':'+str(locus.start())+'-'+str(locus.end())

        command = '%s view %s %s' % (samtoolsString,self._bam,locusLine)
        if printCommand:
            print(command)
        getReads = subprocess.Popen(command,stdin = subprocess.PIPE,stderr = subprocess.PIPE,stdout = subprocess.PIPE,shell = True)
        reads = getReads.communicate()
        reads = reads[0].split('\n')[:-1]
        reads = [read.split('\t') for read in reads]
        if includeJxnReads == False:
            reads = filter(lambda x: x[5].count('N') < 1,reads)

        #convertDict = {'16':'-','0':'+','64':'+','65':'+','80':'-','81':'-','129':'+','145':'-'}
        #convertDict = {'16':'-','0':'+','64':'+','65':'+','80':'-','81':'-','129':'+','145':'-','256':'+','272':'-','99':'+','147':'-'}


        #BJA added 256 and 272, which correspond to 0 and 16 for multi-mapped reads respectively:
        #http://onetipperday.blogspot.com/2012/04/understand-flag-code-of-sam-format.html
        #convert = string.maketrans('160','--+')
        keptReads = []
        seqDict = defaultdict(int)
        if sense == '-':
            strand = ['+','-']
            strand.remove(locus.sense())
            strand = strand[0]
        else:
            strand = locus.sense()
        for read in reads:
            #readStrand = read[1].translate(convert)[0]
            #print read[1], read[0]
            #readStrand = convertDict[read[1]]
            readStrand = convertBitwiseFlag(read[1])

            if sense == 'both' or sense == '.' or readStrand == strand:

                if unique and seqDict[read[9]] == 0:
                    keptReads.append(read)
                elif not unique:
                    keptReads.append(read)
            seqDict[read[9]]+=1

        return keptReads

    def readsToLoci(self,reads,IDtag = 'sequence,seqID,none'):
        '''
        takes raw read lines from the bam and converts them into loci
        '''
        loci = []
        ID = ''
        if IDtag == 'sequence,seqID,none':
            print('please specify one of the three options: sequence, seqID, none')
            return
        #convert = string.maketrans('160','--+')
        #convertDict = {'16':'-','0':'+','64':'+','65':'+','80':'-','81':'-','129':'+','145':'-'}
        #convertDict = {'16':'-','0':'+','64':'+','65':'+','80':'-','81':'-','129':'+','145':'-','256':'+','272':'-'}

        #BJA added 256 and 272, which correspond to 0 and 16 for multi-mapped reads respectively:
        #http://onetipperday.blogspot.com/2012/04/understand-flag-code-of-sam-format.html
        #convert = string.maketrans('160','--+')
        numPattern = re.compile('\d*')
        for read in reads:
            chrom = read[2]
            #strand = read[1].translate(convert)[0]
            #strand = convertDict[read[1]]
            strand = convertBitwiseFlag(read[1])
            if IDtag == 'sequence':
                ID = read[9]
            elif IDtag == 'seqID':
                ID = read[0]
            else:
                ID = ''

            length = len(read[9])
            start = int(read[3])
            if read[5].count('N') == 1:
                #this awful oneliner first finds all of the numbers in the read string
                #then it filters out the '' and converts them to integers
                #only works for reads that span one junction

                [first,gap,second] = [int(x) for x in filter(lambda x: len(x) > 0, re.findall(numPattern,read[5]))][0:3]
                if IDtag == 'sequence':
                    loci.append(Locus(chrom,start,start+first,strand,ID[0:first]))
                    loci.append(Locus(chrom,start+first+gap,start+first+gap+second,strand,ID[first:]))
                else:
                    loci.append(Locus(chrom,start,start+first,strand,ID))
                    loci.append(Locus(chrom,start+first+gap,start+first+gap+second,strand,ID))
            elif read[5].count('N') > 1:
                continue
            else:
                loci.append(Locus(chrom,start,start+length,strand,ID))
        return loci

    def getReadsLocus(self,locus,sense = 'both',unique = True,IDtag = 'sequence,seqID,none',includeJxnReads = False):
        '''
        gets all of the reads for a given locus
        '''
        reads = self.getRawReads(locus,sense,unique,includeJxnReads)

        loci = self.readsToLoci(reads,IDtag)

        return loci

    def getReadSequences(self,locus,sense = 'both',unique = True,includeJxnReads = False):

        reads = self.getRawReads(locus,sense,unique,includeJxnReads)

        return [read[9] for read in reads]

    def getReadStarts(self,locus,sense = 'both',unique = False,includeJxnReads = False):
        reads = self.getRawReads(locus,sense,unique,includeJxnReads)

        return [int(read[3]) for read in reads]


    def getReadCount(self,locus,sense = 'both',unique = True,includeJxnReads = False):
        reads = self.getRawReads(locus,sense,unique,includeJxnReads)

        return len(reads)

    def liquidateLocus(self,locus,sense='.'):

           bamliquidatorCmd = 'bamliquidator %s %s %s %s %s 1 200' % (self._bam, locus.chr(),
                                                                     str(locus.start()), str(locus.end()),
                                                                     sense)

           bamliquidatorOut = subprocess.Popen(bamliquidatorCmd, stdout = subprocess.PIPE, shell=True)
           score = bamliquidatorOut.communicate()[0]

           return score

#==================================================================
#========================MISC FUNCTIONS============================
#==================================================================




#uniquify function
#by Peter Bengtsson
#Used under a creative commons license
#sourced from  here: http://www.peterbe.com/plog/uniqifiers-benchmark

def uniquify(seq, idfun=None):
    # order preserving
    if idfun is None:
        def idfun(x): return x
    seen = {}
    result = []
    for item in seq:
        marker = idfun(item)
        # in old Python versions:
        # if seen.has_key(marker)
        # but in new ones:
        if marker in seen: continue
        seen[marker] = 1
        result.append(item)
    return result


#082009
#taken from http://code.activestate.com/recipes/491268/

def order(x, NoneIsLast = True, decreasing = False):
    """
    Returns the ordering of the elements of x. The list
    [ x[j] for j in order(x) ] is a sorted version of x.

    Missing values in x are indicated by None. If NoneIsLast is true,
    then missing values are ordered to be at the end.
    Otherwise, they are ordered at the beginning.
    """
    omitNone = False
    if NoneIsLast == None:
        NoneIsLast = True
        omitNone = True

    n  = len(x)
    ix = range(n)
    if None not in x:
        ix.sort(reverse = decreasing, key = lambda j : x[j])
    else:
        # Handle None values properly.
        def key(i, x = x):
            elem = x[i]
            # Valid values are True or False only.
            if decreasing == NoneIsLast:
                return not(elem is None), elem
            else:
                return elem is None, elem
        ix = range(n)
        ix.sort(key=key, reverse=decreasing)

    if omitNone:
        n = len(x)
        for i in range(n-1, -1, -1):
            if x[ix[i]] == None:
                n -= 1
        return ix[:n]
    return ix



#==================================================================
#=======================SEQUENCE FUNCTIONS=========================
#==================================================================



#10/22/08
#fetchSeq
#function that fetches a sequence from a genome directory
#directory that contains individual chrom fasta files

def fetchSeq(directory,chrom,start,end,UCSC=False,lineBreaks=True,header = True):
    fn = directory + chrom + '.fa'
    fh = open(fn,'r')
    headerOffset = 0
    nStart = 0
    nEnd = 0
    if header:
        fh.seek(0)
        headerOffset = len(fh.readline())
    if lineBreaks:

        nStart = (start-1)/50
        nEnd = (end-1)/50
    if UCSC:
        fh.seek((start+nStart+headerOffset))
    else:
        fh.seek((start-1+nStart+headerOffset))
    span = ((end+nEnd-1)-(start+nStart-1))
    # if UCSC:
    #     span+=1

    read = fh.read(span)
    if lineBreaks:
        read = read.replace('\n','')
    #print(headerOffset,nStart,nEnd)
    return read
    fh.close()



#10/22/08
#gffToFasta
#function that writes a fasta file from a gff file
#directory is the genome directory with the chromosome folders
def gffToFasta(genome,directory,gff,UCSC = True,useID=False):
    fastaList = []

    ticker = 0
    if type(gff) == str:
        gff = parseTable(gff,'\t')
    for line in gff:
        try:
            sequence = fetchSeq(directory,line[0],int(line[3]),int(line[4]),UCSC)
        except:
            continue
        if ticker%1000 == 0: print(ticker)

        if useID:
            name = '>' + line[1]
        else:
            name = '>'+ join([genome.lower(),line[0],str(line[3]),str(line[4]),line[6]],'|')
        fastaList.append(name)
        if line[6] == '-':
            #print(line[3])
            #print(line[4])
            fastaList.append(revComp(sequence))
        else:
            fastaList.append(sequence)
        ticker += 1
    return fastaList



#10/23/08
#pair(nuc)
#returns the basepair of a nucleotide

def pair(nuc):
    pairDict = {'A':'T','C':'G','G':'C','T':'A','U':'A','a':'t','c':'g','g':'c','t':'a','u':'a'}
    if pairDict.has_key(nuc):
        return pairDict[nuc]
    else:
        return nuc

#10/23/08
#revComp(seq)
#takes in a sequence and spits out the reverse Complement

def revComp(seq,rev = True, RNA=False):
    if rev:
        revComp = join(map(pair,seq[::-1]),'')
    else:
        revComp = join(map(pair,seq),'')
    if RNA:
        revComp = revComp.replace('T','U') # Not sure this is defined (!)
        revComp = revComp.replace('t','u')
    else:
        revComp = revComp.replace('U','T')
        revComp = revComp.replace('u','t')
    return revComp
