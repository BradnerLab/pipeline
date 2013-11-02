#!/usr/bin/python
#GPL16043.py

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

#scripts and packages for dealing with the human primeview array

#==================================================================
#=========================MODULES==================================
#==================================================================

import string
from collections import defaultdict
#==================================================================
#==================HELPER FUNCTIONS================================
#==================================================================

def parseTable(fn, sep, header = False,excel = False):

    '''
    table parser helper function
    '''
    fh = open(fn)
    lines = fh.readlines()
    fh.close()
    if excel:
        lines = lines[0].split('\r')
    if lines[0].count('\r') > 0:
        lines = lines[0].split('\r')
    table = []
    if header == True:
        lines =lines[1:]
    for i in lines:
        table.append(i[:-1].split(sep))

    return table

def unParseTable(table, output, sep):
    fh_out = open(output,'w')
    if len(sep) == 0:
        for i in table:
            fh_out.write(str(i))
            fh_out.write('\n')
    else:
        for line in table:
            line = [str(x) for x in line]
            line = string.join(line,sep)

            fh_out.write(line)
            fh_out.write('\n')

    fh_out.close()


def importRefseq(refseqFile, returnMultiples = False):
    '''
    parser for standard refseq tables downloaded from UCSC
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

def getTSSs(geneList,refseqTable,refseqDict):
    '''
    returns all TSSs for each gene in a geneList
    good for finding multiple starts
    '''
    def refseqFromKey(refseqKeyList,refseqDict,refseqTable):
        typeRefseq = []
        for name in refseqKeyList:
            if refseqDict.has_key(name):
                typeRefseq.append(refseqTable[refseqDict[name][0]])
        return typeRefseq

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


def makeStartDict(annotFile,geneList = []):
    '''makes a dictionary of gene coordinates and names'''
    if type(geneList) == str:
        geneList = parseTable(geneList,'\t')
        geneList = [line[0] for line in geneList]
            
    if string.upper(annotFile).count('REFSEQ') == 1:
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


probeAnnotationFile  = '/nfs/young_ata4/gff/tables/platforms/GPL16043_probe_annotation.txt'

annotFile = '/nfs/young_ata4/gff/annotation/hg18/hg18_refseq.ucsc'

probeExpressionFile = '/nfs/young_ata4/pipeline/mm1s_new/expression/11162012_mas5_allProbes.txt'

output = '/nfs/young_ata4/pipeline/mm1s_new/expression/11162012_mas5_allGenes.txt'


#==================================================================
#==================GENE TABLE FUNCTION=============================
#==================================================================


def makeGeneTable(annotFile,probeAnnotationFile,probeExpressionFile,output):

    startDict= makeStartDict(annotFile)
    refIDList = startDict.keys()
    refIDList.sort()

    geneNameList = [startDict[x]['name'] for x in refIDList]

    probeAnnotationTable = parseTable(probeAnnotationFile,'\t')
    #make a refseq based probe dict
    probeDict = defaultdict(list)
    print('MAKING PROBE ANNOTATION DICTIONARY\nNUMBER OF PROBES PROCESSED')
    ticker = 0
    for line in probeAnnotationTable[1:]:
        if ticker%10000 == 0:
            print(ticker)
        ticker+=1
        try:
            geneID = line[1]
        except IndexError:
            #print(line)
            continue
        #if line[0] == '11759571_a_at':
        #    print(line)
        if line[1][0:2] == 'NM':
            geneID = line[1].split('.')[0]
            probeDict[line[0]] = [geneID]
        else:
            geneName = line[3]
            indices = [i for i, x in enumerate(geneNameList) if x == geneName]
            probeDict[line[0]] = [refIDList[index] for index in indices]



    print('FOUND %s PROBES CORRESPONDING TO REFSEQ GENES' % (len(probeDict.keys())))
    expDict = defaultdict(list)

    probeExpressionTable = parseTable(probeExpressionFile,'\t')

    if len(probeExpressionTable[0]) == len(probeExpressionTable[1]):
        sampleHeader = probeExpressionTable[0][1:]
    else:
        sampleHeader = probeExpressionTable[0]
    print('MAKING EXPRESSION DICT')
    for line in probeExpressionTable[1:]:
        probeID = line[0]
        if probeDict.has_key(probeID):
            geneIDList = probeDict[probeID]
            for geneID in geneIDList:
                expDict[geneID].append([float(x) for x in line[1:]])
            
    print('FOUND EXPRESSION VALUES FOR %s GENES' % (len(expDict.keys())))

    #now for each gene, take the expression probe with the highest average levels

    
    expTable =[['GENE_ID','NAME']+sampleHeader]

    geneList = expDict.keys()
    geneList.sort()

    print('MAKING GENE EXPRESSION TABLE\nNUMBER OF GENES PROCESSED')
    ticker = 0
    for geneID in geneList:
        if ticker % 10000 == 0:
            print(ticker)
        ticker+=1
        meanVector = [float(sum(x))/len(x) for x in expDict[geneID]]
        
        maxIndex = meanVector.index(max(meanVector))

        expression = expDict[geneID][maxIndex]

        try:
            geneName = startDict[geneID]['name']
        except KeyError:
            continue
        expTable.append([geneID,geneName]+expression)

    print('MAKING EXPRESSION TABLE OF %s GENES WITH CORRECT ANNOTATION OUT OF %s TOTAL PROBES IN ARRAY AND WRITING OUTPUT TO %s' % (len(expTable) -1,len(probeDict.keys()),output))
    unParseTable(expTable,output,'\t')

#==================================================================
#======================MAIN METHOD=================================
#==================================================================

def main():

    '''
    main run function
    '''
    from optparse import OptionParser

    usage = "usage: %prog [options] -a [REFSEQ_ANNOTATION] -p [PROBE_ANNOTATION] -i [PROBE_EXPRESSION_TABLE]"
    parser = OptionParser(usage = usage)
    #required flags
    parser.add_option("-a","--annot", dest="annot",nargs = 1, 
                      default='/ark/home/cl512/ressrv19/annotations/hg18_refseq.ucsc',
                      help = "Enter a refseq annotation file to be used. Default is /ark/home/cl512/ressrv19/annotations/hg18_refseq.ucsc")
    parser.add_option("-p","--probe", dest="probe",nargs = 1, 
                      default='/ark/home/cl512/ressrv19/annotations/platforms/GPL16043/GPL16043_probe_annotation.txt', 
                      help = "Enter a refseq annotation file to be used. Default is /ark/home/cl512/ressrv19/annotations/GPL16043/GPL16043_probe_annotation.txt")
    parser.add_option("-i","--input", dest="input",nargs = 1, 
                      default=None, 
                      help = "Enter a expression table indexed by probe ID")
    #optional flags
    parser.add_option("-o","--output", dest="output",nargs = 1, 
                      default=None, 
                      help = "Enter a folder to store output. Default will store output in same folder as input expression table")

    (options,args) = parser.parse_args()

    #print(options)
    #print(args)
    
    if options.input:

        #organize the input
        annotFile = options.annot
        probeAnnotationFile = options.probe
        probeExpressionFile = options.input
        
        if options.output:
            outputFolder = options.output
        else:
            #use the same folder as expression input
            outputFolder = string.join(probeExpressionFile.split('/')[0:-1],'/')

        #make sure outputFolder has a '/' at the end
        if outputFolder[-1] != '/':
            outputFolder += '/'
        #get the name of the input
        #input file name will always be in the form of [foldername]_mas5_spikeNorm_allProbes.txt
        inputName = probeExpressionFile.split('/')[-1]
        outputName = inputName.replace('Probes','Genes')
        outputName = inputName.replace('probe','gene')
        outputFile = outputFolder + outputName
        print('GETTING EXPRESSION DATA FROM %s\n AND WRITING OUTPUT TO %s' % (probeExpressionFile,outputFile))


        makeGeneTable(annotFile,probeAnnotationFile,probeExpressionFile,outputFile)



    else:
        parser.print_help()
        exit()

if __name__ == "__main__":
    main()

# #=============================
# #MM1S MAS5 raw
# probeAnnotationFile  = '/nfs/young_ata4/gff/tables/platforms/GPL16043_probe_annotation.txt'

# annotFile = '/nfs/young_ata4/gff/annotation/hg18/hg18_refseq.ucsc'

# probeExpressionFile = '/nfs/young_ata4/pipeline/mm1s_new/expression/11162012_mas5_allProbes.txt'

# output = '/nfs/young_ata4/pipeline/mm1s_new/expression/11162012_mas5_allGenes.txt'


# makeGeneTable(annotFile,probeAnnotationFile,probeExpressionFile,output)


# #=============================
# #MM1S MAS5 spike in normalized
# probeAnnotationFile  = '/nfs/young_ata4/gff/tables/platforms/GPL16043_probe_annotation.txt'

# annotFile = '/nfs/young_ata4/gff/annotation/hg18/hg18_refseq.ucsc'

# probeExpressionFile = '/nfs/young_ata4/pipeline/mm1s_new/expression/11162012_mas5_norm_allProbes.txt'

# output = '/nfs/young_ata4/pipeline/mm1s_new/expression/11162012_mas5_norm_allGenes.txt'


# makeGeneTable(annotFile,probeAnnotationFile,probeExpressionFile,output)


#=============================
#U87, H23, and H2171  MAS5 spike in normalized
#probeAnnotationFile  = '/nfs/young_ata4/gff/tables/platforms/GPL16043_probe_annotation.txt'

#annotFile = '/nfs/young_ata4/gff/annotation/hg18/hg18_refseq.ucsc'

#probeExpressionFile = '/nfs/young_ata4/pipeline/mm1s_new/expression/02072013_mas5_spikeNorm_allProbes.txt'

#output = '/nfs/young_ata4/pipeline/mm1s_new/expression/02072013_mas5_spikeNorm_allGenes.txt'


#makeGeneTable(annotFile,probeAnnotationFile,probeExpressionFile,output)
