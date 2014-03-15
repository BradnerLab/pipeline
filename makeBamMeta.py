#!/usr/bin/python
#makeBamMeta.py

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

#makes a plus and minus strand meta from 3 gffs (tss,txn,ttr) and a bam


#parameters to set here
#Size of bins in the upstream region (bp)
tssBins = 50

#Number of bins in the txn region (nBins)
txnBins = 200

#Size of bins in the ttr region (bp)

ttrBins =50
import sys
print(sys.version)


sys.path.append('/ark/home/cl512/pipeline/')


from utils import *
import time

def main():

    from optparse import OptionParser
    usage ="usage: %prog [options] -g [COMMA_SEPARATED_GFFS_IN_ORDER] -b [SORTED_BAM_FILE] -o [OUTPUT_FOLDER]"
    parser = OptionParser(usage = usage)

    parser.add_option("-g", "--gff", dest = "gffList", nargs = 1,default = None,
                      help = "Comma separated list of gffs in order e.g. tssGFF,txnGFF,ttr,GFF")
    parser.add_option("-b", "--bam", dest = "bam", nargs = 1, default = None,
                      help = "Sorted bam file with a .bai in the same folder")
    parser.add_option("-o", "--out", dest = "output", nargs =1, default = None,
                      help = "Path of the output FOLDER")
    parser.add_option("-e", "--ext", dest = "ext", nargs =1, default = 200,
                      help = "Extension of reads")

    parser.add_option("-n", "--name", dest = "name", nargs = 1,default = None,
                      help = "Specify a name for the output files. default uses the bam name")
    parser.add_option("-f", "--finish", dest = "finish", action='store_true',default = False,
                      help = "If true, tries to finish the script and stitch gffs together")
    parser.add_option("-c", "--combine", dest = "combine", action = 'store_true',default = False,
                      help = "If true, combines the sense and antisense output into one single combined meta")
    (options,args) = parser.parse_args()

    print(options)
    print(args)

    if options.bam and options.gffList and options.output and not options.finish:

        jobIDRoot = options.name+ '_meta_'+join(str(time.time()).split('.'),'')
        bamFile = options.bam
        gffList = options.gffList
        outFolder = options.output
        extension = options.ext
        combine = options.combine

        gffList = gffList.split(',')
        if len(gffList) != 3:
            print('Must give 3 gffs in order TSS,TXN,TTR')
            exit()

        try:
            #check to make sure the output directory exists
            foo = os.listdir(outFolder)
        except OSError:
            print('OUTPUT DIRECTORY DOES NOT EXIST. CREATING IT NOW')
            os.system('mkdir %s' % (outFolder))

        [tssGFF,txnGFF,ttrGFF] = gffList
        bamName = join(bamFile.split('/')[-1].split('.')[0:2],'.')

        #next set get the name
        if options.name:
            name = options.name
        else:
            #use the txn gff as the root
            name = '%s_%s' % (bamName,txnGFF.split('/')[-1])
        print('using %s as the file name root' % (name))
        
        metaSettings = [['%s_JOB_ID' % (name)]]



        bashFileName = '%s/mapBamMeta_%s.sh' % (outFolder,jobIDRoot)
        bashFile = open(bashFileName,'w')
    
        sense = '+'
        tssOutfile  = '%s%s_%s_%s' % (outFolder,bamName,sense,tssGFF.split('/')[-1])
        txnOutfile  = '%s%s_%s_%s' % (outFolder,bamName,sense,txnGFF.split('/')[-1]) 
        ttrOutfile  = '%s%s_%s_%s' % (outFolder,bamName,sense,ttrGFF.split('/')[-1]) 
        for outFile in [tssOutfile,txnOutfile,ttrOutfile]:
            metaSettings.append([outFile.split('/')[-1]])



        job1ID = jobIDRoot + '_1'
        job2ID = jobIDRoot + '_2'
        job3ID = jobIDRoot + '_3'


        
        cmd1 = "python /ark/home/cl512/pipeline/bamToGFF_turbo.py -r -e %s -c %s -s %s -b %s -i %s  -o %s &" % (extension,tssBins,sense,bamFile,tssGFF,tssOutfile)
        cmd2 =  "python /ark/home/cl512/pipeline/bamToGFF_turbo.py -r -e %s -m %s -s %s -b %s -i %s  -o %s &" % (extension,txnBins,sense,bamFile,txnGFF,txnOutfile)
        cmd3 =  "python /ark/home/cl512/pipeline/bamToGFF_turbo.py -r -e %s -c %s -s %s -b %s -i %s  -o %s &" % (extension,ttrBins,sense,bamFile,ttrGFF,ttrOutfile)

        for cmd in [cmd1,cmd2,cmd3]:
            bashFile.write(cmd)
            bashFile.write('\n')        

        #launch jobs in the antisense strand

        sense = '-'
        tssOutfile  = '%s%s_%s_%s' % (outFolder,bamName,sense,tssGFF.split('/')[-1])
        txnOutfile  = '%s%s_%s_%s' % (outFolder,bamName,sense,txnGFF.split('/')[-1]) 
        ttrOutfile  = '%s%s_%s_%s' % (outFolder,bamName,sense,ttrGFF.split('/')[-1]) 
        for outFile in [tssOutfile,txnOutfile,ttrOutfile]:
            metaSettings.append([outFile.split('/')[-1]])

        job4ID = jobIDRoot + '_4'
        job5ID = jobIDRoot + '_5'
        job6ID = jobIDRoot + '_6'


        cmd4 = "python /ark/home/cl512/pipeline/bamToGFF_turbo.py -r -e %s -c %s -s %s -b %s -i %s  -o %s &" % (extension,tssBins,sense,bamFile,tssGFF,tssOutfile)
        cmd5 =  "python /ark/home/cl512/pipeline/bamToGFF_turbo.py -r -e %s -m %s -s %s -b %s -i %s  -o %s &" % (extension,txnBins,sense,bamFile,txnGFF,txnOutfile)
        cmd6 =  "python /ark/home/cl512/pipeline/bamToGFF_turbo.py -r -e %s -c %s -s %s -b %s -i %s  -o %s &" % (extension,ttrBins,sense,bamFile,ttrGFF,ttrOutfile)



        for cmd in [cmd4,cmd5,cmd6]:
            bashFile.write(cmd)
            bashFile.write('\n')

        
        unParseTable(metaSettings,outFolder+'%s_metaSettings.txt' % (name),'\t')

        #now launch the finishing jobs. waits until all 6 jobs are done.


        if combine:
            finishCommand = "python /ark/home/cl512/pipeline/makeBamMeta.py -c -f -n %s -o %s" % (name,options.output)
        else:
            finishCommand = "python /ark/home/cl512/pipeline/makeBamMeta.py -f -n %s -o %s" % (name,options.output)
        bashFile.write(finishCommand)
        bashFile.close()
        print bashFileName
        os.system('bash %s &' % (bashFileName))
    elif options.finish and options.name and options.output:
        #this is the finishing job
        #want to stitch two gffs together
        print('finishing and writing combined output')
        bamFile = options.bam
        gffList = options.gffList
        outFolder = options.output
        name = options.name
        combine = options.combine
        #make the sense output
        metaSettings = parseTable(outFolder+'%s_metaSettings.txt' % (name),'\t')
        
        [job1ID,job2ID,job3ID,job4ID,job5ID,job6ID]= [line[0] for line in metaSettings[1:]]

        if combine:

            tssOutSense = parseTable(outFolder + job1ID,'\t')
            txnOutSense = parseTable(outFolder + job2ID,'\t')
            ttrOutSense = parseTable(outFolder + job3ID,'\t')

            tssOutAnti = parseTable(outFolder + job4ID,'\t')
            txnOutAnti = parseTable(outFolder + job5ID,'\t')
            ttrOutAnti = parseTable(outFolder + job6ID,'\t')

            combinedOut = []
            header = tssOutSense[0] + txnOutSense[0][2:] + ttrOutSense[0][2:]
            combinedOut.append(header)

            for i in range(1,len(tssOutSense),1):
                headerLine = tssOutSense[i][0:2]
                try:
                    tssLine = [float(tssOutSense[i][j]) + float(tssOutAnti[i][j]) for j in range(2,len(tssOutSense[i]),1)]
                    txnLine = [float(txnOutSense[i][j]) + float(txnOutAnti[i][j]) for j in range(2,len(txnOutSense[i]),1)]
                    ttrLine = [float(ttrOutSense[i][j]) + float(ttrOutAnti[i][j]) for j in range(2,len(ttrOutSense[i]),1)]
                
                    newLine = headerLine + tssLine +txnLine + ttrLine
                except ValueError:
                    newLine = headerLine + ['NA']*320

                combinedOut.append(newLine)

            unParseTable(combinedOut,outFolder+'%s_meta.txt' % (name),'\t')

        else:

            
            tssOut = parseTable(outFolder + job1ID,'\t')
            txnOut = parseTable(outFolder + job2ID,'\t')
            ttrOut = parseTable(outFolder + job3ID,'\t')
            senseOut = []
            for i in range(len(tssOut)):

                newLine = tssOut[i] + txnOut[i][2:] + ttrOut[i][2:]
                senseOut.append(newLine)

            unParseTable(senseOut,outFolder+'%s_+_meta.txt' % (name),'\t')



            #make the antisensesense output

            tssOut = parseTable(outFolder + job4ID,'\t')
            txnOut = parseTable(outFolder + job5ID,'\t')
            ttrOut = parseTable(outFolder + job6ID,'\t')
            antiSenseOut = []
            for i in range(len(tssOut)):

                newLine = tssOut[i] + txnOut[i][2:] + ttrOut[i][2:]
                antiSenseOut.append(newLine)

            unParseTable(antiSenseOut,outFolder+'%s_-_meta.txt' % (name),'\t')
        
    else:
        parser.print_help()
        exit()

if __name__ == '__main__':
    main()
        
