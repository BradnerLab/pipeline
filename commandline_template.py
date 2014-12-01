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



#================================================================================
#============================GLOBAL PARAMETERS===================================
#================================================================================

#add locations of files and global parameters in this section


dataFile ='/location/file.txt'
genome = 'hg18'
projectFolder = '/grail/projects/gordon/'

#================================================================================
#===================================CLASSES======================================
#================================================================================

#user defined classes here

#================================================================================
#=================================FUNCTIONS======================================
#================================================================================

#write your specific functions here


def returnGenome(genome):

    '''
    prints the genome being used
    '''

    print "Using genome %s for analysis" % (genome)





#================================================================================
#===============================MAIN RUN=========================================
#================================================================================


#write the actual script here


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


    #optional arguments
    parser.add_option("-N","--mismatch",dest="mismatchN",nargs =1, default = 0,
                      help = "Specify 0 or 1 for allowed mismatches")
    parser.add_option("-r","--rpm", dest="rpm",action = 'store_true', default=True,
                      help = "Normalizes density to reads per million (rpm) Default is True")


    (options,args) = parser.parse_args()

    if not options.fastq:
        parser.print_help()
        exit()




if __name__ == "__main__":
    main()



