#!/usr/bin/python

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
sys.path.append('/home/cl512/src/pipeline/')
import utils



#================================================================================
#============================GLOBAL PARAMETERS===================================
#================================================================================

#add locations of files and global parameters in this section


dataFile ='/location/file.txt'
genome = 'hg18'


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
    this is the main run function for the script
    all of the work should occur here, but no functions should be defined here
    '''


    returnGenome(genome)






main()

