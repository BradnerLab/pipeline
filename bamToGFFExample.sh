#!/usr/bin/bash
cd /home/cl512/pipeline/

# The MIT License (MIT)

# Copyright (c) 2013 Charles Lin

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

#bamToGFF_turbo.py quickly identifies read density in specified genomic regions and provides output in 
#normalized units of reads per million per basepair

#only works on sorted bam files
#to sort a bam file
#run these commands
#samtools -sort unsorted.bam sorted.


#HOW TO USE: 

#1. Make a copy of this script and rename it for each job.
#2. Edit variables in the parameters
#3. cd to folder containing script
#4. Use only code from 1 example.  Use # to comment out lines from a different example
#5. run by typing 'bash ./script.sh where script.sh is the name of your script




#Example 1 - determining density using fixed bin sizes

#====================================================
#====================PARAMETERS======================
#====================================================

#replace BAM1 with file path of a sorted bam file
BAM1='sample1.sorted.bam'



#edit these variables to specify input regions etc...
GFF='regions.gff' #visit https://genome.ucsc.edu/FAQ/FAQformat.html#format3 for a description of gff format
EXTENSION='200' #length each read is extended 
OUTPUT='OUTPUT_FILE' 
SENSE='both' #sense of reads plotted. use either +,-,both
BINSIZE=50 #size in bp of bins.  E.g. a 200bp gff region will be broken up into 4x50bp bins


echo
echo Calling bamToGFF_turbo.py using $BAM1 on gff $GFF and directing output to $OUTPUT
echo
#COMMAND
echo Running the following command:
echo python bamToGFF_turbo.py -b $BAM1 -i $GFF -o $OUTPUT -s $SENSE -e $EXTENSION -c $BINSIZE -r
python bamToGFF_turbo.py -b $BAM1 -i $GFF -o $OUTPUT -s $SENSE -e $EXTENSION -c $BINSIZE -r








#Example 2 - determining density using variable binning
#In this version, regions are broken up into a fixed number of bins regardless of size

#====================================================
#====================PARAMETERS======================
#====================================================

#replace BAM1 with file path of a sorted bam file
BAM1='sample1.sorted.bam'



#edit these variables to specify input regions etc...
GFF='regions.gff' #visit https://genome.ucsc.edu/FAQ/FAQformat.html#format3 for a description of gff format
EXTENSION='200' #length each read is extended 
OUTPUT='OUTPUT_FILE' 
SENSE='both' #sense of reads plotted. use either +,-,both
NBINS=10 #breaks each region into N bins. E.g. a 2000bp region with NBINS=10 will be broken up into 10 200bp bins


echo
echo Calling bamToGFF_turbo.py using $BAM1 on gff $GFF and directing output to $OUTPUT
echo
#COMMAND
echo Running the following command:
echo python bamToGFF_turbo.py -b $BAM1 -i $GFF -o $OUTPUT -s $SENSE -e $EXTENSION -m $NBINS -r
python bamToGFF_turbo.py -b $BAM1 -i $GFF -o $OUTPUT -s $SENSE -e $EXTENSION -m $NBINS -r


