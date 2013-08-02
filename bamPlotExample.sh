#!/usr/bin/bash
cd /home/cl512/pipeline/

#bamPlot.py makes publication quality tracks of ChIP-Seq data that can be edited in adobe illustrator


#HOW TO USE: 

#1. Make a copy of this script and rename it for each job.
#2. Edit variables in the parameters
#3. cd to folder containing script
#4. run by typing 'bash ./script.sh where script.sh is the name of your script


#====================================================
#====================PARAMETERS======================
#====================================================

#replace BAM1,BAM2 with file paths of sorted bam files
#can add BAM3 and NAME3  etc... using same convention

BAM1='sample1.sorted.bam'
BAM2='sample2.sorted.bam'

#Edit the names to give each bam a title
NAME1='sample1'
NAME2='sample2'

#edit these variables to specify region, genome, title etc...
REGION='chr1:+:1001000-1002000'  
GENOME='HG18'
OUTPUT='OUTPUTFOLDER'
TITLE='TITLE'
YAXIS='UNIFORM'    #use either UNIFORM or RELATIVE
SENSE='both' #sense of reads plotted. use either 
COLOR='255,0,0:0,0,255' #use a colon separated list of RGB colors w/ values from 0 to 255


echo
echo Using $BAM1 as BAM1
echo Using $BAM1 as BAM1
echo Calling bamPlot.py on region $REGION in genome $GENOME and directing output to $TITLE in $OUTPUT
echo
echo Running the following command:
#COMMAND
echo python bamPlot.py -b $BAM1,$BAM2 -i $REGION -t $TITLE -y UNIFORM -o $OUTPUT -n $NAME1,$NAME2 -p single -c $COLOR -s $SENSE
python bamPlot.py -b $BAM1,$BAM2 -i $REGION -t $TITLE -y UNIFORM -o $OUTPUT -n $NAME1,$NAME2 -p single -c $COLOR -s $SENSE
