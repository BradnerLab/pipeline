#!/bin/bash
#
# BAM TO GFF to detect Read density across gene transcripting regions
#
# Edited from Original version 12/11/2019

##############################################################
# ##### Please replace PATHTO with your own directory ###### #
##############################################################
# NOTE: We now assume the PATH and PYTHONPATH has been set up before 
# entering this script
#PATHTO=/PATH/TO/BAM2GFF
#PYTHONPATH=$PATHTO/lib
#export PYTHONPATH
#export PATH=$PATH:$PATHTO/src

if [ $# -lt 3 ]; then
  echo ""
  echo 1>&2 Usage: $0 ["GTF file"] ["feature type"] ["BAM file"] ["CHROM SIZES"] ["SAMPLENAME"]
  echo ""
  exit 1
fi

#================================================================================
#Parameters for running

# GTF files
GTFFILE=$1

#FEATURE TYPE
FEATURE=$2
FEATURE=${FEATURE:=gene}

# BAM file
BAMFILE=$3

# CHROM SIZES
CHROMSIZES=$4

#sample name
SAMPLENAME=$5
SAMPLENAME=${SAMPLENAME:=BAMFILE##*/}

echo "#############################################"
echo "######            BAM2GFF v1           ######"
echo "#############################################"

echo "BAM file: $BAMFILE"
echo "FEATURE type: $FEATURE"
echo "Sample Name: $SAMPLENAME"
#================================================================================
#
# GENERATING GFF files for each genomic region
#
mkdir -p annotation
echo "BAM2GFF_gtftogenes.py -g $GTFFILE -f $FEATURE -c $CHROMSIZES"
BAM2GFF_gtftogenes.py -g $GTFFILE -f $FEATURE -c $CHROMSIZES
echo

#
# BAM TO GFF main code
#
mkdir -p matrix
echo "Working on GeneBody Region\nBAM2GFF_main.py -b $BAMFILE -i annotation/genes.gff -m 100 -o matrix/genebody.txt"
BAM2GFF_main.py -b $BAMFILE -i annotation/genes.gff -m 100 -o matrix/genebody.txt
echo
echo "Working on Upstream Region\nBAM2GFF_main.py -b $BAMFILE -i annotation/upstream.gff -m 50 -o matrix/upstream.txt"
BAM2GFF_main.py -b $BAMFILE -i annotation/upstream.gff -m 50 -o matrix/upstream.txt
echo
echo "Working on Downstream Region\nBAM2GFF_main.py -b $BAMFILE -i annotation/downstream.gff -m 50 -o matrix/downstream.txt"
BAM2GFF_main.py -b $BAMFILE -i annotation/downstream.gff -m 50 -o matrix/downstream.txt
echo
echo "Working on Promoter Region\nBAM2GFF_main.py -b $BAMFILE -i annotation/promoters.gff -m 100 -o matrix/promoters.txt"
BAM2GFF_main.py -b $BAMFILE -i annotation/promoters.gff -m 100 -o matrix/promoters.txt
echo

#
# PLOTS
#
echo "BAM2GFF_plots.R $SAMPLENAME"
BAM2GFF_plots.R $SAMPLENAME

echo "Done!"
