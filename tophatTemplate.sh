#!/usr/bin/bash

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


PROJECT_DIR='/ark/home/cl512/projects/ezhi/'
NAME='K422_cmpd5_rep1'
FASTQ='/ark/home/cl512/ressrv19/raw/130410Bra/D13-1625/130410Bra_D13-1625_2_sequence.fastq'
GENOME='hg18'

mkdir $PROJECT_DIR$NAME
tophat -p 4 -o $PROJECT_DIR$NAME/ --transcriptome-index=/ark/home/cl512/ressrv19/genomes/transcriptome_data/hg18_genes /ark/home/cl512/ressrv19/genomes/hg18_withERCC/hg18_ercc_noRand $FASTQ

samtools sort $PROJECT_DIR$NAME/accepted_hits.bam $PROJECT_DIR$NAME/$NAME.$GENOME.tophat.sorted 
samtools index $PROJECT_DIR$NAME/$NAME.$GENOME.tophat.sorted.bam



python /usr/local/bin/RPKM_count.py -r '/ark/home/cl512/ressrv19/genomes/ERCC_Technical_Data/ERCC92.bed' -i $PROJECT_DIR$NAME/$NAME.$GENOME.tophat.sorted.bam -e -o $PROJECT_DIR$NAME/$NAME\_ERCC

cufflinks -p 4 -G /ark/home/cl512/ressrv19/annotations/hg18_genes.gtf -o $PROJECT_DIR$NAME $PROJECT_DIR$NAME/$NAME.$GENOME.tophat.sorted.bam