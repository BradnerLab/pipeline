#!/usr/bin/env python

from os import sys, path

import bamliquidator_batch as blb
import normalize_plot_and_summarize as nps

import os
import shutil
import subprocess
import tables
import tempfile
import unittest

def create_single_full_read_bam(dir_path, chromosome, sequence):
    # create a sam file, based on instructions at http://genome.ucsc.edu/goldenPath/help/bam.html
    # and http://samtools.github.io/hts-specs/SAMv1.pdf
    sam_file_path = os.path.join(dir_path, 'single.sam') 
    with open(sam_file_path, 'w') as sam_file:
        length = len(sequence)
        sequence_header = '@SQ\tSN:%s\tLN:%d\n' % (chromosome, length)
        qual = '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
        sam_file.write(sequence_header)
        #       qname      chr    quality   next read name                                                   
        #       |      flag|   pos|    CIGAR|  next read pos
        #       |      |   |   |  |    |    |  |  template length 
        #       |      |   |   |  |    |    |  |  |  sequence
        #       |      |   |   |  |    |    |  |  |  |   QUAL           
        #       |      |   |   |  |    |    |  |  |  |   |   distance to ref
        sam_file.write('read1\t16\t%s\t1\t255\t50M\t*\t0\t0\t%s\t%s\tNM:i:0\n' % (chromosome, sequence, qual))
   
    # create bam file
    bam_file_path = os.path.join(dir_path, 'single.bam')
    subprocess.check_call(['samtools', 'view', '-S', '-b', '-o', bam_file_path, sam_file_path])

    # skipping sorting since it is already sorted

    subprocess.check_call(['samtools', 'index', bam_file_path])

    return bam_file_path

def create_single_region_file(dir_path, chromosome, start, stop, strand='.'):
    region_file_path = os.path.join(dir_path, 'single.gff') 
    with open(region_file_path, 'w') as region_file:
        region_file.write('%s\tregion1\t\t%d\t%d\t\t%s\t\tregion1\n' % (chromosome, start, stop, strand))
    return region_file_path


class SingleFullReadBamTest(unittest.TestCase):
    def setUp(self):
        self.dir_path = tempfile.mkdtemp()
        self.chromosome = 'chr1'
        self.sequence = 'ATTTAAAAATTAATTTAATGCTTGGCTAAATCTTAATTACATATATAATT'
        self.bam_file_path = create_single_full_read_bam(self.dir_path, self.chromosome, self.sequence)

    def tearDown(self):
        shutil.rmtree(self.dir_path)

    def test_bin_liquidation(self):
        liquidator = blb.BinLiquidator(bin_size = len(self.sequence),
                                       output_directory = os.path.join(self.dir_path, 'output'),
                                       bam_file_path = self.bam_file_path)
        liquidator.batch(extension = 0, sense = '.')

        with tables.open_file(liquidator.counts_file_path) as counts:
            self.assertEqual(1, len(counts.root.lengths)) # 1 since only a single bam file
            self.assertEqual(1, counts.root.lengths[0]["length"]) # 1 since only a single read 

            self.assertEqual(1, len(counts.root.bin_counts)) # 1 since 1 bin accommodates full sequence 
           
            record = counts.root.bin_counts[0]
            self.assertEqual(0, record["bin_number"])
            self.assertEqual(self.chromosome, record["chromosome"])
            self.assertEqual(len(self.sequence), record["count"]) # count represents how many base pair reads intersected
                                                                  # the bin

    def test_region_liquidation(self):
        start = 1
        stop  = 8
        regions_file_path = create_single_region_file(self.dir_path, self.chromosome, start, stop)

        liquidator = blb.RegionLiquidator(regions_file = regions_file_path,
                                          output_directory = os.path.join(self.dir_path, 'output'),
                                          bam_file_path = self.bam_file_path)
        liquidator.batch(extension = 0, sense = '.')

        with tables.open_file(liquidator.counts_file_path) as counts:
            self.assertEqual(1, len(counts.root.lengths)) # 1 since only a single bam file
            self.assertEqual(1, counts.root.lengths[0]["length"]) # 1 since only a single read 
            self.assertEqual(1, len(counts.root.region_counts)) # 1 since only a single region 
           
            record = counts.root.region_counts[0]
            self.assertEqual(start, record["start"])
            self.assertEqual(stop,  record["stop"])
            self.assertEqual(stop-start, record["count"]) # count represents how many base pair reads intersected
                                                          # the region

    def test_region_with_no_reads(self):
        start = len(self.sequence) + 10
        stop = start + 10
        regions_file_path = create_single_region_file(self.dir_path, self.chromosome, start, stop)

        liquidator = blb.RegionLiquidator(regions_file = regions_file_path,
                                          output_directory = os.path.join(self.dir_path, 'output'),
                                          bam_file_path = self.bam_file_path)
        liquidator.batch(extension = 0, sense = '.')

        with tables.open_file(liquidator.counts_file_path) as counts:
            self.assertEqual(1, len(counts.root.lengths)) # 1 since only a single bam file
            self.assertEqual(1, counts.root.lengths[0]["length"]) # 1 since only a single read 

            record = counts.root.region_counts[0]
            self.assertEqual(start, record["start"])
            self.assertEqual(stop,  record["stop"])
            self.assertEqual(0, record["count"]) # 0 since region doesn't intersect sequence

    def test_region_with_wrong_chromosome(self):
        start = len(self.sequence) + 10
        stop = start + 10
        regions_file_path = create_single_region_file(self.dir_path, self.chromosome + '0', start, stop)

        liquidator = blb.RegionLiquidator(regions_file = regions_file_path,
                                          output_directory = os.path.join(self.dir_path, 'output'),
                                          bam_file_path = self.bam_file_path)
        liquidator.batch(extension = 0, sense = '.')

        with tables.open_file(liquidator.counts_file_path) as counts:
            self.assertEqual(1, len(counts.root.lengths)) # 1 since only a single bam file
            self.assertEqual(1, counts.root.lengths[0]["length"]) # 1 since only a single read 

            record = counts.root.region_counts[0]
            self.assertEqual(start, record["start"])
            self.assertEqual(stop,  record["stop"])
            self.assertEqual(0, record["count"]) # 0 since region doesn't intersect sequence

if __name__ == '__main__':
    unittest.main()

##################################################################################
# The MIT License (MIT)
#
# Copyright (c) 2013 John DiMatteo (jdimatteo@gmail.com)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the "Software"), to deal in
# the Software without restriction, including without limitation the rights to
# use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
# the Software, and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
# FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
# COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
# IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
##################################################################################

