#!/usr/bin/env python

from __future__ import division

from os import sys, path
from itertools import izip

import bamliquidator_batch as blb
import normalize_plot_and_summarize as nps

import os
import shutil
import subprocess
import tables
import tempfile
import unittest

# one full read for each chromosome
def create_bam(dir_path, chromosomes, sequence, file_name='single.bam', flag=16):
    # create a sam file, based on instructions at http://genome.ucsc.edu/goldenPath/help/bam.html
    # and http://samtools.github.io/hts-specs/SAMv1.pdf
    sam_file_path = os.path.join(dir_path, 'single.sam') 
    with open(sam_file_path, 'w') as sam_file:
        length = len(sequence)
        
        # sequence headers for all chromosomes go on top
        for chromosome in chromosomes: 
            sequence_header = '@SQ\tSN:%s\tLN:%d\n' % (chromosome, length)
            sam_file.write(sequence_header)

        for chromosome in chromosomes:
            qual = len(sequence)*'<'
            #               qname      chr    quality   next read name                                                   
            #               |      flag|   pos|    CIGAR|  next read pos
            #               |      |   |   |  |    |    |  |  template length 
            #               |      |   |   |  |    |    |  |  |  sequence
            #               |      |   |   |  |    |    |  |  |  |   QUAL           
            #               |      |   |   |  |    |    |  |  |  |   |   distance to ref
            sam_file.write('read1\t%d\t%s\t1\t255\t%dM\t*\t0\t0\t%s\t%s\tNM:i:0\n' % (flag, chromosome, len(sequence), sequence, qual))
   
    # create bam file
    bam_file_path = os.path.join(dir_path, file_name)
    subprocess.check_call(['samtools', 'view', '-S', '-b', '-o', bam_file_path, sam_file_path])

    # skipping sorting since it is already sorted

    subprocess.check_call(['samtools', 'index', bam_file_path])

    return bam_file_path

def create_single_region_gff_file(dir_path, chromosome, start, stop, strand='.', file_name = 'single.gff', region_name='region1'):
    region_file_path = os.path.join(dir_path, file_name) 
    with open(region_file_path, 'w') as region_file:
        region_file.write('%s\t%s\t\t%d\t%d\t\t%s\t\tregion1\n' % (chromosome, region_name, start, stop, strand))
    return region_file_path

def create_single_region_bed_file(dir_path, chromosome, start, stop, region_name='region1'):
    region_file_path = os.path.join(dir_path, 'single.bed') 
    with open(region_file_path, 'w') as region_file:
        region_file.write('%s\t%d\t%d\t%s\t3156.59\n' % (chromosome, start, stop, region_name))
    return region_file_path

class TempDirTest(unittest.TestCase):
    def setUp(self):
        self.dir_path = tempfile.mkdtemp(prefix='blt_')

    def tearDown(self):
        #print 'tearing down', self.dir_path
        shutil.rmtree(self.dir_path)

class SingleFullReadBamTest(TempDirTest):
    def setUp(self):
        super(SingleFullReadBamTest, self).setUp()
        self.chromosome = 'chr1'
        self.sequence = 'ATTTAAAAATTAATTTAATGCTTGGCTAAATCTTAATTACATATATAATT'
        self.bam_file_path = create_bam(self.dir_path, [self.chromosome], self.sequence)

    def test_bin_liquidation(self):
        bin_size = len(self.sequence)
        liquidator = blb.BinLiquidator(bin_size = bin_size,
                                       output_directory = os.path.join(self.dir_path, 'output'),
                                       bam_file_path = self.bam_file_path)

        liquidator.flatten()

        with tables.open_file(liquidator.counts_file_path) as counts:
            self.assertEqual(1, len(counts.root.files)) # 1 since only a single bam file
            file_record = counts.root.files[0] 
            self.assertEqual(1, file_record['length']) # 1 since only a single read 
            self.assertEqual(1, file_record['key'])

            self.assertEqual(1, len(counts.root.bin_counts)) # 1 since 1 bin accommodates full sequence 
           
            record = counts.root.bin_counts[0]
            self.assertEqual(0, record['bin_number'])
            self.assertEqual(self.chromosome, record['chromosome'])
            self.assertEqual(len(self.sequence), record['count']) # count represents how many base pair reads 
                                                                  # intersected the bin

    def test_bin_liquidation_zero_bin_size(self):
        with self.assertRaises(Exception):
            liquidator = blb.BinLiquidator(bin_size = 0,
                                           output_directory = os.path.join(self.dir_path, 'output'),
                                           bam_file_path = self.bam_file_path)
            liquidator.batch(extension = 0, sense = '.')

    def test_region_liquidation(self):
        start = 1
        stop  = 8
        region_name = 'region_f'
        regions_file_path = create_single_region_gff_file(self.dir_path, self.chromosome, start, stop, region_name=region_name)

        liquidator = blb.RegionLiquidator(regions_file = regions_file_path,
                                          output_directory = os.path.join(self.dir_path, 'output'),
                                          bam_file_path = self.bam_file_path)
        liquidator.flatten()

        matrix_path = os.path.join(self.dir_path, 'matrix.gff')
        blb.write_bamToGff_matrix(matrix_path, liquidator.counts_file_path)

        with tables.open_file(liquidator.counts_file_path) as counts:
            self.assertEqual(1, len(counts.root.files)) # 1 since only a single bam file
            self.assertEqual(1, counts.root.files[0]['length']) # 1 since only a single read 
            self.assertEqual(1, len(counts.root.region_counts)) # 1 since only a single region 
           
            record = counts.root.region_counts[0]
            self.assertEqual(region_name, record['region_name'])
            self.assertEqual(start, record['start'])
            self.assertEqual(stop,  record['stop'])
            self.assertEqual(stop-start, record['count']) # count represents how many base pair reads intersected
                                                          # the region

            # todo: add normalization record checks

        with open(matrix_path, 'r') as matrix_file:
           matrix_lines = matrix_file.readlines()
           self.assertEqual(2, len(matrix_lines))

           header_cols = matrix_lines[0].split('\t')
           self.assertEqual(3, len(header_cols))
           self.assertEqual('GENE_ID', header_cols[0])
           self.assertEqual('locusLine', header_cols[1])
           self.assertEqual('bin_1_%s\n' % 'single.bam', header_cols[2])

           data_cols = matrix_lines[1].split('\t')
           self.assertEqual(3, len(data_cols))
           self.assertEqual(region_name, data_cols[0])
           self.assertEqual('chr1(.):1-8', data_cols[1]) # todo: don't hardcode these values 
           self.assertEqual('1000000.0\n', data_cols[2])

    def test_region_liquidation_with_optional_bed_columns(self):
        start = 1
        stop  = 8
        region_name = 'region_f'
        score = 9.9
        strand = "."

        columns = ["chr1", str(start), str(stop), region_name, str(score), strand,
                   "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts"]

        def create_single_region_bed_file(dir_path, columns):
            region_file_path = os.path.join(dir_path, "single.bed") 
            with open(region_file_path, 'w') as region_file:
                first_column_written = False
                for column in columns:
                    if first_column_written:
                        region_file.write('\t')
                    region_file.write('%s' % column)
                    first_column_written = True
                region_file.write('\n')
            return region_file_path
            
        for number_of_columns in range(0, 9):
            regions_file_path = create_single_region_bed_file(self.dir_path, columns[:number_of_columns])

            try:
                liquidator = blb.RegionLiquidator(regions_file = regions_file_path,
                                                  output_directory = os.path.join(self.dir_path, 'output'),
                                                  bam_file_path = self.bam_file_path)
            except Exception:
                self.assertLess(number_of_columns, 3)
                continue
            else:
                self.assertGreaterEqual(number_of_columns, 3)

            liquidator.flatten()

            matrix_path = os.path.join(self.dir_path, 'matrix.gff')
            blb.write_bamToGff_matrix(matrix_path, liquidator.counts_file_path)

            with tables.open_file(liquidator.counts_file_path) as counts:
                self.assertEqual(1, len(counts.root.files)) # 1 since only a single bam file
                self.assertEqual(1, counts.root.files[0]['length']) # 1 since only a single read 
                self.assertEqual(1, len(counts.root.region_counts)) # 1 since only a single region 
               
                record = counts.root.region_counts[0]
                if number_of_columns >= 4:
                   self.assertEqual(region_name, record['region_name'])
                else:
                   self.assertEqual("", record['region_name'])
                self.assertEqual(start, record['start'])
                self.assertEqual(stop,  record['stop'])
                self.assertEqual(stop-start, record['count']) # count represents how many base pair reads intersected
                                                              # the region

            with open(matrix_path, 'r') as matrix_file:
               matrix_lines = matrix_file.readlines()
               self.assertEqual(2, len(matrix_lines))

               header_cols = matrix_lines[0].split('\t')
               self.assertEqual(3, len(header_cols))
               self.assertEqual('GENE_ID', header_cols[0])
               self.assertEqual('locusLine', header_cols[1])
               self.assertEqual('bin_1_%s\n' % 'single.bam', header_cols[2])

               data_cols = matrix_lines[1].split('\t')
               self.assertEqual(3, len(data_cols))
               if number_of_columns >= 4:
                   self.assertEqual(region_name, data_cols[0])
               else:
                   self.assertEqual("", data_cols[0])
               self.assertEqual('chr1(.):1-8', data_cols[1]) # todo: don't hardcode these values 
               self.assertEqual('1000000.0\n', data_cols[2])

    def test_out_of_range_region(self):
        start = len(self.sequence) + 10
        stop = start + 10
        regions_file_path = create_single_region_gff_file(self.dir_path, self.chromosome, start, stop)

        liquidator = blb.RegionLiquidator(regions_file = regions_file_path,
                                          output_directory = os.path.join(self.dir_path, 'output'),
                                          bam_file_path = self.bam_file_path)

        with tables.open_file(liquidator.counts_file_path) as counts:
            self.assertEqual(1, len(counts.root.files)) # 1 since only a single bam file
            self.assertEqual(1, counts.root.files[0]['length']) # 1 since only a single read 
            
            self.assertEqual(0, len(counts.root.region_counts))

    def test_empty_region_file(self):
        empty_file_path = os.path.join(self.dir_path, 'empty.gff')
        open(empty_file_path, 'w').close()

        liquidator = blb.RegionLiquidator(regions_file = empty_file_path,
                                          output_directory = os.path.join(self.dir_path, 'output'),
                                          bam_file_path = self.bam_file_path)
        liquidator.batch(extension = 0, sense = '.')
        

    def test_region_with_wrong_chromosome(self):
        start = len(self.sequence) + 10
        stop = start + 10
        regions_file_path = create_single_region_gff_file(self.dir_path, self.chromosome + '0', start, stop)

        liquidator = blb.RegionLiquidator(regions_file = regions_file_path,
                                          output_directory = os.path.join(self.dir_path, 'output'),
                                          bam_file_path = self.bam_file_path)

        with tables.open_file(liquidator.counts_file_path) as counts:
            self.assertEqual(1, len(counts.root.files)) # 1 since only a single bam file
            self.assertEqual(1, counts.root.files[0]['length']) # 1 since only a single read 
            self.assertEqual(0, len(counts.root.region_counts)) # no valid regions

    def test_region_with_long_name(self):
        start = 1
        stop  = 8
        region_name = 'r'*64
        truncated_region_name = 'r'*63
        regions_file_path = create_single_region_gff_file(self.dir_path, self.chromosome, start, stop, region_name=region_name)

        liquidator = blb.RegionLiquidator(regions_file = regions_file_path,
                                          output_directory = os.path.join(self.dir_path, 'output'),
                                          bam_file_path = self.bam_file_path)
        liquidator.flatten()

        matrix_path = os.path.join(self.dir_path, 'matrix.gff')
        blb.write_bamToGff_matrix(matrix_path, liquidator.counts_file_path)

        with tables.open_file(liquidator.counts_file_path) as counts:
            self.assertEqual(1, len(counts.root.files)) # 1 since only a single bam file
            self.assertEqual(1, counts.root.files[0]['length']) # 1 since only a single read 
            self.assertEqual(1, len(counts.root.region_counts)) # 1 since only a single region 
           
            record = counts.root.region_counts[0]
            self.assertEqual(truncated_region_name, record['region_name'])
            self.assertEqual(start, record['start'])
            self.assertEqual(stop,  record['stop'])
            self.assertEqual(stop-start, record['count']) # count represents how many base pair reads intersected
                                                          # the region
    
    def test_region_with_really_long_name(self):
        start = 1
        stop  = 8
        region_name = 'r'*84
        truncated_region_name = 'r'*63
        regions_file_path = create_single_region_gff_file(self.dir_path, self.chromosome, start, stop, region_name=region_name)

        liquidator = blb.RegionLiquidator(regions_file = regions_file_path,
                                          output_directory = os.path.join(self.dir_path, 'output'),
                                          bam_file_path = self.bam_file_path)
        liquidator.flatten()

        matrix_path = os.path.join(self.dir_path, 'matrix.gff')
        blb.write_bamToGff_matrix(matrix_path, liquidator.counts_file_path)

        with tables.open_file(liquidator.counts_file_path) as counts:
            self.assertEqual(1, len(counts.root.files)) # 1 since only a single bam file
            self.assertEqual(1, counts.root.files[0]['length']) # 1 since only a single read 
            self.assertEqual(1, len(counts.root.region_counts)) # 1 since only a single region 
           
            record = counts.root.region_counts[0]
            self.assertEqual(truncated_region_name, record['region_name'])
            self.assertEqual(start, record['start'])
            self.assertEqual(stop,  record['stop'])
            self.assertEqual(stop-start, record['count']) # count represents how many base pair reads intersected
                                                          # the region

    def helper_check_region_with_chromosome(self, chromosome):
        start = 1 
        stop = 8 
        regions_file_path = create_single_region_gff_file(self.dir_path, chromosome, start, stop)
        bam_file_path = create_bam(self.dir_path, [chromosome], self.sequence)

        liquidator = blb.RegionLiquidator(regions_file = regions_file_path,
                                          output_directory = os.path.join(self.dir_path, 'output'),
                                          bam_file_path = bam_file_path)

        with tables.open_file(liquidator.counts_file_path) as counts:
            self.assertEqual(1, len(counts.root.files)) # 1 since only a single bam file
            self.assertEqual(1, counts.root.files[0]['length']) # 1 since only a single read 

            record = counts.root.region_counts[0]
            self.assertEqual(start, record['start'])
            self.assertEqual(stop,  record['stop'])
            self.assertEqual(stop-start, record['count']) # count represents how many base pair reads intersected
                                                          # the region
            counts.root.files[0]['length']
            factor = (1 / (stop-start)) * (1 / (1 / 10**6))
            self.assertEqual(record['count'] * factor, record['normalized_count'])

    def test_region_with_long_chromosome(self):
        self.helper_check_region_with_chromosome('a'*63)

    def test_region_with_too_long_chromosome(self):
        with self.assertRaises(Exception):
            self.helper_check_region_with_chromosome('a'*64)

    def test_region_with_non_canonical_chromosome(self):
        self.helper_check_region_with_chromosome('chr7_blah_a')

    def test_region_with_black_listed_chromosome_pattern(self):
        self.helper_check_region_with_chromosome('chr7_random')

    def test_bin_long_bam_file_name(self):
        long_file_name = 'x' * 65 # more than Float64Col 
        long_file_path = os.path.join(self.dir_path, long_file_name)
        shutil.copyfile(self.bam_file_path, long_file_path) 
        shutil.copyfile(self.bam_file_path + '.bai', long_file_path + '.bai') 
        bin_liquidator = blb.BinLiquidator(bin_size = len(self.sequence),
                                           output_directory = os.path.join(self.dir_path, 'bin_output'),
                                           bam_file_path = long_file_path)

    def test_region_long_bam_file_name(self):
        long_file_name = 'x' * 65 # more than Float64Col 
        long_file_path = os.path.join(self.dir_path, long_file_name)
        shutil.copyfile(self.bam_file_path, long_file_path) 
        shutil.copyfile(self.bam_file_path + '.bai', long_file_path + '.bai') 

        start = 1
        stop  = len(self.sequence) 
        regions_file_path = create_single_region_gff_file(self.dir_path, self.chromosome, start, stop)
        region_liquidator = blb.RegionLiquidator(regions_file = regions_file_path,
                                                 output_directory = os.path.join(self.dir_path, 'region_output'),
                                                 bam_file_path = long_file_path) 

    def test_region_other_extension(self):
        start = 1
        stop  = len(self.sequence) 
        regions_file_path = create_single_region_gff_file(self.dir_path, self.chromosome, start, stop,
                                                          file_name = 'single.txt')
        region_liquidator = blb.RegionLiquidator(regions_file = regions_file_path,
                                                 region_format = 'gff',
                                                 output_directory = os.path.join(self.dir_path, 'region_output'),
                                                 bam_file_path = self.bam_file_path)
    
class MultipleChromosomeTest(TempDirTest):
    def testBinLiquidation(self):
        chromosomes = ['chr1', 'chr2']
        sequence = 'ATTTAAAAATTAATTTAATGCTTGGCTAAATCTTAATTACATATATAATT'
        bam_file_path = create_bam(self.dir_path, chromosomes, sequence, file_name='multiple.bam')
        bin_size = len(sequence)

        liquidator = blb.BinLiquidator(bin_size = bin_size,
                                       output_directory = os.path.join(self.dir_path, 'output'),
                                       bam_file_path = bam_file_path)

        with tables.open_file(liquidator.counts_file_path) as counts:
            self.assertEqual(1, len(counts.root.files)) # 1 since only a single bam file
            file_record = counts.root.files[0] 
            self.assertEqual(2, file_record['length']) # 1 read for each chromosome
            self.assertEqual(1, file_record['key'])

            self.assertEqual(2, len(counts.root.bin_counts)) # 1 for each chromosome since bin accommodates full sequence
           
            for record_index, chromosome in enumerate(chromosomes):
                record = counts.root.bin_counts[record_index]
                self.assertEqual(0, record['bin_number'])
                self.assertEqual(chromosome, record['chromosome'])
                self.assertEqual(len(sequence), record['count']) # count represents how many base pair reads 
                                                                      # intersected the bin

    def testDefaultBlackList(self):
        chromosomes_that_should_be_liquidated = ['chr1', 'chr2']
        all_chromosomes = chromosomes_that_should_be_liquidated + ['chr2_random']

        sequence = 'ATTTAAAAATTAATTTAATGCTTGGCTAAATCTTAATTACATATATAATT'
        bam_file_path = create_bam(self.dir_path, all_chromosomes, sequence, file_name='multiple.bam')
        bin_size = len(sequence)

        liquidator = blb.BinLiquidator(bin_size = bin_size,
                                       output_directory = os.path.join(self.dir_path, 'output'),
                                       bam_file_path = bam_file_path)

        with tables.open_file(liquidator.counts_file_path) as counts:
            self.assertEqual(1, len(counts.root.files)) # 1 since only a single bam file
            file_record = counts.root.files[0] 
            self.assertEqual(3, file_record['length']) # 1 read for each chromosome
            self.assertEqual(1, file_record['key'])

            self.assertEqual(2, len(counts.root.bin_counts)) # 1 for each chromosome that should be liquidated
           
            for record_index, chromosome in enumerate(chromosomes_that_should_be_liquidated):
                record = counts.root.bin_counts[record_index]
                self.assertEqual(0, record['bin_number'])
                self.assertEqual(chromosome, record['chromosome'])
                self.assertEqual(len(sequence), record['count']) # count represents how many base pair reads 

    def testOverridingBlackList(self):
        chromosomes = ['chr1', 'chr2', 'chr2_random']

        sequence = 'ATTTAAAAATTAATTTAATGCTTGGCTAAATCTTAATTACATATATAATT'
        bam_file_path = create_bam(self.dir_path, chromosomes, sequence, file_name='multiple.bam')
        bin_size = len(sequence)

        liquidator = blb.BinLiquidator(bin_size = bin_size,
                                       output_directory = os.path.join(self.dir_path, 'output'),
                                       bam_file_path = bam_file_path,
                                       blacklist = [])

        with tables.open_file(liquidator.counts_file_path) as counts:
            self.assertEqual(1, len(counts.root.files)) # 1 since only a single bam file
            file_record = counts.root.files[0] 
            self.assertEqual(3, file_record['length']) # 1 read for each chromosome
            self.assertEqual(1, file_record['key'])

            self.assertEqual(3, len(counts.root.bin_counts)) # 1 for each chromosome that should be liquidated
           
            for record_index, chromosome in enumerate(chromosomes):
                record = counts.root.bin_counts[record_index]
                self.assertEqual(0, record['bin_number'])
                self.assertEqual(chromosome, record['chromosome'])
                self.assertEqual(len(sequence), record['count']) # count represents how many base pair reads 

class AppendingTest(TempDirTest):
    def setUp(self):
        super(AppendingTest, self).setUp()
        self.chromosome = 'chr1'
        self.sequence1 = 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'
        self.bam1_file_path = create_bam(self.dir_path, [self.chromosome], self.sequence1, 'single1.bam')
        self.sequence2 = 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'
        self.bam2_file_path = create_bam(self.dir_path, [self.chromosome], self.sequence2, 'single2.bam')

    def testBin(self):
        # first do all liquidation together at once, then do it in two runs appending, and verify everything matches
        bin_size = len(self.sequence1)
        together_dir_path = os.path.join(self.dir_path, 'together')
        liquidator = blb.BinLiquidator(bin_size = bin_size,
                                       output_directory = together_dir_path, 
                                       bam_file_path = self.dir_path)

        appending_dir = os.path.join(self.dir_path, 'appending')
        #print 'liquidating bams at path:', self.bam1_file_path
        liquidator = blb.BinLiquidator(bin_size = bin_size,
                                       output_directory = appending_dir,
                                       bam_file_path = self.bam1_file_path)

        appending_h5_path = os.path.join(appending_dir, 'counts.h5')
        liquidator = blb.BinLiquidator(bin_size = bin_size,
                                       output_directory = os.path.join(self.dir_path, 'appending_extra_without_h5_file'),
                                       bam_file_path = self.bam2_file_path,
                                       counts_file_path = appending_h5_path) 

        with tables.open_file(os.path.join(together_dir_path, 'counts.h5')) as together_h5:
            with tables.open_file(appending_h5_path) as appending_h5:
                self.assertEqual(str(together_h5.root.bin_counts[:]), str(appending_h5.root.bin_counts[:]))
                self.assertEqual(str(together_h5.root.normalized_counts[:]), str(appending_h5.root.normalized_counts[:]))
                self.assertEqual(str(together_h5.root.summary[:]), str(appending_h5.root.summary[:]))
                self.assertEqual(str(together_h5.root.sorted_summary[:]), str(appending_h5.root.sorted_summary[:]))

class LiquidateBamInDifferentDirectories(unittest.TestCase):
    def setUp(self):
        self.dir_before = os.getcwd()
        self.chromosome = 'chr1'
        self.sequence = 'ATTTAAAAATTAATTTAATGCTTGGCTAAATCTTAATTACATATATAATT'
        self.dir_path = None

    def tearDown(self):
        #print 'tearing down', self.dir_path
        os.chdir(self.dir_before)
        shutil.rmtree(self.dir_path)

    def test_liquidation_in_current_directory(self):
        self.dir_path = tempfile.mkdtemp(prefix='blt_')
        bam_file_name = "current.bam"
        self.bam_file_path = create_bam(self.dir_path, [self.chromosome], self.sequence, bam_file_name)
        os.chdir(self.dir_path)
        bin_size = len(self.sequence)
        liquidator = blb.BinLiquidator(bin_size = bin_size,
                                       output_directory = os.path.join(self.dir_path, 'output'),
                                       bam_file_path = bam_file_name)

        with tables.open_file(liquidator.counts_file_path) as counts:
            self.assertEqual(1, len(counts.root.files)) # 1 since only a single bam file
            file_record = counts.root.files[0] 
            self.assertEqual(1, file_record['length']) # 1 since only a single read 
            self.assertEqual(1, file_record['key'])

            self.assertEqual(1, len(counts.root.bin_counts)) # 1 since 1 bin accommodates full sequence 
           
            record = counts.root.bin_counts[0]
            self.assertEqual(0, record['bin_number'])
            self.assertEqual("-", record['cell_type'])
            self.assertEqual(self.chromosome, record['chromosome'])
            self.assertEqual(len(self.sequence), record['count']) # count represents how many base pair reads 
                                                                  # intersected the bin

    def test_liquidation_in_long_directory(self):
        self.dir_path = tempfile.mkdtemp(prefix='blt_' + 'a'*16)
        truncated_cell_type = os.path.basename(self.dir_path)[:15]
        self.bam_file_path = create_bam(self.dir_path, [self.chromosome], self.sequence)
        bin_size = len(self.sequence)
        liquidator = blb.BinLiquidator(bin_size = bin_size,
                                       output_directory = os.path.join(self.dir_path, 'output'),
                                       bam_file_path = self.bam_file_path)

        with tables.open_file(liquidator.counts_file_path) as counts:
            self.assertEqual(1, len(counts.root.files)) # 1 since only a single bam file
            file_record = counts.root.files[0] 
            self.assertEqual(1, file_record['length']) # 1 since only a single read 
            self.assertEqual(1, file_record['key'])

            self.assertEqual(1, len(counts.root.bin_counts)) # 1 since 1 bin accommodates full sequence 
           
            record = counts.root.bin_counts[0]
            self.assertEqual(0, record['bin_number'])
            self.assertEqual(truncated_cell_type, record['cell_type'])
            self.assertEqual(self.chromosome, record['chromosome'])
            self.assertEqual(len(self.sequence), record['count']) # count represents how many base pair reads 
                                                                  # intersected the bin

def number_hits(motif_liquidator_output):
    hits = {}
    for line in motif_liquidator_output.split('\n'):
        split = line.split()
        if line.startswith('# total hits: '):
            hits['total'] = int(split[3]) 
        elif line.startswith('# (mapped hit) / (mapped reads) = '):
            hits['mapped'] = int(split[7].split('/')[0])
        elif line.startswith('# (unmapped hit) / (unmapped reads) = '):
            hits['unmapped'] = int(split[7].split('/')[0])
    return hits 

def fimo_style_scores(fimo_style_output):
    scores = []
    for line in fimo_style_output.split('\n'):
        split = line.split('\t')
        if len(split) == 9 and split[0][:1] != '#':
            scores.append({'pattern name': split[0],
                           'sequence name': split[1],
                           'start': int(split[2]),
                           'stop': int(split[3]),
                           'strand': split[4],
                           'score': float(split[5]),
                           'p-value': float(split[6]),
                           'q-value': split[7],
                           'matched sequence': split[8]})
    return scores

class MotifLiquidatorTest(TempDirTest):
    def setUp(self):
        super(MotifLiquidatorTest, self).setUp()
        self.pwm_10a_path = self.create_pwm('10a', ((1,0,0,0),)*10, 18)
        self.pwm_10g_path = self.create_pwm('10g', ((0,0,1,0),)*10, 18)
        self.pwm_10t_path = self.create_pwm('10t', ((0,0,0,1),)*10, 18)
        self.executable_path = os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), 'motif_liquidator')
        self.bam_a = create_bam(self.dir_path, ['chr1'], 10*'A')
        self.expected_perfect_score  = 18.6139  # value calculated by fimo 
        self.expected_perfect_pvalue = 2.43e-06 # value calculated by fimo 
        self.p65_floats = ((0.000000, 0.222222, 0.611111, 0.166667),
                           (0.000000, 0.000000, 0.944444, 0.055556),
                           (0.000000, 0.000000, 1.000000, 0.000000),
                           (0.611111, 0.000000, 0.388889, 0.000000),
                           (0.555556, 0.166667, 0.222222, 0.055556),
                           (0.111111, 0.000000, 0.000000, 0.888889),
                           (0.000000, 0.000000, 0.000000, 1.000000),
                           (0.000000, 0.111111, 0.000000, 0.888889),
                           (0.000000, 1.000000, 0.000000, 0.000000),
                           (0.000000, 1.000000, 0.000000, 0.000000))
        self.p65_pwm_path = self.create_pwm('p65', self.p65_floats, 18)
        self.fasta_p65_match = os.path.join(self.dir_path, 'p65_match.fasta')
        self.first_fasta_name = 'fasta_1'
        with open(self.fasta_p65_match, 'w') as fasta_file:
            fasta_file.write('>%s\nGGGAATTTCC\n' % self.first_fasta_name)

    def create_pwm(self, name, acgt_float_tuple_list, nsites=None):
        path = os.path.join(self.dir_path, name + '_pwm.txt')
        with open(path, 'w') as pwm_file:
            pwm_file.write('MEME version 4\n')
            pwm_file.write('MOTIF ' + name + '\n')
            nsites_text = ''
            if nsites:
                nsites_text = ' nsites= %d' % nsites
            pwm_file.write('letter-probability matrix:' + nsites_text + '\n')
            for acgt in acgt_float_tuple_list:
                pwm_file.write('%f\t%f\t%f\t%f\n' % acgt)
        return path
    
    def test_single_matching_read_fimo_output(self):
        out_bam = os.path.join(self.dir_path, '10a.bam')
        output = subprocess.check_output([self.executable_path, '-p', 'fimo', '-o', out_bam, self.pwm_10a_path, self.bam_a])
        self.assertEqual({'total':1, 'mapped':1, 'unmapped':0}, number_hits(output))
        scores = fimo_style_scores(output)
        self.assertEqual(1, len(scores))
        actual_score = scores[0]
        expected_score = {'pattern name': '10a', 
                          'sequence name': 'read1',
                          'start': 1, 
                          'stop': 10, 
                          'strand': '+',
                          'score': self.expected_perfect_score,
                          'p-value': self.expected_perfect_pvalue, 
                          'q-value': '',
                          'matched sequence': 10*'A'}
        self.assertEqual(expected_score, actual_score)
        sam_out = subprocess.check_output(['samtools', 'view', out_bam])
        self.assertEqual('read1	16	chr1	1	255	10M	*	0	0	AAAAAAAAAA	<<<<<<<<<<	NM:i:0\n',
                         sam_out)

    def test_single_matching_read_mapped_output(self):
        out_bam = os.path.join(self.dir_path, '10a.bam')
        output = subprocess.check_output([self.executable_path, '-p', 'mapped-fimo', '-o', out_bam, self.pwm_10a_path, self.bam_a])
        self.assertEqual({'total':1, 'mapped':1, 'unmapped':0}, number_hits(output))
        scores = fimo_style_scores(output)
        self.assertEqual(1, len(scores))
        actual_score = scores[0]
        expected_score = {'pattern name': '10a', 
                          'sequence name': 'mapped:chr1:read1',
                          'start': 1, 
                          'stop': 10, 
                          'strand': '+',
                          'score': self.expected_perfect_score,
                          'p-value': self.expected_perfect_pvalue, 
                          'q-value': '',
                          'matched sequence': 10*'A'}
        self.assertEqual(expected_score, actual_score)
        sam_out = subprocess.check_output(['samtools', 'view', out_bam])
        self.assertEqual('read1	16	chr1	1	255	10M	*	0	0	AAAAAAAAAA	<<<<<<<<<<	NM:i:0\n',
                         sam_out)

    def test_single_mismatching_read(self):
        out_bam = os.path.join(self.dir_path, '10a.bam')
        output = subprocess.check_output([self.executable_path, '-p', 'mapped-fimo', '-o', out_bam, self.pwm_10g_path, self.bam_a])
        self.assertEqual({'total':0, 'mapped':0, 'unmapped':0}, number_hits(output))
        scores = fimo_style_scores(output)
        self.assertEqual(0, len(scores))
        sam_out = subprocess.check_output(['samtools', 'view', out_bam])
        self.assertEqual('', sam_out)

    def test_single_matching_read_filtered_out_due_to_unmapped_arg(self):
        out_bam = os.path.join(self.dir_path, '10a.bam')
        output = subprocess.check_output([self.executable_path, '-u', '-p', 'mapped-fimo', '-o', out_bam, self.pwm_10a_path, self.bam_a])
        self.assertEqual({'total':0, 'unmapped':0}, number_hits(output))
        scores = fimo_style_scores(output)
        self.assertEqual(0, len(scores))
        sam_out = subprocess.check_output(['samtools', 'view', out_bam])
        self.assertEqual('', sam_out)

    def test_single_matching_read_filtered_out_due_to_bed_region(self):
        out_bam = os.path.join(self.dir_path, '10a.bam')
        region_file = create_single_region_bed_file(self.dir_path, 'chr1', 500, 600)
        output = subprocess.check_output([self.executable_path, '-r', region_file, '-p', 'mapped-fimo', '-o', out_bam, self.pwm_10a_path, self.bam_a])
        self.assertEqual({'total':0, 'mapped':0, 'unmapped':0}, number_hits(output))
        scores = fimo_style_scores(output)
        self.assertEqual(0, len(scores))
        sam_out = subprocess.check_output(['samtools', 'view', out_bam])
        self.assertEqual('', sam_out)

    def helper_check_single_matching_region_read(self, region_file):
        out_bam = os.path.join(self.dir_path, '10a.bam')
        output = subprocess.check_output([self.executable_path, '-r', region_file, '-p', 'mapped-fimo', '-o', out_bam, self.pwm_10a_path, self.bam_a])
        self.assertEqual({'total':1, 'mapped':1, 'unmapped':0}, number_hits(output))
        scores = fimo_style_scores(output)
        self.assertEqual(1, len(scores))
        expected_score = {'pattern name': '10a', 
                          'sequence name': 'mapped:chr1:read1',
                          'start': 1,
                          'stop': 10, 
                          'strand': '+',
                          'score': self.expected_perfect_score,
                          'p-value': self.expected_perfect_pvalue, 
                          'q-value': '',
                          'matched sequence': 10*'A'}
        self.assertEqual(expected_score, scores[0])
        sam_out = subprocess.check_output(['samtools', 'view', out_bam])
        self.assertEqual('read1	16	chr1	1	255	10M	*	0	0	AAAAAAAAAA	<<<<<<<<<<	NM:i:0\n',
                         sam_out)

    def test_single_matching_bed_region_read(self):
        self.helper_check_single_matching_region_read(create_single_region_bed_file(self.dir_path, 'chr1', 0, 100))

    def test_single_matching_gff_region_read(self):
        self.helper_check_single_matching_region_read(create_single_region_gff_file(self.dir_path, 'chr1', 0, 100))

    def test_single_reverse_matching_read(self):
        output = subprocess.check_output([self.executable_path, '-p', 'mapped-fimo', self.pwm_10t_path, self.bam_a])
        self.assertEqual({'total':1, 'mapped':1, 'unmapped':0}, number_hits(output))
        scores = fimo_style_scores(output)
        expected_score = {'pattern name': '10t', 
                          'sequence name': 'mapped:chr1:read1',
                          'start': 1,
                          'stop': 10, 
                          'strand': '-',
                          'score': self.expected_perfect_score,
                          'p-value': self.expected_perfect_pvalue, 
                          'q-value': '',
                          'matched sequence': 10*'T'}
        self.assertEqual(expected_score, scores[0])
        self.assertEqual(1, len(scores))

    def test_single_unmapped_read(self):
        bam_unmapped_a = create_bam(self.dir_path, ['*'], 10*'A', flag=4)
        out_bam = os.path.join(self.dir_path, 'unmapped_10a.bam')
        output = subprocess.check_output([self.executable_path, '-u', '-p', 'mapped-fimo', '-o', out_bam, self.pwm_10a_path, bam_unmapped_a])
        self.assertEqual({'total':1, 'unmapped':1}, number_hits(output))
        scores = fimo_style_scores(output)
        self.assertEqual(1, len(scores))
        actual_score = scores[0]
        actual_name = actual_score.pop('sequence name')
        self.assertIn(actual_name, ['unmapped:*:read1','unmapped:=:read1']) # different versions of samtools seem to handle this differently
        expected_score = {'pattern name': '10a', 
                          'start': 1,
                          'stop': 10, 
                          'strand': '+',
                          'score': self.expected_perfect_score,
                          'p-value': self.expected_perfect_pvalue, 
                          'q-value': '',
                          'matched sequence': 10*'A'}
        self.assertEqual(expected_score, actual_score)
        sam_out = subprocess.check_output(['samtools', 'view', out_bam])
        self.assertIn(sam_out,
                      ['read1	4	*	1	255	10M	*	0	0	AAAAAAAAAA	<<<<<<<<<<	NM:i:0\n',
                       'read1	4	*	1	255	10M	=	0	0	AAAAAAAAAA	<<<<<<<<<<	NM:i:0\n'])

    def test_default_bg_score(self):
        out_fimo_style_path = os.path.join(self.dir_path, 'fimo_out.txt')
        output = subprocess.check_output([self.executable_path, '-o', out_fimo_style_path, self.p65_pwm_path, self.fasta_p65_match])
        self.assertEqual('', output)
        with open(out_fimo_style_path, 'r') as out_fimo_style: 
            scores = fimo_style_scores(out_fimo_style.read())
        self.assertEqual(1, len(scores))
        expected_score = {'pattern name': 'p65',
                          'sequence name': self.first_fasta_name,
                          'start': 1, 
                          'stop': 10, 
                          'strand': '+',
                          'score': 17.3265,    # actual fimo value
                          'p-value': 9.09e-07, # actual fimo value 
                          'q-value': '',
                          'matched sequence': 'GGGAATTTCC'}
        self.assertEqual(expected_score, scores[0])

    def test_not_averaged_nor_normalized_bg_score(self):
        bg_path = os.path.join(self.dir_path, 'not_averaged_nor_normalized_bg.txt')
        with open(bg_path, 'w') as bg: 
            bg.write('A 0.7e-01\nC 0.4e-01\nG 0.5e-01\nT 0.6e-01\n')
        out_fimo_style_path = os.path.join(self.dir_path, 'fimo_out.txt')
        output = subprocess.check_output([self.executable_path, '-b', bg_path, '-o', out_fimo_style_path, self.p65_pwm_path, self.fasta_p65_match])
        self.assertEqual('', output)
        with open(out_fimo_style_path, 'r') as out_fimo_style: 
            scores = fimo_style_scores(out_fimo_style.read())
        self.assertEqual(1, len(scores))
        expected_score = {'pattern name': 'p65',
                          'sequence name': self.first_fasta_name,
                          'start': 1, 
                          'stop': 10, 
                          'strand': '+',
                          'score': 17.4691,    # actual fimo value
                          'p-value': 8.06e-07, # actual fimo value 
                          'q-value': '',
                          'matched sequence': 'GGGAATTTCC'}
        self.assertEqual(expected_score, scores[0])

    @unittest.skip('psite parsing not implemented yet')
    def test_default_sites_score(self):
        pwm_default_sites = self.create_pwm('p65', self.p65_floats)
        out_fimo_style_path = os.path.join(self.dir_path, 'fimo_out.txt')
        output = subprocess.check_output([self.executable_path, '-o', out_fimo_style_path, pwm_default_sites, self.fasta_p65_match])
        with open(out_fimo_style_path, 'r') as out_fimo_style: 
            scores = fimo_style_scores(out_fimo_style.read())
        self.assertEqual('', output)
        expected_score = {'pattern name': 'p65',
                          'sequence name': self.first_fasta_name,
                          'start': 1, 
                          'stop': 10, 
                          'strand': '+',
                          'score': 17.3367,    # actual fimo value
                          'p-value': 9.09e-07, # actual fimo value 
                          'q-value': '',
                          'matched sequence': 'GGGAATTTCC'}
        self.assertEqual(expected_score, scores[0])

    @unittest.skip('psite parsing not implemented yet')
    def test_specified_sites_score(self):
        pwm_6_sites = self.create_pwm('p65', self.p65_floats, 6)
        out_fimo_style_path = os.path.join(self.dir_path, 'fimo_out.txt')
        output = subprocess.check_output([self.executable_path, '-o', out_fimo_style_path, pwm_6_sites, self.fasta_p65_match])
        with open(out_fimo_style_path, 'r') as out_fimo_style: 
            scores = fimo_style_scores(out_fimo_style.read())
        self.assertEqual('', output)
        expected_score = {'pattern name': 'p65',
                          'sequence name': self.first_fasta_name,
                          'start': 1, 
                          'stop': 10, 
                          'strand': '+',
                          'score': 17.2213,    # actual fimo value
                          'p-value': 9.09e-07, # actual fimo value 
                          'q-value': '',
                          'matched sequence': 'GGGAATTTCC'}
        self.assertEqual(expected_score, scores[0])

    @unittest.skip('todo')
    def test_single_strand(self):
        self.assertEqual(1, 2)

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

