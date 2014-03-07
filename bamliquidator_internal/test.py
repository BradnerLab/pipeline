#!/usr/bin/env python

# unfortunately these tests are currently defunct, as they weren't updated when the normalization
# method was switched over to samtoolds idxstats

from os import sys, path
sys.path.append(path.dirname(path.dirname(path.abspath(__file__)))) # so we can import parent directory

import bamliquidator_batch as blb
import normalize_plot_and_summarize as nps

import unittest
import tables

million = 10**6 

class TestNormalizationAndPercentiles(unittest.TestCase):

    def setUp(self):
        # create an in memory H5 file
        self.h5file = tables.open_file("test.h5", "a", driver="H5FD_CORE",
                                       driver_core_backing_store=0)
        self.genome_1 = "genome1"
        self.genome_2 = "genome2"
        self.genome_3 = "genome3"
        self.cell_type_1 = "celltype1"
        self.cell_type_2 = "celltype2"
        self.counts = blb.create_count_table(self.h5file)
        self.normalized_counts = nps.create_normalized_counts_table(self.h5file)

        self.bin_size = 1000 

        # genome 1, cell type 1:
        for i in xrange(10):
            self.counts.row["bin_number"] = i 
            self.counts.row["cell_type"] = self.cell_type_1 
            self.counts.row["chromosome"] = "chr1" 
            self.counts.row["file_name"] = self.genome_1
            self.counts.row["count"] = million
            self.counts.row.append()
        for i in xrange(6):
            self.counts.row["bin_number"] = i 
            self.counts.row["cell_type"] = self.cell_type_1
            self.counts.row["chromosome"] = "chr2" 
            self.counts.row["file_name"] = self.genome_1
            self.counts.row["count"] = i*million 
            self.counts.row.append()

        # genome 2, cell type 1:
        for i in xrange(16):
            if i < 10:
                self.counts.row["bin_number"] = i 
                self.counts.row["chromosome"] = "chr1" 
            else:
                self.counts.row["bin_number"] = i-10 
                self.counts.row["chromosome"] = "chr2" 
            self.counts.row["cell_type"] = self.cell_type_1
            self.counts.row["file_name"] = self.genome_2
            self.counts.row["count"] = million 
            self.counts.row.append()
            
        # genome 3, cell type 2:
        for i in xrange(16):
            if i < 10:
                self.counts.row["bin_number"] = i 
                self.counts.row["chromosome"] = "chr1" 
            else:
                self.counts.row["bin_number"] = i-10 
                self.counts.row["chromosome"] = "chr2" 
            self.counts.row["cell_type"] = self.cell_type_2
            self.counts.row["file_name"] = self.genome_3
            self.counts.row["count"] = 2*million 
            self.counts.row.append()

        self.h5file.flush()


    def tearDown(self):
        self.h5file.close()

    def test_normalization_of_single_genome(self):
        nps.populate_normalized_counts(self.normalized_counts, self.counts, self.genome_1, self.bin_size)

        self.assertEqual(16, len(self.normalized_counts))

        # genome 1:
        # 10 million reads in chr1, 0+1+2+3+4+5 = 15 million reads in chr2, 25 million reads total
        # normalized units are reads per million per basepair
        # e.g. chr1 bin 0 normalized is:  1,000,000 / 1,000 / 25 =  40 rpm/bp

        normalized_count_per_million = 40 
        
        for i in xrange(16):
            self.assertEqual(self.counts[i]["bin_number"], self.normalized_counts[i]["bin_number"])
            self.assertEqual(self.counts[i]["cell_type"], self.normalized_counts[i]["cell_type"])
            self.assertEqual(self.counts[i]["chromosome"], self.normalized_counts[i]["chromosome"])
            self.assertEqual(self.counts[i]["file_name"], self.normalized_counts[i]["file_name"])

        for i in xrange(10):
            self.assertAlmostEqual(normalized_count_per_million, self.normalized_counts[i]["count"])
         
        for i in xrange(10, 16):
            self.assertAlmostEqual((i-10)*normalized_count_per_million, self.normalized_counts[i]["count"])

    def test_file_names_in_cell_type(self):
        file_names = nps.file_names_in_cell_type(self.normalized_counts, self.cell_type_1)
        self.assertEqual(0, len(file_names))

        nps.populate_normalized_counts(self.normalized_counts, self.counts, self.genome_1, self.bin_size)
        nps.populate_normalized_counts(self.normalized_counts, self.counts, self.genome_2, self.bin_size)
        nps.populate_normalized_counts_for_cell_type(self.normalized_counts, self.cell_type_1,
                                                     [self.genome_1, self.genome_2])
        file_names = nps.file_names_in_cell_type(self.normalized_counts, self.cell_type_1)
        self.assertEqual(2, len(file_names))
        self.assertTrue(   (self.genome_1 == file_names[0] and self.genome_2 == file_names[1])
                        or (self.genome_1 == file_names[1] and self.genome_2 == file_names[0]))

    def test_normalization_of_cell_type(self):
        nps.populate_normalized_counts(self.normalized_counts, self.counts, self.genome_1, self.bin_size)
        nps.populate_normalized_counts(self.normalized_counts, self.counts, self.genome_2, self.bin_size)
        nps.populate_normalized_counts_for_cell_type(self.normalized_counts, self.cell_type_1,
                                                     [self.genome_1, self.genome_2])

        # genome 1: 40 rpm/bp for each 1 million reads in a bin (see notes from above test) 

        # genome 2: 1 million reads in each bin, 16 bins, so 16 million reads total
        # so reads per million per basepair is 1,000,000 / 1,000 / 16 = 62.5 rpm/bp

        # average across cell types just sum and divide by total number of genomes in the cell type
        # e.g. so chr1 bin [0-9] has (40 + 62.5)/2 = 51.25

        self.assertEqual(16*3, len(self.normalized_counts))

        condition = "(file_name == '*') & (cell_type == '%s')" % self.cell_type_1
        self.assertEqual(16, len(self.normalized_counts.read_where(condition)))

        for i, row in enumerate(self.normalized_counts.where(condition + " & (chromosome == 'chr1')")):
            self.assertEqual("*", row["file_name"])
            self.assertEqual(i, row["bin_number"])
            self.assertEqual(self.cell_type_1, row["cell_type"])
            self.assertEqual("chr1", row["chromosome"])
            self.assertAlmostEqual(51.25, row["count"])

        for i, row in enumerate(self.normalized_counts.where(condition + " & (chromosome == 'chr2')")):
            self.assertEqual("*", row["file_name"])
            self.assertEqual(i, row["bin_number"])
            self.assertEqual(self.cell_type_1, row["cell_type"])
            self.assertEqual("chr2", row["chromosome"])
            if i == 1:
                self.assertAlmostEqual(51.25, row["count"])

    def test_normalized_counts_add_up_to_the_same(self):
        nps.populate_normalized_counts(self.normalized_counts, self.counts, self.genome_1, self.bin_size)
        nps.populate_normalized_counts(self.normalized_counts, self.counts, self.genome_2, self.bin_size)
        nps.populate_normalized_counts_for_cell_type(self.normalized_counts, self.cell_type_1,
                                                     [self.genome_1, self.genome_2])
        
        def count_sum(file_name):
            return sum(self.normalized_counts.read_where("file_name == '%s'" % file_name, field="count"))

        genome_1_sum = count_sum(self.genome_1) 
        genome_2_sum = count_sum(self.genome_2) 
        cell_type_sum = count_sum("*") 
        self.assertAlmostEqual(genome_1_sum, genome_2_sum, cell_type_sum)


    def test_percentiles_of_single_genome(self):
        '''      Percentiles
          200                  o
          160                 o 
   rpm/bp 120                o  
           80               o   
           40    ooooooooooo      
            0   o                 
               -|--|---|---|---|
                0  25  50  75  100
        '''

        nps.populate_normalized_counts(self.normalized_counts, self.counts, self.genome_1, self.bin_size)
        nps.populate_percentiles(self.normalized_counts, self.cell_type_1, self.genome_1)
        self.assertEqual(16, len(self.normalized_counts))

        for i in xrange(10):
            self.assertEqual(self.normalized_counts[i]["percentile"], 40) 
        
        self.assertEqual(self.normalized_counts[10]["percentile"], 0)
        self.assertEqual(self.normalized_counts[11]["percentile"], 40)
        self.assertEqual(self.normalized_counts[15]["percentile"], 100)


    def test_percentiles_of_cell_type(self):
        '''      Percentiles
          130                  o
          110                 o 
           90                o  
   rpm/bp  70               o   
           50    ooooo*ooooo      
           30   o               

               -|--|---|---|---|
                0  25  50  75  100
        '''
        nps.populate_normalized_counts(self.normalized_counts, self.counts, self.genome_1, self.bin_size)
        nps.populate_normalized_counts(self.normalized_counts, self.counts, self.genome_2, self.bin_size)
        nps.populate_normalized_counts_for_cell_type(self.normalized_counts, self.cell_type_1,
                                                     [self.genome_1, self.genome_2])
        nps.populate_percentiles(self.normalized_counts, self.cell_type_1, "*")

        self.assertEqual(16, len(self.normalized_counts.read_where("file_name == '*'")))

        for row in self.normalized_counts.where("file_name == '*'"):
            if 50 < row["count"] < 70:
                self.assertTrue(35 < row["percentile"] < 50) 
            elif row["count"] < 40:
                self.assertEqual(0, row["percentile"])
            elif row["count"] > 125:
                self.assertEqual(100, row["percentile"])


class TestSummary(unittest.TestCase):

    def setUp(self):
        # create an in memory H5 file
        self.h5file = tables.open_file("test.h5", "a", driver="H5FD_CORE",
                                       driver_core_backing_store=0)
        self.genome_1 = "genome1"
        self.genome_2 = "genome2"
        self.genome_3 = "genome3"
        self.cell_type_1 = "celltype1"
        self.cell_type_2 = "celltype2"
        self.normalized_counts = nps.create_normalized_counts_table(self.h5file)
        self.summary = nps.create_summary_table(self.h5file)

        # first just set everything to the 50th percentile
        for file_name in [self.genome_1, self.genome_2, self.genome_3, self.cell_type_1, self.cell_type_2]:
            for i in xrange(16):
                if i < 10:
                    self.normalized_counts.row["bin_number"] = i 
                    self.normalized_counts.row["chromosome"] = "chr1" 
                else:
                    self.normalized_counts.row["bin_number"] = i-10 
                    self.normalized_counts.row["chromosome"] = "chr2" 
                if file_name in [self.cell_type_1, self.cell_type_2]:
                    self.normalized_counts.row["cell_type"] = file_name 
                    self.normalized_counts.row["file_name"] = "*" 
                else:
                    self.normalized_counts.row["cell_type"] = (
                        self.cell_type_2 if file_name == self.genome_3 else self.cell_type_1)
                    self.normalized_counts.row["file_name"] = file_name 
                self.normalized_counts.row["percentile"] = 50 
                self.normalized_counts.row["count"] = -1 
                self.normalized_counts.row.append()
        self.h5file.flush()

        # then bin 0 of chromosome 1 in cell type 1 to 100th percentile,
        #  and bin 1 of chromosome 2 in cell type 1 to 0th percentile 
        for row in self.normalized_counts: 
            if row["chromosome"] == "chr1" and row["cell_type"] == self.cell_type_1:
                if row["bin_number"] == 0:
                    row["percentile"] = 100
                elif row["bin_number"] == 1:
                    row["percentile"] = 0 
                row.update()

    def tearDown(self):
        self.h5file.close()


    def test_summary(self):
        nps.populate_summary(self.summary, self.normalized_counts, "chr1")
        nps.populate_summary(self.summary, self.normalized_counts, "chr2")
        self.assertEqual(16, len(self.summary.read()))
        
        for i, row in enumerate(self.summary):
            self.assertEqual('chr1' if i < 10 else 'chr2', row["chromosome"])
            if row["chromosome"] == "chr1":
                self.assertEqual(i, row["bin_number"])
                if row["bin_number"] == 0:
                    self.assertAlmostEqual(75, row["avg_cell_type_percentile"])
                    self.assertEqual(1, row["cell_types_gte_95th_percentile"])
                    self.assertEqual(1, row["cell_types_lt_95th_percentile"])
                    self.assertEqual(2, row["lines_gte_95th_percentile"])
                    self.assertEqual(1, row["lines_lt_95th_percentile"])
                    self.assertEqual(2, row["cell_types_gte_5th_percentile"])
                    self.assertEqual(0, row["cell_types_lt_5th_percentile"])
                    self.assertEqual(3, row["lines_gte_5th_percentile"])
                    self.assertEqual(0, row["lines_lt_5th_percentile"])
                    continue
                elif row["bin_number"] == 1:
                    self.assertAlmostEqual(25, row["avg_cell_type_percentile"])
                    self.assertEqual(0, row["cell_types_gte_95th_percentile"])
                    self.assertEqual(2, row["cell_types_lt_95th_percentile"])
                    self.assertEqual(0, row["lines_gte_95th_percentile"])
                    self.assertEqual(3, row["lines_lt_95th_percentile"])
                    self.assertEqual(1, row["cell_types_gte_5th_percentile"])
                    self.assertEqual(1, row["cell_types_lt_5th_percentile"])
                    self.assertEqual(1, row["lines_gte_5th_percentile"])
                    self.assertEqual(2, row["lines_lt_5th_percentile"])
                    continue 
            self.assertEqual(50, row["avg_cell_type_percentile"])
            self.assertEqual(0, row["cell_types_gte_95th_percentile"])
            self.assertEqual(2, row["cell_types_lt_95th_percentile"])
            self.assertEqual(0, row["lines_gte_95th_percentile"])
            self.assertEqual(3, row["lines_lt_95th_percentile"])
            self.assertEqual(2, row["cell_types_gte_5th_percentile"])
            self.assertEqual(0, row["cell_types_lt_5th_percentile"])
            self.assertEqual(3, row["lines_gte_5th_percentile"])
            self.assertEqual(0, row["lines_lt_5th_percentile"])

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

