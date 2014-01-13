#!/usr/bin/env python

import tables
import argparse

# creates empty file, overwriting any prior existing files
def create_count_table(h5file):
    class BinCount(tables.IsDescription):
        bin_number = tables.UInt32Col(    pos=0)
        cell_type  = tables.StringCol(16, pos=1)
        chromosome = tables.StringCol(16, pos=2)
        count      = tables.UInt64Col(    pos=3)
        file_name  = tables.StringCol(64, pos=4)

    table = h5file.create_table("/", "counts", BinCount, "bin counts")

    table.flush()

    return table

'''
todo

 * testing
     * check why top couple bins don't match in plot
     * check some unnormalized count rows
     * check some normalized count rows
     * check some summary table rows
 * csv file replacement
     * add a sample script to generate csv files from summary table rows,
       with comments so that people could create similar scripts and hopefully
       not use csv files when pytables is a nicer interface
 * master script
     * do 3 steps outlined in comment
 * add arguments to this script
     * output directory
        * probably exit with error code if directory already exists
          (I don't want to overwrite any files)
     * number of bins
     * option to add files to prior counts instead of start from scratch (maybe by supplying path to prior h5 file?)
     * directory to search for individual bam files, or a specific bam file
 * intended usage:
   make
   -- compiles misc code in bamliquidator_internal and places bamliquidator in current directory (and bamliquidator_batch under internal directory)
'''

def main():
    parser = argparse.ArgumentParser(description='Count the number of base pair reads in each bin of each chromosome '
                                                 'in the bam file(s) at the given directory, and then normalize, plot, '
                                                 'and summarize the counts in the output directory.  For additional '
                                                 'help, please see the wiki at '
                                                 ' https://github.com/BradnerLab/pipeline/wiki , review related '
                                                 'support tickets, and/or create a new ticket at '
                                                 'https://github.com/BradnerLab/pipeline/issues?labels=support .')
    parser.add_argument('--output_directory', default='output',
                        help='Directory to create and output the h5 and/or html files to (aborts if already exists).')
    parser.add_argument('--bin_counts_file',
                        help='HDF5 counts file from a prior run to be appended to.  If unspecified, defaults to '
                             'creating a new file "bin_counts.h5" in the output directory.')
    parser.add_argument('--ucsc_chrom_sizes',
                        help='Tab delimitted text file with the first column naming the chromosome (e.g. chr1), the '
                             'third column naming the genome type (e.g. mm8), and the fifth column naming the number '
                             'of base pairs in the reference chromosome.')
    parser.add_argument('bam_file_path', 
                        help='The directory to recursively search for .bam files for counting.  Every .bam file must '
                             'have a corresponding .bai file at the same location.  To count just a single file, '
                             'provide the .bam file path instead of a directory.  The parent directory of each .bam '
                             'file is interpretted as the cell type.   If your .bam files are not in this format, '
                             'please consider creating a directory of sym links to your actual .bam and .bai files. '
                             'If the .bam file already has 1 or more reads in the HDF5 counts file, then the .bam file '
                             'is skipped.')
    args = parser.parse_args()

    assert(tables.__version__ >= '3.0.0')

    os.mkdir(args.output_directory)

    if args.bin_counts_file is None:
        counts_file = tables.open_file(output_directory + "/bin_counts.h5", mode = "w",
                                       title = "bam liquidator genome bin read counts")
        counts = create_count_table(counts_file)
    else:
        counts_file = tables.open_file(args.bin_counts_h5_file, "r+")
        counts = counts_file.root.counts

    
   
    #create_count_table("bin_counts.h5")

if __name__ == "__main__":
    main()

'''
   The MIT License (MIT) 

   Copyright (c) 2013 John DiMatteo 

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
