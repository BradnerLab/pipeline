#!/usr/bin/env python

from bamliquidator_internal.normalize_plot_and_summarize import normalize_plot_and_summarize

import argparse
import os
import subprocess
import tables
import datetime
from time import time 

# creates empty file, overwriting any prior existing files
def create_count_table(h5file):
    class BinCount(tables.IsDescription):
        bin_number = tables.UInt32Col(    pos=0)
        cell_type  = tables.StringCol(16, pos=1)
        chromosome = tables.StringCol(16, pos=2)
        count      = tables.UInt64Col(    pos=3)
        # todo: rename file_name to genome or line
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
'''

def all_bam_files_in_directory(bam_directory):
    bam_files = []
    for dirpath, _, files in os.walk(bam_directory, followlinks=True):
        for file_ in files:
            if file_.endswith(".bam"):
                bam_files.append(os.path.join(dirpath, file_))
    return bam_files

def bam_files_with_no_counts(counts, bam_files):
    with_no_counts = []

    for bam_file in bam_files:
        count_found = False
        condition = "file_name == '%s'" % os.path.basename(bam_file)

        for _ in counts.where(condition):
            count_found = True

        if not count_found:
            with_no_counts.append(bam_file)

    return with_no_counts

def liquidate(bam_files, output_directory, ucsc_chrom_sizes, bin_size, bin_counts_file_path, executable_path):
    for i, bam_file in enumerate(bam_files):
        print "Liquidating %s (file %d of %d, %s)" % (
            bam_file, i+1, len(bam_files), datetime.datetime.now().strftime('%H:%M:%S'))
        cell_type = os.path.basename(os.path.dirname(bam_file))
        args = [executable_path, cell_type, str(bin_size), ucsc_chrom_sizes, bam_file, bin_counts_file_path]
        start = time()
        return_code = subprocess.call(args)
        end = time()
        print "Liquidation completed in %f seconds" % (end - start)

        if return_code != 0:
            print "%s failed with exit code %d" % (executable_path, return_code)
            exit(return_code)

    counts_file = tables.open_file(bin_counts_file_path, mode = "r")
    counts = counts_file.root.counts

    normalize_plot_and_summarize(counts, output_directory, bin_size) 
    counts_file.close()

def main():
    parser = argparse.ArgumentParser(description='Count the number of base pair reads in each bin of each chromosome '
                                                 'in the bam file(s) at the given directory, and then normalize, plot, '
                                                 'and summarize the counts in the output directory.  For additional '
                                                 'help, please see https://github.com/BradnerLab/pipeline/wiki')
    parser.add_argument('--output_directory', default='output',
                        help='Directory to create and output the h5 and/or html files to (aborts if already exists). '
                             'Default is "./output".')
    parser.add_argument('--bin_counts_file',
                        help='HDF5 counts file from a prior run to be appended to.  If unspecified, defaults to '
                             'creating a new file "bin_counts.h5" in the output directory.')
    parser.add_argument('--bin_size', type=int, default=100000,
                        help="Number of base pairs in each bin -- the smaller the bin size the longer the runtime and "
                             "the larger the data files (default is 100000).")
    parser.add_argument('ucsc_chrom_sizes',
                        help='Tab delimited text file with the first column naming the chromosome (e.g. chr1), the '
                             'third column naming the genome type (e.g. mm8), and the fifth column naming the number '
                             'of base pairs in the reference chromosome.')
    parser.add_argument('bam_file_path', 
                        help='The directory to recursively search for .bam files for counting.  Every .bam file must '
                             'have a corresponding .bai file at the same location.  To count just a single file, '
                             'provide the .bam file path instead of a directory.  The parent directory of each .bam '
                             'file is interpreted as the cell type (e.g. mm1s might be an appropriate directory '
                             'name).  Bam files in the same directory are grouped together for plotting. Plots use '
                             'normalized counts, such that all .bam files in the same directory have bin '
                             'counts that add up to 1 for each chromosome.  The .bam file name is also required to '
                             'contain the genome type so that the corresponding entries in the ucsc_chrom_sizes file '
                             'can be used.  If your .bam files are not in this directory format, please consider '
                             'creating a directory of sym links to your actual .bam and .bai files. If the .bam file '
                             'already has 1 or more reads in the HDF5 counts file, then that .bam file is skipped.')
    args = parser.parse_args()

    assert(tables.__version__ >= '3.0.0')

    executable_path = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                   "bamliquidator_internal", "bamliquidator_batch")
    if not os.path.isfile(executable_path):
        print "%s is missing -- try cd'ing into the directory and running 'make'" % executable_path
        exit(1)

    os.mkdir(args.output_directory)

    if args.bin_counts_file is None:
        args.bin_counts_file = os.path.join(args.output_directory, "bin_counts.h5")
        counts_file = tables.open_file(args.bin_counts_file, mode = "w",
                                       title = "bam liquidator genome bin read counts")
        counts = create_count_table(counts_file)
    else:
        counts_file = tables.open_file(args.bin_counts_file, "r+")
        counts = counts_file.root.counts

    if os.path.isdir(args.bam_file_path):
        bam_files = all_bam_files_in_directory(args.bam_file_path)
    else:
        bam_files = [args.bam_file_path]
   
    bam_files = bam_files_with_no_counts(counts, bam_files)

    counts_file.close() # The bamliquidator_internal/bamliquidator_batch will open this file and modify it,
                        # so it is best that we not hold an out of sync reference to it

    liquidate(bam_files, args.output_directory, args.ucsc_chrom_sizes, args.bin_size, args.bin_counts_file,
              executable_path)

if __name__ == "__main__":
    main()

'''
   The MIT License (MIT) 

   Copyright (c) 2013 John DiMatteo (jdimatteo@gmail.com)

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
