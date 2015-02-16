#!/usr/bin/env python

from __future__ import division

import argparse
import itertools
import os
import tables

import common_util as util
from flattener import write_tab_for_all 

def create_ratio_table(h5file):
    class Ratio(tables.IsDescription):
        bin_number = tables.UInt32Col(pos=0)
        chromosome = tables.StringCol(util.chromosome_name_length, pos=1)
        ratio      = tables.Float64Col(pos=2)

    table = h5file.create_table('/', 'ratios', Ratio, 'filtered count / baseline count for each bin')
    table.flush()

    return table

def populate_ratios(baseline_counts, baseline_file_key, filtered_counts, filtered_file_key, ratios):
    for baseline_row, filtered_row in itertools.izip_longest(baseline_counts.where('file_key == %d' % baseline_file_key),
                                                             filtered_counts.where('file_key == %d' % filtered_file_key)):
        if baseline_row is None or filtered_row is None:
            raise RuntimeError('mismatched number of bin counts between baseline and filtered bams -- were inconsistent bin sizes used?')

        assert(baseline_row['bin_number'] == filtered_row['bin_number'])
        assert(baseline_row['chromosome'] == filtered_row['chromosome'])

        ratios.row['bin_number'] = baseline_row['bin_number']
        ratios.row['chromosome'] = baseline_row['chromosome']
        ratios.row['ratio'] = 0 if baseline_row['count'] == 0 else filtered_row['count'] / baseline_row['count']
        ratios.row.append()

    ratios.flush()
    ratios.cols.bin_number.create_csindex()
    ratios.cols.chromosome.create_csindex()
    ratios.cols.ratio.create_csindex()

    sorted_ratios = ratios.copy(newname='sorted_ratios', sortby=ratios.cols.ratio,
                                step=-1, checkCSI=True,
                                title='ratio table sorted in decreasing ratio column order')
    sorted_ratios.flush()

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('motif')

    parser.add_argument('bam_file')

    default_pass_through_args='-b 5000 --skip_plot'
    parser.add_argument('--pass_through', default=default_pass_through_args,
                        help='Override arguments to pass through to bamliquidator_batch.py. Should not include '
                             '-o/--ouput nor positional input .bam arguments. '
                             'Default is "%s".' % default_pass_through_args)

    parser.add_argument('--baseline', help='Baseline counts.h5 file from bamliquidator_batch.py. Bin size should '
                                           'match bin size in --pass_through argument.  If not provided, then '
                                           'it is generated and stored in output directory.')

    default_output_directory='./output'
    parser.add_argument('-o', '--output_directory', default=default_output_directory,
                        help='Directory to store results (and intermediate files). Creates directory if necessary. '
                              'Default is "%s".' % default_output_directory)

    parser.add_argument('-f', '--flatten', action='store_true',
                        help='flatten all HDF5 tables into tab delimited text files in the output directory, one for each '
                              'chromosome (note that HDF5 files can be efficiently queried and used directly -- e.g. please '
                              'see http://www.pytables.org/ for easy to use Python APIs and '
                              'http://www.hdfgroup.org/products/java/hdf-java-html/hdfview/ for an easy to use GUI for '
                              'browsing HDF5 files)')

    args = parser.parse_args()
    
    assert(tables.__version__ >= '3.0.0')
 
    util.mkdir_if_not_exists(args.output_directory)

    # call motif_liquidator on args.bam_file

    # call bamliquidator as necessary on both baseline and filtered
    baseline_counts_file_path = '../../20110819_580_hg19_5k_bamliquidator_output/counts.h5'
    filtered_counts_file_path = '../../20110819_580_hg19_5k_TGGGAA_bamliquidator_output/counts.h5'
    baseline_file_key = 1
    filtered_file_key = 1

    ratio_file_path = os.path.join(args.output_directory, 'ratios.h5')

    with tables.open_file(baseline_counts_file_path, mode = 'r') as baseline_counts_h5:
        with tables.open_file(filtered_counts_file_path, mode = 'r') as filtered_counts_h5:
            with tables.open_file(ratio_file_path, 'w') as ratios_h5:
                ratios = create_ratio_table(ratios_h5)

                baseline_counts = baseline_counts_h5.root.bin_counts
                filtered_counts = filtered_counts_h5.root.bin_counts

                populate_ratios(baseline_counts, baseline_file_key, filtered_counts, filtered_file_key, ratios)

                if args.flatten:
                    write_tab_for_all(ratios_h5, args.output_directory, log=False)

if __name__ == '__main__':
    main()

'''
   The MIT License (MIT) 

   Copyright (c) 2015 John DiMatteo (jdimatteo@gmail.com)

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
