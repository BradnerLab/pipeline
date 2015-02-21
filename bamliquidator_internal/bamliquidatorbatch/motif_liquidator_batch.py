#!/usr/bin/env python

from __future__ import division

import argparse
import itertools
import logging
import os
import subprocess
import sys
import tables
from time import time 

import common_util as util
from flattener import write_tab_for_all 

__version__ = util.version

def create_bin_ratio_table(h5file):
    class Ratio(tables.IsDescription):
        chromosome = tables.StringCol(util.chromosome_name_length, pos=0)
        bin_number = tables.UInt32Col(pos=1)
        ratio      = tables.Float64Col(pos=2)

    table = h5file.create_table('/', 'bin_ratios', Ratio, 'filtered count / baseline count for each bin')
    table.flush()

    return table

def create_region_ratio_table(h5file):
    class Ratio(tables.IsDescription):
        chromosome  = tables.StringCol(util.chromosome_name_length, pos=0)
        region_name = tables.StringCol(64, pos=1)
        start       = tables.UInt64Col(pos=2)
        stop        = tables.UInt64Col(pos=3)
        strand      = tables.StringCol(1, pos=4)
        ratio       = tables.Float64Col(pos=5)

    table = h5file.create_table('/', 'region_ratios', Ratio, 'filtered count / baseline count for each bin')
    table.flush()

    return table

def populate_bin_ratios(baseline_counts, baseline_file_key, filtered_counts, filtered_file_key, ratios):
    for baseline_row, filtered_row in itertools.izip_longest(baseline_counts.where('file_key == %d' % baseline_file_key),
                                                             filtered_counts.where('file_key == %d' % filtered_file_key)):
        if baseline_row is None or filtered_row is None:
            raise RuntimeError('mismatched number of bin counts between baseline and filtered bams -- were inconsistent bin sizes used?')

        assert(baseline_row['chromosome'] == filtered_row['chromosome'])
        assert(baseline_row['bin_number'] == filtered_row['bin_number'])

        ratios.row['chromosome'] = baseline_row['chromosome']
        ratios.row['bin_number'] = baseline_row['bin_number']
        ratios.row['ratio'] = 0 if baseline_row['count'] == 0 else filtered_row['count'] / baseline_row['count']
        ratios.row.append()

    ratios.flush()
    ratios.cols.bin_number.create_csindex()
    ratios.cols.chromosome.create_csindex()
    ratios.cols.ratio.create_csindex()

    sorted_ratios = ratios.copy(newname='sorted_bin_ratios', sortby=ratios.cols.ratio,
                                step=-1, checkCSI=True,
                                title='ratio table sorted in decreasing ratio column order')
    sorted_ratios.flush()

def populate_region_ratios(baseline_counts, baseline_file_key, filtered_counts, filtered_file_key, ratios):
    for baseline_row, filtered_row in itertools.izip_longest(baseline_counts.where('file_key == %d' % baseline_file_key),
                                                             filtered_counts.where('file_key == %d' % filtered_file_key)):
        if baseline_row is None or filtered_row is None:
            raise RuntimeError('mismatched number of region counts between baseline and filtered bams -- were inconsistent regions used?')

        assert(baseline_row['chromosome'] == filtered_row['chromosome'])
        assert(baseline_row['start'] == filtered_row['start'])
        assert(baseline_row['stop'] == filtered_row['stop'])
        assert(baseline_row['strand'] == filtered_row['strand'])

        ratios.row['chromosome'] = baseline_row['chromosome']
        ratios.row['region_name'] = baseline_row['region_name']
        ratios.row['start'] = baseline_row['start']
        ratios.row['stop'] = baseline_row['stop']
        ratios.row['strand'] = baseline_row['strand']
        ratios.row['ratio'] = 0 if baseline_row['count'] == 0 else filtered_row['count'] / baseline_row['count']
        ratios.row.append()

    ratios.flush()
    # if anybody actually queries the hdf5 files directly, they are welcome to add additional indices here
    ratios.cols.ratio.create_csindex()

    sorted_ratios = ratios.copy(newname='sorted_region_ratios', sortby=ratios.cols.ratio,
                                step=-1, checkCSI=True,
                                title='ratio table sorted in decreasing ratio column order')
    sorted_ratios.flush()

def main():
    parser = argparse.ArgumentParser(description='Calculates bin or region ratios between the motif matching read counts and '
                                     'the total read counts.  May use results from prior runs of %s, bamliquidator_batch, or '
                                     'motif_liquidator to skip regenerating intermediate results.' % sys.argv[0])

    parser.add_argument('bam_file', nargs='?', help='The bam file to search for the motif for and run bamliquidator on. '
                                                    'If not provided, then both --prior_filtered_bam and --baseline args '
                                                    'must be provided.')

    motif_or_filtered_grp = parser.add_mutually_exclusive_group(required=True)
    motif_or_filtered_grp.add_argument('-m', '--motif', help='The short string to search for in reads in the bam, e.g. "TGGGAA".  '
                                              'If not provided, then a previously filtered bam file must be provided '
                                              'with the --prior_filtered_bam arg.')
    motif_or_filtered_grp.add_argument('-p', '--prior_filtered_bam',
                                       help='A bam already filtered by the motif, e.g. from a prior %s or motif_liquidator '
                                            'run.' % sys.argv[0])

    bin_or_region_grp = parser.add_mutually_exclusive_group()
    default_bin_size = 5000
    bin_or_region_grp.add_argument('-b', '--bin_size', type=int, default=default_bin_size,
                                   help='Number of base pairs in each bin -- the smaller the bin size the longer the runtime and '
                                        'the larger the data files (default is %d).  Not compatible with --regions_file '
                                        'argument.' % default_bin_size)
    bin_or_region_grp.add_argument('-r', '--regions_file',
                                   help='Region file in either .gff or .bed format.  Not compatible with --bin_size argument.')

    parser.add_argument('-e', '--extension', type=int, default=0,
                        help='Extends reads by n bp (default is 0)')
    parser.add_argument('--sense', default=None, choices=['+', '-', '.'],
                        help="Map to '+' (forward), '-' (reverse) or '.' (both) strands. For gff regions, default is to use "
                             "the sense specified by the gff file; otherwise, default maps to both.")

    parser.add_argument('--baseline', help='Baseline counts.h5 file from prior %s or bamliquidator_batch.py run. Bin size or '
                                           'regions should match the --bins or --regions_file argument.  If not provided, then '
                                           'it is generated and stored in output directory.' % sys.argv[0])
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
    parser.add_argument('--version', action='version', version='%s %s' % (os.path.basename(sys.argv[0]), __version__))

    args = parser.parse_args()

    if args.baseline is None and args.bam_file is None:
        logging.error('Either the --baseline arg or the bam_file positional arg must be provided.')
        sys.exit(1)
    
    assert(tables.__version__ >= '3.0.0')
 
    util.mkdir_if_not_exists(args.output_directory)

    util.configure_logging(args, args.output_directory, quiet=False)
    # todo, log file doesn't contain child process logging

    if args.prior_filtered_bam:
        filtered_bam_path = args.prior_filtered_bam
    else:
        if args.bam_file is None or args.motif is None:
            logging.error('Unless --prior_filtered_bam arg is provided, both the --motif arg and the bam_file positional arg are required.')
            sys.exit(1)
        logging.info('Running motif_liquidator to get filtered bam')
        start = time()
        motif_executable = util.most_appropriate_executable_path('motif_liquidator')
        filtered_bam_path = os.path.join(args.output_directory, 'filtered.' + os.path.basename(args.bam_file))
        subprocess.check_call([motif_executable, args.bam_file, args.motif, filtered_bam_path])
        duration = time() - start
        logging.info('Running motif_liquidator took %f seconds' % duration)

    bamliquidator_batch_args = [os.path.join(os.path.dirname(os.path.realpath(__file__)), 'bamliquidator_batch.py'),
                                '--skip_plot',
                                '--regions_file=%s' % args.regions_file if args.regions_file else '--bin_size=%d' % args.bin_size,
                                '--extension=%s' % args.extension]
    if args.sense:
        bamliquidator_batch_args.append('--sense=%s' % args.sense)
    if args.flatten:
        bamliquidator_batch_args.append('--flatten')

    logging.info('Running bamliquidator_batch on filtered bam')
    start = time()
    filtered_bamliquidator_output_path = os.path.join(args.output_directory, 'filtered_bamliquidator_output')
    subprocess.check_call(bamliquidator_batch_args + ['--output_directory=%s' % filtered_bamliquidator_output_path, filtered_bam_path])
    filtered_counts_file_path = os.path.join(filtered_bamliquidator_output_path, 'counts.h5')
    filtered_file_key = 1
    duration = time() - start
    logging.info('Running bamliquidator_batch on filtered bam took %f seconds' % duration)

    if args.baseline:
        baseline_counts_file_path = args.baseline
    else:
        logging.info('Running bamliquidator_batch on (unfiltered) bam')
        start = time()
        baseline_bamliquidator_output_path = os.path.join(args.output_directory, 'baseline_bamliquidator_output')
        subprocess.check_call(bamliquidator_batch_args + ['--output_directory=%s' % baseline_bamliquidator_output_path, args.bam_file])
        baseline_counts_file_path = os.path.join(baseline_bamliquidator_output_path, 'counts.h5')
        duration = time() - start
        logging.info('Running bamliquidator_batch on (unfiltered) bam took %f seconds' % duration)
    baseline_file_key = 1 # todo: support args.baseline with file_key not equal to 1

    ratio_file_path = os.path.join(args.output_directory, 'ratios.h5')

    with tables.open_file(baseline_counts_file_path, mode = 'r') as baseline_counts_h5:
        with tables.open_file(filtered_counts_file_path, mode = 'r') as filtered_counts_h5:
            with tables.open_file(ratio_file_path, 'w') as ratios_h5:
                logging.info('Calculating ratios between filtered and unfiltered bin counts')
                start = time()
                if args.regions_file:
                    ratios = create_region_ratio_table(ratios_h5)

                    baseline_counts = baseline_counts_h5.root.region_counts
                    filtered_counts = filtered_counts_h5.root.region_counts

                    populate_region_ratios(baseline_counts, baseline_file_key, filtered_counts, filtered_file_key, ratios)
                else:
                    ratios = create_bin_ratio_table(ratios_h5)

                    baseline_counts = baseline_counts_h5.root.bin_counts
                    filtered_counts = filtered_counts_h5.root.bin_counts

                    populate_bin_ratios(baseline_counts, baseline_file_key, filtered_counts, filtered_file_key, ratios)
                    
                duration = time() - start
                logging.info('Calculating ratios took %f seconds' % duration)

                if args.flatten:
                    logging.info('Flattening ratios tables')
                    start = time()
                    write_tab_for_all(ratios_h5, args.output_directory, log=False)
                    duration = time() - start
                    logging.info('Flattening took %f seconds' % duration)

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
