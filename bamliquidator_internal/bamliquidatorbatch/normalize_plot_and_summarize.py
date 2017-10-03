#!/usr/bin/env python

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

from __future__ import division

import sys
import os
import argparse
import tables
import scipy.stats as stats
import collections
import logging

bokeh_install_command = 'sudo pip install bokeh==0.9.3 "openpyxl>=1.6.1,<2.0.0"'

try:
    try:
        import bokeh.plotting as bp
    except:
        bokeh_import_error = 'Bokeh module not found; consider running the following command to install:\n%s' % (
                bokeh_install_command)
        raise

    try:
        # bamliquidatorbatch was originally developed for bokeh 0.4.4, but
        # around version 0.7 the bokeh API had some breaking changes (like requiring the explicit use of vplot).
        # I confirmed things work with version 0.9.3, and I hope future versions maintain compatability.
        bp.vplot
    except:
        bokeh_import_error = 'Bokeh version is incompatible; consider running the following command to upgrade:\n%s' % (
                bokeh_install_command)
        raise
except:
    bp = None 

# note that my initial version didn't do any flush calls, which lead to bogus rows being added
# to the normalized_counts table (which was evident when the normalized counts <= 95 + > 95 didn't add up right).
# -- I should probably look into why flush was necessary and/or file a bug with pytables

# I also found that create_index doesn't always work (this was causing where statements to not work)
# -- I don't know if this was my fault or a bug in pytables, but I just always use create_csindex instead

chromosome_name_length = 64 # Includes 1 for null terminator, so really max of 63 characters.
                            # Note that changing this value requires updating C++ code as well.

def delete_all_but_bin_counts_and_files_table(h5file):
    for table in h5file.root:
        if table.name != "bin_counts" and table.name != "files" and table.name != "file_names":
            for index in table.colindexes.values():
                index.column.remove_index()
            table.remove()

def create_normalized_counts_table(h5file):
    class BinCount(tables.IsDescription):
        bin_number = tables.UInt32Col(    pos=0)
        cell_type  = tables.StringCol(16, pos=1)
        chromosome = tables.StringCol(chromosome_name_length, pos=2)
        count      = tables.Float64Col(   pos=3)
        percentile = tables.Float64Col(   pos=4)
        file_key   = tables.UInt32Col(    pos=5)

    table = h5file.create_table("/", "normalized_counts", BinCount, "normalized bin counts")

    table.flush()

    return table

def all_cell_types(counts):
    types = set()

    for row in counts:
        types.add(row["cell_type"])

    return types 

def all_chromosomes(counts):
    chromosomes = collections.OrderedDict() 

    for row in counts:
        chromosomes[row["chromosome"]] = None

    return chromosomes.keys() 

# todo: if this used the files table and we added the cell_type to the files table, this would be much faster,
#       but it is probably necessary to leave cell_type in counts table as well (for queries)
file_keys_memo = {}
def file_keys(counts, cell_type):
    if not cell_type in file_keys_memo:
        file_keys = set() 
       
        logging.debug("Getting file keys for cell type %s", cell_type)
        for row in counts.where("cell_type == '%s'" % cell_type):
            file_keys.add(row["file_key"])

        file_keys_memo[cell_type] = file_keys

        logging.debug("memoizing files for %s: %s", cell_type, str(file_keys_memo[cell_type]))
        
    return file_keys_memo[cell_type] 

def plot_summaries(output_directory, normalized_counts, chromosomes):
    bp.output_file(output_directory + "/summary.html", mode="cdn")
    plots = []
    for chromosome in chromosomes:
        plot = plot_summary(normalized_counts, chromosome)
        if plot:
            plots.append(plot)
    bp.save(bp.vplot(*plots))

def plot_summary(normalized_counts, chromosome):
    logging.debug(" - plotting %s summary", chromosome)

    condition = "(file_key == 0) & (chromosome == '%s')" % chromosome

    chromosome_count_by_bin = collections.defaultdict(int) 
    for row in normalized_counts.where(condition):
        chromosome_count_by_bin[row["bin_number"]] += row["count"]
  
    num_bins = len(chromosome_count_by_bin)
    if num_bins < 2:
        logging.info("-- skipping plotting %s because not enough bins (only %d)", chromosome, num_bins)
        return

    plot = bp.figure()
    plot.title = chromosome + " counts per bin across all bam files"
    plot.scatter(chromosome_count_by_bin.keys(), chromosome_count_by_bin.values())
    return plot

def plot(output_directory, normalized_counts, chromosome, cell_types):
    bp.output_file(output_directory + "/" + chromosome + ".html", mode="cdn")
    plots = []

    summary = plot_summary(normalized_counts, chromosome)
    if summary:
        plots.append(summary)
    else:
        # if summary can't be plotted, then rest probably can't be plotted either,
        # so return without even saving the file (the file is never created on disk if not saved)
        return
    for cell_type in cell_types:
        logging.debug(" - plotting %s", cell_type)

        bin_number = [] 
        count = [] 
        
        condition = "(file_key == 0) & (chromosome == '%s') & (cell_type == '%s')" % (chromosome, cell_type)

        for row in normalized_counts.where(condition):
            bin_number.append(row["bin_number"])
            count.append(row["count"])

        plot = bp.figure()
        plot.title = "%s counts per bin" % cell_type
        plot.scatter(bin_number, count)
        plots.append(plot)
    bp.save(bp.vplot(*plots))

def populate_normalized_counts(normalized_counts, counts, file_key, bin_size, files):
    total_count = length_for_file_key(files, file_key)

    '''
    Excerpt from Feb 13, 2014 email from Charles Lin:

    We typically report read density in units of reads per million per basepair

    bamliquidator reports counts back in total read positions per bin.  To convert that 
    into reads per million per basepair, we first need to divide by the total million 
    number of reads in the bam.  Then we need to divide by the size of the bin

    So for instance if you have a 1kb bin and get 2500 counts from a bam with 30 million
    reads you would calculate density as 2500/1000/30 = 0.083rpm/bp
    '''
    factor = (1 / bin_size) * (1 / (total_count / 10**6))

    for count_row in counts.where("file_key == %d" % file_key):
        normalized_counts.row["bin_number"] = count_row["bin_number"]
        normalized_counts.row["cell_type"] = count_row["cell_type"] 
        normalized_counts.row["chromosome"] = count_row["chromosome"] 
        assert file_key == count_row["file_key"]
        normalized_counts.row["file_key"] = file_key
        normalized_counts.row["count"] = count_row["count"] * factor 
        normalized_counts.row["percentile"] = -1
        normalized_counts.row.append()

    normalized_counts.flush()
  

def length_for_file_key(files, file_key):
    file_rows = files.read_where("key == %d" % file_key)
    assert len(file_rows) == 1
    return file_rows[0]["length"]

def normalize_regions(region_counts, files):
    logging.info("Normalizing")
     
    file_key = None

    for row in region_counts:
        if row["file_key"] != file_key:
            file_key = row["file_key"]
            total_count = length_for_file_key(files, file_key)
        
        region_size = row["stop"] - row["start"]
        factor = (1 / region_size) * (1 / (total_count / 10**6))

        row["normalized_count"] = row["count"] * factor 
        row.update()

    region_counts.flush()

# leave off file_key argument to calculate percentiles for the cell_type averaged normalized counts
def populate_percentiles(normalized_counts, cell_type, file_key = 0):
    bin_numbers = []
    normalized_count_list = []

    condition = "(cell_type == '%s') & (file_key == %d)" % (cell_type, file_key)

    for row in normalized_counts.where(condition):
        bin_numbers.append(row["bin_number"])
        normalized_count_list.append(row["count"])

    percentiles = (stats.rankdata(normalized_count_list) - 1) / (len(normalized_count_list)-1) * 100
    # percentiles calculated in bulk as suggested at 
    # http://grokbase.com/t/python/python-list/092235vj27/faster-scipy-percentileofscore

    for i, row in enumerate(normalized_counts.where(condition)):
        assert bin_numbers[i] == row["bin_number"]
        row["percentile"] = percentiles[i]
        row.update()
    normalized_counts.flush()

# the cell type normalized counts are the averages of the genomes in the cell type
def populate_normalized_counts_for_cell_type(normalized_counts, cell_type, file_keys):
    processed_a_single_file = False
    chromosome_to_summed_counts = collections.OrderedDict() 

    for file_key in file_keys:
        condition = "(file_key == %d) & (cell_type == '%s')" % (file_key, cell_type)
        for row in normalized_counts.where(condition):
            if processed_a_single_file:
                chromosome_to_summed_counts[row["chromosome"]][row["bin_number"]] += row["count"]
            else:
                if not chromosome_to_summed_counts.has_key(row["chromosome"]):
                    chromosome_to_summed_counts[row["chromosome"]] = []
                chromosome_to_summed_counts[row["chromosome"]].append(row["count"])
        processed_a_single_file = True
            
    cell_type_condition = "(file_key == 0) & (cell_type == '%s')" % cell_type

    len_file_keys = len(file_keys)

    for chromosome, summed_counts in chromosome_to_summed_counts.iteritems():
        for i, summed_count in enumerate(summed_counts):
            normalized_counts.row["bin_number"] = i
            normalized_counts.row["cell_type"] = cell_type 
            normalized_counts.row["chromosome"] = chromosome 
            normalized_counts.row["file_key"] = 0 
            normalized_counts.row["count"] = chromosome_to_summed_counts[chromosome][i] / len_file_keys
            normalized_counts.row["percentile"] = -1
            normalized_counts.row.append()

    normalized_counts.flush()

def create_summary_table(h5file):
    class Summary(tables.IsDescription):
        bin_number = tables.UInt32Col(                    pos=0)
        chromosome = tables.StringCol(chromosome_name_length, pos=2)
        avg_cell_type_percentile = tables.Float64Col(     pos=1)
        cell_types_gte_95th_percentile = tables.UInt32Col(pos=2)
        cell_types_lt_95th_percentile = tables.UInt32Col( pos=3)
        lines_gte_95th_percentile = tables.UInt32Col(     pos=4)
        lines_lt_95th_percentile = tables.UInt32Col(      pos=5)
        cell_types_gte_5th_percentile = tables.UInt32Col( pos=6)
        cell_types_lt_5th_percentile = tables.UInt32Col(  pos=7)
        lines_gte_5th_percentile = tables.UInt32Col(      pos=8)
        lines_lt_5th_percentile = tables.UInt32Col(       pos=9)

    table = h5file.create_table("/", "summary", Summary, "bin count summary")

    table.flush()

    return table
   

def populate_summary(summary, normalized_counts, chromosome):
    high = 95 # 95th percentile
    low  = 5  # 5th percentile

    summed_cell_type_percentiles_by_bin = collections.defaultdict(float) 
    cell_types_gte_high_percentile_by_bin = collections.defaultdict(int)
    cell_types_lt_high_percentile_by_bin = collections.defaultdict(int)
    lines_gte_high_percentile_by_bin = collections.defaultdict(int)
    lines_lt_high_percentile_by_bin = collections.defaultdict(int)
    cell_types_gte_low_percentile_by_bin = collections.defaultdict(int)
    cell_types_lt_low_percentile_by_bin = collections.defaultdict(int)
    lines_gte_low_percentile_by_bin = collections.defaultdict(int)
    lines_lt_low_percentile_by_bin = collections.defaultdict(int)
    lines = set()
    cell_types = set()
    max_bin = 0

    # note populating the dictionaries this way is much faster than looping through
    # each bin and finding the matching fraction rows
    for row in normalized_counts.where("chromosome == '%s'" % chromosome):
        bin_number = row["bin_number"]
        max_bin = max(max_bin, bin_number)
        percentile = row["percentile"]

        if row["file_key"] == 0:
            cell_types.add(row["cell_type"])
            summed_cell_type_percentiles_by_bin[bin_number] += percentile
            if percentile >= high:
                cell_types_gte_high_percentile_by_bin[bin_number] += 1
            else:
                cell_types_lt_high_percentile_by_bin[bin_number] += 1

            if percentile >= low:
                cell_types_gte_low_percentile_by_bin[bin_number] += 1
            else:
                cell_types_lt_low_percentile_by_bin[bin_number] += 1
        else:
            lines.add(row["file_key"])
            if percentile >= high:
                lines_gte_high_percentile_by_bin[bin_number] += 1
            else:
                lines_lt_high_percentile_by_bin[bin_number] += 1

            if percentile >= low:
                lines_gte_low_percentile_by_bin[bin_number] += 1
            else:
                lines_lt_low_percentile_by_bin[bin_number] += 1

    logging.debug(" - populating summary table with calculated summaries")

    for bin_number in xrange(max_bin+1):
        summary.row["bin_number"] = bin_number
        summary.row["chromosome"] = chromosome
        summary.row["avg_cell_type_percentile"] = summed_cell_type_percentiles_by_bin[bin_number] / len(cell_types)
        summary.row["cell_types_gte_95th_percentile"] = cell_types_gte_high_percentile_by_bin[bin_number]
        summary.row["cell_types_lt_95th_percentile"] = cell_types_lt_high_percentile_by_bin[bin_number]
        summary.row["lines_gte_95th_percentile"] = lines_gte_high_percentile_by_bin[bin_number]
        summary.row["lines_lt_95th_percentile"] = lines_lt_high_percentile_by_bin[bin_number]
        summary.row["cell_types_gte_5th_percentile"] = cell_types_gte_low_percentile_by_bin[bin_number]
        summary.row["cell_types_lt_5th_percentile"] = cell_types_lt_low_percentile_by_bin[bin_number]
        summary.row["lines_gte_5th_percentile"] = lines_gte_low_percentile_by_bin[bin_number]
        summary.row["lines_lt_5th_percentile"] = lines_lt_low_percentile_by_bin[bin_number]
        summary.row.append()
    summary.flush()

def normalize_plot_and_summarize(counts_file, output_directory, bin_size, skip_plot):
    delete_all_but_bin_counts_and_files_table(counts_file)

    # recreating the entirity of the remaining tables is quick and easier than updating prior records correctly

    counts = counts_file.root.bin_counts
    files = counts_file.root.files
    normalized_counts = create_normalized_counts_table(counts_file)
    summary = create_summary_table(counts_file)

    cell_types = all_cell_types(counts)
    chromosomes = all_chromosomes(counts)

    logging.info("Cell Types: %s", ", ".join(cell_types))

    for cell_type in cell_types:
        logging.info("Normalizing and calculating percentiles for cell type %s", cell_type)
        current_file_keys = file_keys(counts, cell_type)
        for file_key in current_file_keys:
           populate_normalized_counts(normalized_counts, counts, file_key, bin_size, files)
           populate_percentiles(normalized_counts, cell_type, file_key)
        populate_normalized_counts_for_cell_type(normalized_counts, cell_type, current_file_keys) 
        populate_percentiles(normalized_counts, cell_type)

    logging.info("Indexing normalized counts")
    normalized_counts.cols.bin_number.create_csindex()
    normalized_counts.cols.percentile.create_csindex()
    normalized_counts.cols.file_key.create_csindex()
    normalized_counts.cols.chromosome.create_csindex()

    if not skip_plot:
        if bp is None:
            logging.error('Skipping plotting because plots require a compatible version of bokeh -- '
                          'see https://github.com/BradnerLab/pipeline/wiki/bamliquidator#Install . %s'
                          % bokeh_import_error)
        else:
            logging.info("Plotting")
            for chromosome in chromosomes:
                plot(output_directory, normalized_counts, chromosome, cell_types)
            plot_summaries(output_directory, normalized_counts, chromosomes)

    logging.info("Summarizing")
    for chromosome in chromosomes:
        populate_summary(summary, normalized_counts, chromosome)
    summary.cols.avg_cell_type_percentile.create_csindex()

    # Iterating over this index in reverse order is hundreds of times slower than iterating
    # in ascending order in my tests, but copying into a reverse sorted table is very fast.
    # So we create a sorted summary table sorted in decreasing percentile order.  If we need to
    # iterate in the reverse sorted order, than this sorted_summary table should be used.
    # Otherwise, we should use the summary table (including the case of ascending percentile
    # order, which is fast since the table is indexed by that column). See
    # https://groups.google.com/d/topic/pytables-users/EKMUxghQiPQ/discussion
    sorted_summary = summary.copy(newname="sorted_summary", sortby=summary.cols.avg_cell_type_percentile,
                                  step=-1, checkCSI=True,
                                  title="Summary table sorted in decreasing percentile order")
    sorted_summary.cols.bin_number.create_csindex()

def debugging_handler(signal, frame):
    import pdb
    pdb.set_trace()

def main():
    parser = argparse.ArgumentParser(description='Calculate and plot normalized bin counts and percentiles. '
        'Normalized counts, percentiles, and summaries are stored in hdf5 tables in the file "normalized_counts.h5". '
        'Plots are stored in .html files.  The hdf5 and html files are stored by default in a new directory "output" '
        '(which can be overridden by argument, see below), and the program aborts if this directory already exists.') 
    parser.add_argument('-o', '--output_directory', default='output',
                        help='directory to create and output the h5 and/or html files to (aborts if already exists)')
    parser.add_argument('-b', '--bin_size', type=int, default=100000,
                        help="Number of base pairs in each bin -- should match the bin size in the bin_counts_h5_file")
    parser.add_argument('-v', '--validate', action='store_true',
                        help='validates the previously generated normalization and/or summary tables, returning '
                             'non-zero if any problems detected')
    parser.add_argument('-d', '--debug', action='store_true',
                        help='enables debugging hooks so ctr-c (SIGINT) enters debugging instead of halting execution')
    parser.add_argument('--skip_plot', action='store_true', help='skip generating plots (this can speed up execution')
    parser.add_argument('bin_counts_h5_file', help='the hdf5 file with a "counts" and "files" tables as generated by ' 
                                                   'bamliquidate_batch')
    args = parser.parse_args()

    if args.debug:
        import signal;
        signal.signal(signal.SIGINT, debugging_handler)

    if args.validate:
        sys.exit(validate(args.bin_counts_h5_file)) 

    os.mkdir(args.output_directory)

    counts_file = tables.open_file(args.bin_counts_h5_file, "r+")

    normalize_plot_and_summarize(counts_file, args.output_directory, args.bin_size, args.skip_plot)

    counts_file.close()

def validate(counts_file_path):
    counts_file = tables.open_file(counts_file_path, "r")

    error_count = 0

    counts = counts_file.root.bin_counts
    cell_types = all_cell_types(counts)
    num_cell_types = len(cell_types)
    num_files = 0
    for cell_type in cell_types:
        num_files += len(file_keys(counts, cell_type))

    logging.info("Verifying that summary files add up to %d and cell types add up to %d", num_files, num_cell_types)
    
    for row in counts_file.root.summary:
        if (num_cell_types != (row["cell_types_gte_95th_percentile"] + row["cell_types_lt_95th_percentile"])
         or num_cell_types != (row["cell_types_gte_5th_percentile"] + row["cell_types_lt_5th_percentile"])
         or num_files      != (row["lines_gte_95th_percentile"] + row["lines_lt_95th_percentile"])
         or num_files      != (row["lines_gte_5th_percentile"] + row["lines_lt_5th_percentile"])):
            error_count += 1
            logging.error("Summary row doesn't add up: %s", row[:])

    counts_file.close()

    if error_count != 0:
        logging.error("%d validation errors", error_count)
        return 1
    return 0
        
if __name__ == "__main__":
    main()

