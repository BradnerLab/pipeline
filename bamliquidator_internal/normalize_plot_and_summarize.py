#!/usr/bin/env python

##################################################################################
# The MIT License (MIT)
#
# Copyright (c) 2013 John DiMatteo 
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

import sys
import os
import argparse
import bokeh.plotting as bp
import tables
import scipy.stats as stats
import collections

# note that my initial version didn't do any flush calls, which lead to bogus rows being added
# to the fractions table (which was evident when the normalized counts <= 95 + > 95 didn't add up right)
# -- I should probably look into why flush was necessary and/or file a bug with pytables

def create_fractions_table(h5file):
    class BinCount(tables.IsDescription):
        bin_number = tables.UInt32Col(    pos=0)
        cell_type  = tables.StringCol(16, pos=1)
        chromosome = tables.StringCol(16, pos=2)
        count      = tables.Float64Col(   pos=3)
        percentile = tables.Float64Col(   pos=4)
        file_name  = tables.StringCol(64, pos=5)

    table = h5file.create_table("/", "counts", BinCount, "bin counts")

    table.flush()

    return table

def all_cell_types(counts):
    types = set()

    for row in counts:
        types.add(row["cell_type"])

    return types 

file_names_memo = {}
def file_names(counts, cell_type):
    if not cell_type in file_names_memo:
        file_names = set() 
       
        #print "Getting file names for cell type " + cell_type
        for row in counts.where("cell_type == '%s'" % cell_type):
            file_names.add(row["file_name"])

        file_names_memo[cell_type] = file_names

        #print "memoizing files for " + cell_type + ": " + str(file_names_memo[cell_type]) 
        
    return file_names_memo[cell_type] 

def plot_summaries(output_directory, fractions, chromosomes):
    bp.output_file(output_directory + "/summary.html")
    
    for chromosome in chromosomes:
        plot_summary(fractions, chromosome)

    bp.save()

def plot_summary(fractions, chromosome):
    print " - plotting " + chromosome + " summary"

    condition = "(file_name == '*') & (chromosome == '%s')" % chromosome

    chromosome_count_by_bin = collections.defaultdict(int) 
    for row in fractions.where(condition):
        chromosome_count_by_bin[row["bin_number"]] += row["count"]
    
    overall = bp.scatter(chromosome_count_by_bin.keys(), chromosome_count_by_bin.values())
    overall.title = chromosome + " counts per bin across all bam files"

def plot(output_directory, fractions, chromosome, cell_types):
    bp.output_file(output_directory + "/" + chromosome + ".html")

    plot_summary(fractions, chromosome)

    for cell_type in cell_types:
        print " - plotting " + cell_type

        bin_number = [] 
        count = [] 
        
        condition = "(file_name == '*') & (chromosome == '%s') & (cell_type == '%s')" % (chromosome, cell_type)

        for row in fractions.where(condition):
            bin_number.append(row["bin_number"])
            count.append(row["count"])

        cell_type_plot = bp.scatter(bin_number, count)
        cell_type_plot.title = "%s counts per bin" % cell_type 

    bp.save()

def total_count_for_file(counts, file_name, chromosome):
    count = 0

    condition = "(file_name == '%s') & (chromosome == '%s')" % (file_name, chromosome)

    for row in counts.where(condition):
        count += row["count"]

    return count 

def populate_count_fractions(fractions, counts, chromosome, cell_types):
    for cell_type in cell_types:
        files = file_names(counts, cell_type)

        processed_a_single_file_for_this_cell_type = False
        cell_type_fractions = []

        for file_name in files:
            print " - populating count fractions for cell_type " + cell_type + " in file " + file_name 

            file_total_count = total_count_for_file(counts, file_name, chromosome)

            #print "There are %d reads in this file for this chromosome" % file_total_count 

            file_divisor = float(file_total_count)
            cell_type_divisor = float(file_total_count * len(files))

            condition = "(file_name == '%s') & (chromosome == '%s')" % (file_name, chromosome)

            for bin_number, count_row in enumerate(counts.where(condition)):
                assert bin_number == count_row["bin_number"]

                file_count_fraction = count_row["count"] / file_divisor
                cell_type_count_fraction = count_row["count"] / cell_type_divisor 

                fractions.row["bin_number"] = bin_number 
                fractions.row["cell_type"] = cell_type 
                fractions.row["chromosome"] = chromosome 
                fractions.row["file_name"] = file_name
                fractions.row["count"] = file_count_fraction 
                fractions.row["percentile"] = -1.0
                fractions.row.append()

                if not processed_a_single_file_for_this_cell_type:
                    cell_type_fractions.append(cell_type_count_fraction)
                else:
                    cell_type_fractions[bin_number] += cell_type_count_fraction

            fractions.flush()
            processed_a_single_file_for_this_cell_type = True

        for bin_number, cell_type_fraction in enumerate(cell_type_fractions):
            fractions.row["bin_number"] = bin_number
            fractions.row["cell_type"] = cell_type 
            fractions.row["chromosome"] = chromosome 
            fractions.row["file_name"] = "*" 
            # using "*" instead of "" due to pytables empty string query bug
            # -- https://github.com/PyTables/PyTables/issues/184
            fractions.row["count"] = cell_type_fraction 
            fractions.row["percentile"] = -1.0
            fractions.row.append()

        fractions.flush()

def populate_count_percentiles_for_cell_types(fractions, counts, chromosome, cell_types):
    for cell_type in cell_types:
        print " - populating percentiles for cell_type " + cell_type
        populate_count_percentiles(fractions, counts, chromosome, cell_type, file_name="*")

        for file_name in file_names(counts, cell_type):
            print " - populating percentiles for file " + file_name + " in cell_type " + cell_type
            populate_count_percentiles(fractions, counts, chromosome, cell_type, file_name)

def populate_count_percentiles(fractions, counts, chromosome, cell_type, file_name):
    condition = "(chromosome == '%s') & (file_name == '%s') & (cell_type == '%s')" % (
        chromosome, file_name, cell_type)

    bin_numbers = []
    normalized_counts = []

    for row in fractions.where(condition):
        bin_numbers.append(row["bin_number"])
        normalized_counts.append(row["count"])

    percentiles = (stats.rankdata(normalized_counts) - 1) / (len(normalized_counts)-1) * 100
    # percentiles calculated in bulk as suggested at 
    # http://grokbase.com/t/python/python-list/092235vj27/faster-scipy-percentileofscore

    for i, row in enumerate(fractions.where(condition)):
        assert bin_numbers[i] == row["bin_number"]
        row["percentile"] = percentiles[i]
        row.update()
    fractions.flush()

def create_summary_table(h5file):
    class Summary(tables.IsDescription):
        bin_number = tables.UInt32Col(                    pos=0)
        chromosome = tables.StringCol(16,                 pos=2)
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
    

def populate_summary(summary, fractions, chromosome, cell_types):
    #print " - calculating summaries"

    condition = "chromosome == '%s'" % chromosome

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
    max_bin = 0

    # note populating the dictionaries this way is much faster than looping through
    # each bin and finding the matching fraction rows
    for row in fractions.where(condition):
        bin_number = row["bin_number"]
        max_bin = max(max_bin, bin_number)
        percentile = row["percentile"]

        if row["file_name"] == "*":
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
            lines.add(row["file_name"])
            if percentile >= high:
                lines_gte_high_percentile_by_bin[bin_number] += 1
            else:
                lines_lt_high_percentile_by_bin[bin_number] += 1

            if percentile >= low:
                lines_gte_low_percentile_by_bin[bin_number] += 1
            else:
                lines_lt_low_percentile_by_bin[bin_number] += 1

    #print " - populating summary table with calculated summaries"

    for bin_number in xrange(max_bin):
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


def normalize_plot_and_summarize(counts, cell_types, output_directory):
    os.mkdir(output_directory)

    skip_plots = False # useful if you are just experimenting with normalization and/or summary tables

    normalized_counts_file = tables.open_file(output_directory + "/normalized_counts.h5", "w",
                                              title = "normalized bam liquidator genome bin read counts")

    fractions = create_fractions_table(normalized_counts_file)
    summary = create_summary_table(normalized_counts_file)

    print "cell_types = %s" % cell_types

    chromosomes = [ 'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr8', 'chr9', 'chrX',
                    'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17',
                    'chr18', 'chr19', 'chr20', 'chr21', 'chr22' ]

    for chromosome in chromosomes:
        print "Processing " + chromosome
        populate_count_fractions(fractions, counts, chromosome, cell_types)
        populate_count_percentiles_for_cell_types(fractions, counts, chromosome, cell_types)
        if not skip_plots: plot(output_directory, fractions, chromosome, cell_types)

    print "Indexing normalized counts"
    fractions.cols.bin_number.create_index()
    fractions.cols.percentile.create_index()
    fractions.cols.file_name.create_index()
    fractions.cols.chromosome.create_index()

    if not skip_plots:
        print "Plotting summaries"
        plot_summaries(output_directory, fractions, chromosomes)

    for chromosome in chromosomes:
        print "Creating table summary for " + chromosome
        populate_summary(summary, fractions, chromosome, cell_types)

    normalized_counts_file.close() 
                
def main():
    parser = argparse.ArgumentParser(description='Calculate and plot normalized bin counts and percentiles. '
        'Normalized counts, percentiles, and summaries are stored in hdf5 tables in the file "normalized_counts.h5". '
        'Plots are stored in .html files.  The hdf5 and html files are stored by default in a new directory "output" '
        '(which can be overridden by argument, see below), and the program aborts if this directory already exists.') 
    parser.add_argument('--cell_type', nargs='*',
                        help='restrict to one or more cell types to process (defaults to processing all cell types)')
    parser.add_argument('--output_directory', default='output',
                        help='directory to create and output the h5 and/or html files to (aborts if already exists)')
    parser.add_argument('bin_counts_h5_file', help='the hdf5 file with a "counts" table with columns ' 
                                                   '"bin_number", "cell_type", "chromosome", "count", and '
                                                   '"file_name", e.g. "bin_counts.hdf5" as generated by '
                                                   'bamliquidate_batch')
    args = parser.parse_args()

    counts_file = tables.open_file(args.bin_counts_h5_file, "r")
    counts = counts_file.root.counts

    if args.cell_type is None:
        args.cell_types = all_cell_types(counts)

    normalize_plot_and_summarize(counts, args.cell_types, args.output_directory)
    counts_file.close()

if __name__ == "__main__":
    main()
