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

# globals
counts    = None
fractions = None
summary   = None
# todo: reverse these two booleans, since I always end up if not'ing them
skip_normalized_count_population = False # useful if you just want to modify the csv and/or plotting 
skip_plots = True # useful if you are just experimenting with summary tables

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

def all_cell_types():
    #print "Getting cell types"

    types = set()

    for row in counts:
        types.add(row["cell_type"])

    return types 

file_names_memo = {}
def file_names(cell_type):
    if not cell_type in file_names_memo:
        file_names = set() 
       
        #print "Getting file names for cell type " + cell_type
        for row in counts.where("cell_type == '%s'" % cell_type):
            file_names.add(row["file_name"])

        file_names_memo[cell_type] = file_names

        #print "memoizing files for " + cell_type + ": " + str(file_names_memo[cell_type]) 
        
    return file_names_memo[cell_type] 

def plot_summaries(chromosomes):
    bp.output_file("summary.html")
    
    for chromosome in chromosomes:
        plot_summary(chromosome)

    bp.save()

def plot_summary(chromosome):
    print " - plotting " + chromosome + " summary"

    condition = "(file_name == '*') & (chromosome == '%s')" % chromosome

    chromosome_count_by_bin = collections.defaultdict(int) 
    for row in fractions.where(condition):
        chromosome_count_by_bin[row["bin_number"]] += row["count"]
    
    overall = bp.scatter(chromosome_count_by_bin.keys(), chromosome_count_by_bin.values())
    overall.title = chromosome + " counts per bin across all bam files"

def plot(chromosome, cell_types):
    bp.output_file(chromosome + ".html")

    plot_summary(chromosome)

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

def total_count_for_file(file_name, chromosome):
    count = 0

    condition = "(file_name == '%s') & (chromosome == '%s')" % (file_name, chromosome)

    for row in counts.where(condition):
        count += row["count"]

    return count 

def populate_count_fractions(chromosome, cell_types):
    for cell_type in cell_types:
        files = file_names(cell_type)

        processed_a_single_file_for_this_cell_type = False
        cell_type_fractions = []

        for file_name in files:
            print " - populating count fractions for cell_type " + cell_type + " in file " + file_name 

            file_total_count = total_count_for_file(file_name, chromosome)

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

def populate_count_percentiles_for_cell_types(chromosome, cell_types):
    for cell_type in cell_types:
        print " - populating percentiles for cell_type " + cell_type
        populate_count_percentiles(chromosome, cell_type, file_name="*")

        for file_name in file_names(cell_type):
            print " - populating percentiles for file " + file_name + " in cell_type " + cell_type
            populate_count_percentiles(chromosome, cell_type, file_name)

def normalized_table(by_file):
    return "normalized_bins_by_file" if by_file else "normalized_bins"

def percentile_column(by_file):
    return "percentile_in_file" if by_file else "percentile_in_cell_type" 

def populate_count_percentiles(chromosome, cell_type, file_name):
    condition = "(chromosome == '%s') & (file_name == '%s') & (cell_type == '%s')" % (
        chromosome, file_name, cell_type)

    normalized_counts = fractions.read_where(condition, field='count')

    #print " -- condition = " + condition
    for row in fractions.where(condition):
        assert row["percentile"] == -1

        # todo: don't calculate percentile inside loop, use numpy.percentile with the axis argument 
        # or try http://grokbase.com/t/python/python-list/092235vj27/faster-scipy-percentileofscore
        row["percentile"] = stats.percentileofscore(normalized_counts, row["count"])
        row.update()

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

    if skip_normalized_count_population:
        # file wasn't created new, so already has a prior summary table that should be dropped
        h5file.remove_node("/", "summary")

    table = h5file.create_table("/", "summary", Summary, "bin count summary")

    table.flush()

    return table
    

def populate_summary(chromosome, cell_types):
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
                
def parse_args():
    parser = argparse.ArgumentParser(description='Calculate and plot normalized bin count percentiles. '
        'Normalized counts and percentiles are stored in the table normalized_bins. Plots are stored in '
        '.html files in the current directory.  Note that this is designed to be run once per cell type '
        'for a single version -- corresponding rows from normalized_bins should be deleted if new files '
        'are added or the same cell type needs to be re-processed for any reason (an error will occur '
        'and nothing will be changed in normalized_bins if there are already entries in this table for '
        'the given cell type and version)')
    parser.add_argument('cell_types', metavar='cell_type', nargs='+',
                        help='a cell type to calculate -- specify "all" to process all cell types')
    args = parser.parse_args()
    if len(args.cell_types) == 1 and args.cell_types[0] == 'all':
        return all_cell_types()
    else:
        return args.cell_types 

def main():
    global counts
    global fractions
    global summary

    counts = tables.open_file("bin_counts.h5", "r").root.counts

    if skip_normalized_count_population:
        normalized_counts_file = tables.open_file("normalized_counts.h5", "r+")
        fractions = normalized_counts_file.root.counts
    else:
        normalized_counts_file = tables.open_file("normalized_counts.h5", "w",
                                                  title = "normalized bam liquidator genome bin read counts")
        fractions = create_fractions_table(normalized_counts_file)

    summary = create_summary_table(normalized_counts_file)

    cell_types = parse_args()

    print "cell_types = %s" % cell_types

    chromosomes = [ 'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr8', 'chr9', 'chrX',
                    'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17',
                    'chr18', 'chr19', 'chr20', 'chr21', 'chr22' ]

    for chromosome in chromosomes:
        print "Processing " + chromosome
        if not skip_normalized_count_population:
            populate_count_fractions(chromosome, cell_types)
            populate_count_percentiles_for_cell_types(chromosome, cell_types)
        if not skip_plots: plot(chromosome, cell_types)

    if not skip_normalized_count_population:
        print "Indexing normalized counts"
        fractions.cols.bin_number.create_index()
        fractions.cols.percentile.create_index()
        fractions.cols.file_name.create_index()
        fractions.cols.chromosome.create_index()

    if not skip_plots:
        print "Plotting summaries"
        plot_summaries(chromosomes)

    for chromosome in chromosomes:
        print "Creating table summary for " + chromosome
        populate_summary(chromosome, cell_types)


if __name__ == "__main__":
    main()
