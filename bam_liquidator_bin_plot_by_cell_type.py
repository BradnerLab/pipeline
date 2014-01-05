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
counts = None
fractions = None
skip_population = False # useful if you just want to modify the csv and/or plotting 

def create_fractions_table(file_name):
    class BinCount(tables.IsDescription):
        bin_number = tables.UInt32Col(    pos=0);
        cell_type  = tables.StringCol(16, pos=1);
        chromosome = tables.StringCol(16, pos=2);
        count      = tables.Float64Col(   pos=3);
        percentile = tables.Float64Col(   pos=4);
        file_name  = tables.StringCol(64, pos=5);

    h5file = tables.open_file(file_name, mode = "w",
                              title = "normalized bam liquidator genome bin read counts")

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

def create_csv_file(chromosome, common_clause, cell_types):
    #csv_file = os.getcwd() + "/" + chromosome + ".csv"
    csv_file = "/tmp/" + chromosome + ".csv"
    
    def subqueries(percentile, by_file):
        label = "lines" if by_file else "cell_types"
        template = ("       (SELECT COUNT(*) FROM " + normalized_table(by_file) + " AS n2\n"
                    "         WHERE n2.bin = n1.bin "
                               "AND " + percentile_column(by_file) + " %s "
                               "AND" + common_clause + "\n"
                    "       ) AS " + label + "_with_%s_percentile")
        higher = template % (">= %d" % percentile, "bin_in_%dth_or_higher" % percentile) 
        lower  = template % ("< %d"  % percentile, "bin_lower_than_%dth"   % percentile)
        return higher + ",\n" + lower 

    sql = ("SELECT 'bin', 'average cell type percentile',\n"
                  "'cell types >= 95th percentile', 'cell types < 95th percentile',\n"
	          "'lines >= 95th percentile', 'lines < 95th percentile',\n" 
                  "'cell types >= 5th percentile', 'cell types < 5th percentile',\n"
	          "'lines >= 5th percentile', 'lines < 5th percentile'\n" 
           " UNION\n"
           "SELECT * FROM\n"
           "(SELECT n1.bin, AVG(n1.percentile_in_cell_type) AS average_cell_type_percentile,\n"
            + subqueries(95, by_file=False) + ",\n"
            + subqueries(95, by_file=True)  + ",\n"
            + subqueries(5,  by_file=False) + ",\n"
            + subqueries(5,  by_file=True)  + "\n"
           "   FROM normalized_bins AS n1\n"
           "  WHERE n1.chromosome = '" + chromosome + "' AND n1.counter_version = " + version + "\n"
           "  GROUP BY n1.bin\n"
           "  ORDER BY average_cell_type_percentile DESC, n1.bin ASC\n"
           "  ) AN_ALIAS \n"
           "  INTO OUTFILE '" + csv_file + "' FIELDS TERMINATED BY ',' LINES TERMINATED BY '\\n'")

    print sql

    cursor = db.cursor()
    cursor.execute(sql)

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

    counts = tables.open_file("bin_counts.h5", "r").root.counts
    if skip_population:
        fractions = tables.open_file("normalized_counts.h5", "r").root.counts
    else:
        fractions = create_fractions_table("normalized_counts.h5")

    cell_types = parse_args()

    print "cell_types = %s" % cell_types

    chromosomes = [ 'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr8', 'chr9', 'chrX',
                    'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17',
                    'chr18', 'chr19', 'chr20', 'chr21', 'chr22' ]

    for chromosome in chromosomes:
        print "Processing " + chromosome
        if not skip_population:
            populate_count_fractions(chromosome, cell_types)
            populate_count_percentiles_for_cell_types(chromosome, cell_types)
        #create_csv_file(chromosome, common_clause, cell_types)
        plot(chromosome, cell_types)

    print "Plotting summaries"
    plot_summaries(chromosomes)

if __name__ == "__main__":
    main()
