#!/usr/bin/env python

import numpy as np

import matplotlib
#from matplotlib.backends.backend_pgf import FigureCanvasPgf
import matplotlib.pyplot as plt

import MySQLdb
import math            

db = MySQLdb.connect(user="counter", db="meta_analysis")
counter_version = " counter_version = 203 "

def main():
    #matplotlib.backend_bases.register_backend('pdf', FigureCanvasPgf)
    cursor = db.cursor()

    cursor.execute("select distinct file_name from run where " + counter_version) # + " limit 25")
    file_names = []
    for row in cursor.fetchall():
        file_names.append(str(row[0]))

    num_subplot_rows = int(math.ceil(math.sqrt(len(file_names))))
    num_subplot_cols = int(num_subplot_rows)
    print num_subplot_rows, num_subplot_cols
    subplot_row = 0
    subplot_col = 0
    print "setting up count subplots"
    count_fig, count_axes = plt.subplots(num_subplot_rows, num_subplot_cols, sharex=True, sharey=True)
    print "setting up histogram subplots"
    hist_fig, hist_axes = plt.subplots(num_subplot_rows, num_subplot_cols, sharex=True, sharey=True)

    for file_name in file_names:
        print "Populating plot %d %d" % (subplot_row, subplot_col)
        chromosome = 'chr1' # todo, loop over chromosomes and plot them all, maybe in a single subplot?

        cursor.execute("select bin, count from counts "
                        "where " + counter_version + 
                          "and file_name = '" + file_name + "' "
                          "and chromosome = '" + chromosome + "' " 
                        "order by bin")

        # todo: I should probably be able to do something like the following without copying: 
        #       bp.scatter(cursor.fetchall())
        # -- see http://stackoverflow.com/questions/7061824/whats-the-most-efficient-way-to-covert-mysql-output-into-a-numpy-array-in-pytho

        bin_number = [] 
        count = []
        
        for row in cursor.fetchall():
            bin_number.append(int(row[0]))
            count.append(int(row[1]))

        # todo: use numpy? bin_counts = np.array([bin_number, count], int)

        #count_plot = counts.add_subplot(num_subplots, 1, subplot)
        count_plot = count_axes[subplot_row, subplot_col]
        count_plot.scatter(bin_number, count, s=1, lw = 0)
        count_plot.axis('off')
        #count_plot.ylabel(file_name + ' - count')
        #count_plot.xlabel('bin #')

        #histogram_plot = histograms.add_subplot(num_subplots, 1, subplot)
        #print subplot
        #histogram_plot.hist(count, 256)
        hist_plot = hist_axes[subplot_row, subplot_col]
        hist_plot.hist(count, 256)
        hist_plot.axis('off')
        #histogram_plot.xlabel('count histogram')

        subplot_col += 1
        if subplot_col >= num_subplot_cols:
            subplot_col= 0
            subplot_row += 1
    
    #plt.savefig('./figure.pdf')
    #plt.show()
    print "writing count plots to disk" 
    count_fig.subplots_adjust(wspace=0, hspace=0)
    count_fig.set_size_inches(num_subplot_rows,num_subplot_cols) # inch per plot
    count_fig.savefig('counts.png', dpi=900)
    
    print "writing histogram plots to disk" 
    hist_fig.subplots_adjust(wspace=0, hspace=0)
    hist_fig.set_size_inches(num_subplot_rows,num_subplot_cols) # inch per plot
    hist_fig.savefig('histograms.png', dpi=900)

if __name__ == "__main__":
    main()
