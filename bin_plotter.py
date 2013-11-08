#!/usr/bin/env python

import bokeh.plotting as bp
import MySQLdb

def main():
    bp.output_file("output.html")

    db = MySQLdb.connect(user="counter", db="meta_analysis")
    cursor = db.cursor()
    cursor.execute("SELECT bin, count FROM chr1_bin_counts");

    # todo: I should probably be able to do something like the following without copying: 
    #       bp.scatter(cursor.fetchall())

    bin_number = [] 
    count = [] 
    
    for row in cursor.fetchall():
        bin_number.append(int(row[0]))
        count.append(int(row[1]))

    overall = bp.scatter(bin_number, count)
    overall.title = "counts per bin accross all bam files (chr1)"

    cursor.execute("SELECT distinct parent FROM chr1_bin_counts_by_parent");
    for row in cursor.fetchall():
        parent = row[0]

        bin_number = [] 
        count = [] 
        
        parent_cursor = db.cursor()
        parent_cursor.execute("SELECT bin, count FROM chr1_bin_counts_by_parent where parent = '%s'" % parent)
        for parent_row in parent_cursor.fetchall():
            bin_number.append(int(parent_row[0]))
            count.append(int(parent_row[1]))

        parent_plot = bp.scatter(bin_number, count)
        parent_plot.title = "%s counts per bin (chr1)" % parent

    bp.save()
    bp.show()

if __name__ == "__main__":
    main()
