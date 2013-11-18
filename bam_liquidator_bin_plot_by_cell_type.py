#!/usr/bin/env python

import bokeh.plotting as bp
import MySQLdb

def main():
    bp.output_file("output.html")

    print "Plotting summary"

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
    overall.canvas_width = 1300

    cursor.execute("SELECT distinct cell_type FROM chr1_bin_counts_by_cell_type");
    for row in cursor.fetchall():
        cell_type = row[0]
        print "Plotting " + cell_type

        bin_number = [] 
        count = [] 
        
        cell_type_cursor = db.cursor()
        cell_type_cursor.execute("SELECT bin, count FROM chr1_bin_counts_by_cell_type where cell_type = '%s'" % cell_type)
        for cell_type_row in cell_type_cursor.fetchall():
            bin_number.append(int(cell_type_row[0]))
            count.append(int(cell_type_row[1]))

        cell_type_plot = bp.scatter(bin_number, count)
        cell_type_plot.title = "%s counts per bin (chr1)" % cell_type 

    bp.save()
    bp.show()

if __name__ == "__main__":
    main()
