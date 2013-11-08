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
    
    for row in cursor.fetchall() :
        bin_number.append(int(row[0]))
        count.append(int(row[1]))

    bp.scatter(bin_number, count)
    bp.show()

if __name__ == "__main__":
    main()
