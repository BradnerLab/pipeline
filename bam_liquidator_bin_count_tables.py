#!/usr/bin/env python

import tables

class BinCount(tables.IsDescription):
    bin_number = tables.UInt32Col(    pos=0);
    cell_type  = tables.StringCol(16, pos=1);
    chromosome = tables.StringCol(16, pos=2);
    count      = tables.UInt64Col(    pos=3);
    file_name  = tables.StringCol(64, pos=4);

def create_table(file_name):
    # creates empty files, overwriting any prior existing files

    h5file = tables.open_file(file_name, mode = "w",
                            title = "bam liquidator genome bin read counts")

    table = h5file.create_table("/", "counts", BinCount, "bin counts")

    table.flush()

    return table

if __name__ == "__main__":
    create_table("bin_counts.h5")

'''
Before pytables, we prototyped with mysql:

CREATE DATABASE meta_analysis;
CREATE USER 'counter'@'localhost';

USE meta_analysis;

CREATE TABLE `counts` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `parent_directory` varchar(15) NOT NULL,
  `file_name` varchar(63) NOT NULL,
  `chromosome` varchar(15) NOT NULL,
  `bin` int(11) NOT NULL,
  `count` int(11) NOT NULL,
  `counter_version` int(11) NOT NULL,
  PRIMARY KEY (`id`),
  UNIQUE KEY `counter_version` (`counter_version`,`chromosome`,`bin`,`file_name`),
  KEY file (counter_version, file_name)
) ENGINE=MyISAM;

CREATE TABLE errors (id INT NOT NULL AUTO_INCREMENT, file_name VARCHAR(63) NOT NULL,
  chromosome VARCHAR(15) NOT NULL, bin INT NOT NULL, error VARCHAR(512) NOT NULL,
  counter_version INT NOT NULL, create_time TIMESTAMP DEFAULT CURRENT_TIMESTAMP, PRIMARY KEY (id));

CREATE TABLE run (id INT NOT NULL AUTO_INCREMENT, file_name VARCHAR(63) NOT NULL,
                         resume_count INT NOT NULL DEFAULT 0, start_time DATETIME NOT NULL,
                         end_time DATETIME, counter_version INT NOT NULL,
                         PRIMARY KEY (id), UNIQUE KEY (file_name, counter_version));

CREATE TABLE normalized_bins (
  id int NOT NULL AUTO_INCREMENT,
  cell_type varchar(15) NOT NULL,
  chromosome varchar(15) NOT NULL,
  bin int NOT NULL,
  count_fraction double NOT NULL,
  percentile_in_cell_type DOUBLE,
  counter_version int NOT NULL,
  PRIMARY KEY (id),
  UNIQUE KEY plotter (counter_version, chromosome, cell_type, bin),
  KEY overall_view (counter_version, chromosome, bin)
) ENGINE=MyISAM;

CREATE TABLE normalized_bins_by_file (
  id int NOT NULL AUTO_INCREMENT,
  cell_type varchar(15) NOT NULL,
  file_name varchar(63) NOT NULL,
  chromosome varchar(15) NOT NULL,
  bin int NOT NULL,
  count_fraction double NOT NULL,
  percentile_in_file DOUBLE,
  counter_version int NOT NULL,
  PRIMARY KEY (id),
  UNIQUE KEY plotter (counter_version, chromosome, cell_type, file_name, bin),
  KEY overall_view (counter_version, chromosome, bin)
) ENGINE=MyISAM;

GRANT SELECT ON meta_analysis.* TO 'counter'@'localhost';
GRANT SELECT,INSERT ON meta_analysis.counts TO 'counter'@'localhost';
GRANT SELECT,INSERT ON meta_analysis.errors TO 'counter'@'localhost';
GRANT SELECT,INSERT,UPDATE ON meta_analysis.run TO 'counter'@'localhost';
GRANT SELECT,INSERT,UPDATE ON meta_analysis.normalized_bins TO 'counter'@'localhost';
GRANT SELECT,INSERT,UPDATE ON meta_analysis.normalized_bins_by_file TO 'counter'@'localhost';
GRANT FILE ON *.* to 'counter'@'localhost';
'''

'''
   The MIT License (MIT) 

   Copyright (c) 2013 John DiMatteo 

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
