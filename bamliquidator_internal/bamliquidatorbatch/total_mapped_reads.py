import argparse
import csv
import os
import subprocess

# returns (total mapped read count, list of pairs of chromosome to lengthes)
def total_mapped_reads(bam_file_path, chromosome_length_pairs=None, max_chromosome_name_length=None):
    chr_col         = 0
    length_col      = 1
    mapped_read_col = 2
    
    output = subprocess.check_output(["samtools", "idxstats", bam_file_path])
    # skip last two lines: the unmapped chromosome line and the empty line
    reader = csv.reader(output.split('\n')[:-2], delimiter='\t')
    file_name = os.path.basename(bam_file_path)
    file_count = 0
    
    chromosome_length_pairs = []
    for row in reader:
        chromosome = row[chr_col]

        if max_chromosome_name_length is not None and len(chromosome) >= max_chromosome_name_length:
            raise RuntimeError('Chromosome name "%s" exceeds the max supported chromosome name length (%d). '
                               'This max chromosome length may be updated in the code if necessary -- please '
                               'contact the bamliquidator developers for additional assistance.'
                               % (chromosome, max_chromosome_name_length))
                                
        file_count += int(row[mapped_read_col])
        chromosome_length_pairs.append((chromosome, int(row[length_col])))

    return (file_count, chromosome_length_pairs)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Counts total number of mapped reads by calling idxstats')
    parser.add_argument('bam_file_path')
    args = parser.parse_args()

    total_number_of_mapped_reads, _ = total_reads(args.bam_file_path)
    print total_number_of_mapped_reads

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
