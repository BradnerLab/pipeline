#!/usr/bin/env python

import normalize_plot_and_summarize as nps 
from flattener import write_tab_for_all 

import argparse
import csv
import datetime
import errno
import os
import subprocess
import tables
import logging
import sys
import abc
import collections
import numpy

from time import time 
from os.path import basename
from os.path import dirname

__version__ = '1.3.8'

default_black_list = ["chrUn", "_random", "Zv9_", "_hap"]

def create_files_table(h5file):
    class Files(tables.IsDescription):
        key       = tables.UInt32Col(    pos=0) # is there an easier way to assign keys?
        length    = tables.UInt64Col(    pos=2)
        # file_name would be included here, but pytables doesn't support variable length strings as table column
        # so it is instead in a vlarray "file_names" 

    table = h5file.create_table("/", "files", Files, "File keys and reference sequence lengths corresponding "
                                                     "to the counts table")
    table.flush()

    return table

def create_file_names_array(h5file):
    # vlarray of strings only supports a single column, so the file_key is implicitly the array index
    array = h5file.create_vlarray("/", "file_names", tables.VLStringAtom(),
                                "File names with index corresponding to Files table key")
    array.append("*") # index/key 0 is reserved for this
    array.flush()

    return array

def all_bam_file_paths_in_directory(bam_directory):
    bam_file_paths = []
    for dirpath, _, files in os.walk(bam_directory, followlinks=True):
        for file_ in files:
            if file_.endswith(".bam"):
                bam_file_paths.append(os.path.join(dirpath, file_))
    return bam_file_paths

def bam_file_paths_with_no_file_entries(file_names, bam_file_paths):
    with_no_counts = []

    for bam_file_path in bam_file_paths:
        if basename(bam_file_path) not in file_names:
            with_no_counts.append(bam_file_path)

    return with_no_counts

# BaseLiquidator is an abstract base class, with concrete classes BinLiquidator and RegionLiquidator
# that implement the abstract methods.
class BaseLiquidator(object):
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def liquidate(self, bam_file_path, extension, sense = None):
        pass

    @abc.abstractmethod
    def normalize(self):
        pass

    @abc.abstractmethod
    def create_counts_table(self, h5file):
        pass

    def __init__(self, executable, counts_table_name, output_directory, bam_file_path,
                 include_cpp_warnings_in_stderr = True, counts_file_path = None, number_of_threads = 0):
        # clear all memoized values from any prior runs
        nps.file_keys_memo = {}

        self.timings = collections.OrderedDict()

        self.output_directory = output_directory
        self.counts_file_path = counts_file_path
        self.include_cpp_warnings_in_stderr = include_cpp_warnings_in_stderr
        self.number_of_threads = number_of_threads
        self.chromosome_patterns_to_skip = [] 

        # This script may be run by either a developer install from a git pipeline checkout,
        # or from a user install so that the exectuable is on the path.  First we try to
        # find the exectuable for a developer install, and if that fails we look on the
        # standard path.
        if basename(dirname(dirname(os.path.realpath(__file__)))) == 'bamliquidator_internal':
            # look for developer executable location 
            self.executable_path = os.path.join(dirname(dirname(os.path.realpath(__file__))), executable)
            if not os.path.isfile(self.executable_path):
                exit("%s is missing -- try cd'ing into the directory and running 'make'" % self.executable_path)
        else:
            # just look on standard path
            self.executable_path = executable 

        mkdir_if_not_exists(output_directory)

        if self.counts_file_path is None:
            self.counts_file_path = os.path.join(output_directory, "counts.h5")
        
            counts_file = tables.open_file(self.counts_file_path, mode = "w",
                                           title = 'bam liquidator genome read counts - version %s' % __version__)
        else:
            counts_file = tables.open_file(self.counts_file_path, "r+")

        try: 
            counts = counts_file.get_node("/", counts_table_name)
            files = counts_file.root.files
            file_names = counts_file.root.file_names
        except:
            counts = self.create_counts_table(counts_file)
            files = create_files_table(counts_file)
            file_names = create_file_names_array(counts_file)

        if os.path.isdir(bam_file_path):
            self.bam_file_paths = all_bam_file_paths_in_directory(bam_file_path)
        else:
            self.bam_file_paths = [bam_file_path]
       
        self.bam_file_paths = bam_file_paths_with_no_file_entries(file_names, self.bam_file_paths)

        self.preprocess(files, file_names)

        counts_file.close() # bamliquidator_bins/bamliquidator_regions will open this file and modify
                            # it, so it is probably best that we not hold an out of sync reference

    
    # adds files being liquidated to the files table and populates the following member dictionaries:
    # 1) file_name -> [(chromosome, sequence length), ...] 
    # 2) file_name -> total mapped count
    # 3) file_name -> file key number
    def preprocess(self, files, file_names):
        self.file_to_chromosome_length_pairs = {}
        self.file_to_count = {}
        self.file_to_key = {}

        chr_col         = 0
        length_col      = 1
        mapped_read_col = 2

        # bam file keys start at 1.
        # key 0 is special and denotes "no specific file", which
        # is used in normalizated_counts tables to mean an average or total for all bam files
        # of a specific cell type.
        next_file_key = 0 # see += 1 below
        for file_record in files:
            next_file_key = max(next_file_key, file_record["key"])
        next_file_key += 1

        for bam_file_path in self.bam_file_paths:
            args = ["samtools", "idxstats", bam_file_path]
            output = subprocess.check_output(args)
            # skip last two lines: the unmapped chromosome line and the empty line
            reader = csv.reader(output.split('\n')[:-2], delimiter='\t')
            file_name = basename(bam_file_path)
            file_count = 0
            
            chromosome_length_pairs = []
            for row in reader:
                chromosome = row[chr_col]
                if len(chromosome) >= nps.chromosome_name_length:
                    raise RuntimeError('Chromosome name "%s" exceeds the max supported chromosome name length (%d). '
                                       'This max chromosome length may be updated in the code if necessary -- please '
                                       'contact the bamliquidator developers for additional assistance.'
                                       % (chromosome, nps.chromosome_name_length))
                                        
                file_count += int(row[mapped_read_col])
                chromosome_length_pairs.append((chromosome, int(row[length_col])))
            
            files.row["key"] = next_file_key
            files.row["length"] = file_count
            files.row.append()
            file_names.append(file_name)

            self.file_to_chromosome_length_pairs[file_name] = chromosome_length_pairs
            self.file_to_count[file_name] = file_count
            self.file_to_key[file_name] = next_file_key

            next_file_key += 1

        files.flush()
        file_names.flush()
        assert(len(file_names) - 1 == len(files))
        assert(len(file_names) == next_file_key)

    def batch(self, extension, sense):
        for i, bam_file_path in enumerate(self.bam_file_paths):
            logging.info("Liquidating %s (file %d of %d)", bam_file_path, i+1, len(self.bam_file_paths))

            return_code = self.liquidate(bam_file_path, extension, sense)
            if return_code != 0:
                raise Exception("%s failed with exit code %d" % (self.executable_path, return_code))

        start = time()
        self.normalize()
        duration = time() - start
        logging.info("Post liquidation processing took %f seconds", duration)
        self.log_time('post_liquidation', duration)

    def flatten(self):
        logging.info("Flattening HDF5 tables into text files")
        start = time()

        with tables.open_file(self.counts_file_path, mode = "r") as counts_file:
            write_tab_for_all(counts_file, self.output_directory)

        duration = time() - start
        logging.info("Flattening took %f seconds" % duration)
        self.log_time('flattening', duration)

    def chromosome_args(self, bam_file_name, skip_non_canonical):
        args = []
        for chromosome, length in self.file_to_chromosome_length_pairs[bam_file_name]:
            if skip_non_canonical:
                if any(pattern in chromosome for pattern in self.chromosome_patterns_to_skip):
                    continue
            args.append(chromosome)
            args.append(str(length))
        return args
        
    def logging_cpp_args(self):
        return [os.path.join(self.output_directory, "log.txt"), "1" if self.include_cpp_warnings_in_stderr else "0"]

    def log_time(self, title, seconds):
        self.timings[title] = seconds

    def write_timings_to_junit_xml(self):
        with open(os.path.join(self.output_directory, 'timings.xml'), 'w') as xml:
            xml.write('<testsuite tests="%d">\n' % len(self.timings.keys()))
            for title in self.timings:
                xml.write('\t<testcase classname="bamliquidator" name="%s" time="%f"/>\n' % (title, self.timings[title]))
            xml.write('</testsuite>\n')

class BinLiquidator(BaseLiquidator):
    def __init__(self, bin_size, output_directory, bam_file_path,
                 counts_file_path = None, extension = 0, sense = '.', skip_plot = False,
                 include_cpp_warnings_in_stderr = True, number_of_threads = 0, blacklist = default_black_list):
        self.bin_size = bin_size
        self.skip_plot = skip_plot
        super(BinLiquidator, self).__init__("bamliquidator_bins", "bin_counts", output_directory, bam_file_path,
                                            include_cpp_warnings_in_stderr, counts_file_path, number_of_threads)
        self.chromosome_patterns_to_skip = blacklist
        self.batch(extension, sense)

    def liquidate(self, bam_file_path, extension, sense = None):
        if sense is None: sense = '.'

        cell_type = basename(dirname(bam_file_path))
        if cell_type == '':
            cell_type = '-'
        bam_file_name = basename(bam_file_path)
        args = [self.executable_path, str(self.number_of_threads), cell_type, str(self.bin_size), str(extension), sense, bam_file_path, 
                str(self.file_to_key[bam_file_name]), self.counts_file_path]
        args.extend(self.logging_cpp_args())
        args.extend(self.chromosome_args(bam_file_name, skip_non_canonical=True))

        start = time()
        return_code = subprocess.call(args)
        duration = time() - start

        reads = self.file_to_count[bam_file_name]
        rate = reads / (10**6) / duration
        logging.info("Liquidation completed: %f seconds, %d reads, %f millions of reads per second", duration, reads, rate)
        self.log_time('liquidation', duration)

        return return_code
       
    def normalize(self):
        with tables.open_file(self.counts_file_path, mode = "r+") as counts_file:
            nps.normalize_plot_and_summarize(counts_file, self.output_directory, self.bin_size, self.skip_plot) 

    def create_counts_table(self, h5file):
        class BinCount(tables.IsDescription):
            bin_number = tables.UInt32Col(    pos=0)
            cell_type  = tables.StringCol(16, pos=1)
            chromosome = tables.StringCol(nps.chromosome_name_length, pos=2)
            count      = tables.UInt64Col(    pos=3)
            file_key   = tables.UInt32Col(    pos=4)

        table = h5file.create_table("/", "bin_counts", BinCount, "bin counts")
        table.flush()
        return table

class RegionLiquidator(BaseLiquidator):
    def __init__(self, regions_file, output_directory, bam_file_path,
                 region_format=None, counts_file_path = None, extension = 0, sense = '.',
                 include_cpp_warnings_in_stderr = True, number_of_threads = 0):
        self.regions_file = regions_file
        self.region_format = region_format
        if self.region_format is None:
            _, self.region_format = os.path.splitext(regions_file)
            if len(self.region_format) > 0 and self.region_format[0] == '.':
                self.region_format = self.region_format[1:]
        if self.region_format not in ("gff", "bed"):
            raise RuntimeError("Only bed and gff region file formats are supported -- %s format specified"
                               % str(self.region_format))

        super(RegionLiquidator, self).__init__("bamliquidator_regions", "region_counts", output_directory, 
                                               bam_file_path, include_cpp_warnings_in_stderr, counts_file_path, number_of_threads)
        
        self.batch(extension, sense)

    def liquidate(self, bam_file_path, extension, sense = None):
        bam_file_name = basename(bam_file_path)
        args = [self.executable_path, str(self.number_of_threads), self.regions_file, str(self.region_format), str(extension), bam_file_path, 
                str(self.file_to_key[bam_file_name]), self.counts_file_path]
        args.extend(self.logging_cpp_args())
        if sense is None:
            args.append('_') # _ means use strand specified in region file (or . if none specified)
        else:
            args.append(sense)
        args.extend(self.chromosome_args(bam_file_name, skip_non_canonical=False))

        start = time()
        return_code = subprocess.call(args)
        duration = time() - start

        logging.info("Liquidation completed: %f seconds", duration)
        self.log_time('liquidation', duration)

        return return_code

    def normalize(self):
        with tables.open_file(self.counts_file_path, mode = "r+") as counts_file:
            nps.normalize_regions(counts_file.root.region_counts, counts_file.root.files)

    def create_counts_table(self, h5file):
        class Region(tables.IsDescription):
            file_key         = tables.UInt32Col(    pos=0)
            chromosome       = tables.StringCol(nps.chromosome_name_length, pos=1)
            region_name      = tables.StringCol(64, pos=2)
            start            = tables.UInt64Col(    pos=3)
            stop             = tables.UInt64Col(    pos=4)
            strand           = tables.StringCol(1,  pos=5)
            count            = tables.UInt64Col(    pos=6)
            normalized_count = tables.Float64Col(   pos=7)

        table = h5file.create_table("/", "region_counts", Region, "region counts")
        table.flush()
        return table

def write_bamToGff_matrix(output_file_path, h5_region_counts_file_path):
    with tables.open_file(h5_region_counts_file_path, "r") as counts_file:
        with open(output_file_path, "w") as output:
            file_keys = []

            output.write("GENE_ID\tlocusLine")
            for file_record in counts_file.root.files:
                file_key = file_record["key"] 
                file_keys.append(file_key)
                output.write("\tbin_1_%s" % counts_file.root.file_names[file_key])
            output.write("\n")

            number_of_files = len(file_keys)
            number_of_regions = counts_file.root.region_counts.nrows / number_of_files 

            # first loop through all but the last file index, storing those counts 
            prior_region_counts = numpy.zeros((number_of_regions,  number_of_files - 1))
            for col, file_key in enumerate(file_keys[:-1]):
                for row, region in enumerate(counts_file.root.region_counts.where("file_key == %d" % file_key)):
                    prior_region_counts[row, col] = region["normalized_count"]

            # then loop through the last index,
            # printing the region columns and the counts for the prior files,
            # along with the count for the last index
            for row, region in enumerate(counts_file.root.region_counts.where("file_key == %d" % file_keys[-1])):
                output.write("%s\t%s(%s):%d-%d" % (region["region_name"], region["chromosome"],
                    region["strand"], region["start"], region["stop"]))
                for col in range(0, number_of_files-1):
                    output.write("\t%s" % round(prior_region_counts[row, col], 4))
                output.write("\t%s\n" % round(region["normalized_count"], 4))

def configure_logging(args):
    # Using root logger so we can just do logging.info/warn/error in this and other files.
    # If people start using bamliquidator_batch as an imported module, then we should probably
    # change this logging to not use the root logger directly.
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    file_handler = logging.FileHandler(os.path.join(args.output_directory, 'log.txt'))
    file_handler.setFormatter(logging.Formatter('%(asctime)s %(levelname)s\t%(message)s',
                                                datefmt='%Y-%m-%d %H:%M:%S'))
    
    logger.addHandler(file_handler)
    # todo: add bamliquidator version to the starting log message
    logging.info("Starting %s %s with args %s", basename(sys.argv[0]), __version__, vars(args))

    # Adding console handler after writing the startup log entry.  The startup log could be useful 
    # in a file that is being appended to from a prior run, but would be annonying on stderr.

    console_handler = logging.StreamHandler()
    if args.quiet:
        console_handler.setLevel(logging.ERROR)
    else:
        console_handler.setLevel(logging.INFO)

    class FormatterNotFormattingInfo(logging.Formatter):
        def __init__(self, fmt):
            logging.Formatter.__init__(self, fmt)

        def format(self, record):
            if record.levelno == logging.INFO:
                return record.getMessage()
            return logging.Formatter.format(self, record)

    console_handler.setFormatter(FormatterNotFormattingInfo('%(levelname)s\t%(message)s'))
    logger.addHandler(console_handler)

def mkdir_if_not_exists(directory):
    try:
        os.mkdir(directory)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

def main():
    parser = argparse.ArgumentParser(description='Count the number of base pair reads in each bin or region '
                                                 'in the bam file(s) at the given directory, and then normalize, plot bins, '
                                                 'and summarize the counts in the output directory.  For additional '
                                                 'help, please see https://github.com/BradnerLab/pipeline/wiki')

    mut_exclusive_group = parser.add_mutually_exclusive_group()
    mut_exclusive_group.add_argument('-b', '--bin_size', type=int, default=100000,
                        help="Number of base pairs in each bin -- the smaller the bin size the longer the runtime and "
                             "the larger the data files (default is 100000)")
    mut_exclusive_group.add_argument('-r', '--regions_file',
                        help='a region file in either .gff or .bed format')

    parser.add_argument('-o', '--output_directory', default='output',
                        help='Directory to output the h5, log, gff, tab, and/or html files to.  Creates directory if necessary.  '
                             'May overwrite prior run results if present. Default is "./output".')
    parser.add_argument('-c', '--counts_file', default=None,
                        help='HDF5 counts file from a prior run to be appended to.  If unspecified, defaults to '
                             'creating a new file "counts.h5" in the output directory.')
    parser.add_argument('-f', '--flatten', action='store_true',
                        help='flatten all HDF5 tables into tab delimited text files in the output directory, one for each '
                              'chromosome (note that HDF5 files can be efficiently queried and used directly -- e.g. please '
                              'see http://www.pytables.org/ for easy to use Python APIs and '
                              'http://www.hdfgroup.org/products/java/hdf-java-html/hdfview/ for an easy to use GUI for '
                              'browsing HDF5 files)')
    parser.add_argument('-e', '--extension', type=int, default=0,
                        help='Extends reads by n bp (default is 0)')
    parser.add_argument('--sense', default=None, choices=['+', '-', '.'],
                        help="Map to '+' (forward), '-' (reverse) or '.' (both) strands. For gff regions, default is to use "
                             "the sense specified by the gff file; otherwise, default maps to both.")
    parser.add_argument('-m', '--match_bamToGFF', default=False, action='store_true',
                        help="match bamToGFF_turbo.py matrix output format, storing the result as matrix.txt in the output folder")
    parser.add_argument('--region_format', default=None, choices=['gff', 'bed'],
                        help="Interpret region file as having the given format.  Default is to deduce format from file extension.")
    parser.add_argument('--skip_plot', action='store_true', help='Skip generating plots.  (This can speed up execution.)')
    parser.add_argument('--black_list', nargs='+', type=str, default=default_black_list,
                        help='One or more (space separated) chromosome patterns to skip during bin liquidation. Default is '
                             'to skip any chromosomes that contain any of the following substrings: %s. ' %  " ".join(default_black_list))
    parser.add_argument('-q', '--quiet', action='store_true',
                        help='Informational and warning output is suppressed so only errors are written to the console (stderr).  '
                             'All bamliquidator logs are still written to log.txt in the output directory.  This also disables '
                             'samtools error messages to stderr, but a corresponding bamliquidator message should still be logged '
                             'in log.txt.')
    parser.add_argument('-n', '--number_of_threads', type=int, default=0,
                        help='Number of threads to run concurrently during liquidation.  Defaults to the total number of logical '
                             'cpus on the system.')
    parser.add_argument('--xml_timings', action='store_true',
                        help='Write performance timings to junit style timings.xml in output folder, which is useful for '
                             'tracking performance over time with automatically generated Jenkins graphs')
    parser.add_argument('--version', action='version', version='%s %s' % (basename(sys.argv[0]), __version__))
    parser.add_argument('bam_file_path', 
                        help='The directory to recursively search for .bam files for counting.  Every .bam file must '
                             'have a corresponding .bai file at the same location.  To count just a single file, '
                             'provide the .bam file path instead of a directory.  The parent directory (up to 16 char) of each '
                             '.bam file is interpreted as the cell type (e.g. mm1s might be an appropriate directory '
                             'name).  Bam files in the same directory are grouped together for plotting. Plots use '
                             'normalized counts, such that all .bam files in the same directory have bin '
                             'counts that add up to 1 for each chromosome.  If your .bam files are not in this '
                             'directory format, please consider creating a directory of sym links to your actual '
                             '.bam and .bai files. If the .bam file already has 1 or more reads in the HDF5 counts file, '
                             'then that .bam file is skipped from liquidation, but is still included in normalization, '
                             'plotting, and summaries.')

    args = parser.parse_args()

    assert(tables.__version__ >= '3.0.0')

    mkdir_if_not_exists(args.output_directory)

    configure_logging(args)

    if args.regions_file is None:
        liquidator = BinLiquidator(args.bin_size, args.output_directory, args.bam_file_path,
                                   args.counts_file, args.extension, args.sense, args.skip_plot,
                                   not args.quiet, args.number_of_threads, args.black_list)
    else:
        if args.counts_file:
            raise Exception("Appending to a prior regions counts.h5 file is not supported at this time -- "
                            "please email the developer if you need this feature")
        # non-exhaustive list of items that would need to be handled to get this working:
        ## review matrix output, specifically the assumption that each file has the exact same regions in the same order
        liquidator = RegionLiquidator(args.regions_file, args.output_directory, args.bam_file_path, 
                                      args.region_format, args.counts_file, args.extension, args.sense,
                                      not args.quiet, args.number_of_threads)

    if args.flatten:
        liquidator.flatten()

    if args.match_bamToGFF:
        if args.regions_file is None:
            logging.warning("Ignoring match_bamToGFF argument (this is only supported if a regions file is provided)")
        else:
            logging.info("Writing bamToGff style matrix.txt file")
            start = time()
            write_bamToGff_matrix(os.path.join(args.output_directory, "matrix.txt"), liquidator.counts_file_path)  
            duration = time() - start
            logging.info("Writing matrix.txt took %f seconds" % duration)
            liquidator.log_time('matrix', duration)

    if args.xml_timings:
        liquidator.write_timings_to_junit_xml()

if __name__ == "__main__":
    main()

'''
   The MIT License (MIT) 

   Copyright (c) 2013 John DiMatteo (jdimatteo@gmail.com)

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
