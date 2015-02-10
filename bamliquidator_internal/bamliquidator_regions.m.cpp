#include "bamliquidator.h"
#include "bamliquidator_util.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include <hdf5.h>
#include <hdf5_hl.h>

#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>
#include <tbb/task_scheduler_init.h>

//#define time_region_parsing
#ifdef time_region_parsing 
#include <boost/timer/timer.hpp>
#endif

const size_t region_name_length = 64;

// this Region must match exactly the structure in HDF5
// -- see bamliquidator_batch.py function create_regions_table
struct Region
{
  uint32_t bam_file_key;
  char chromosome[64];
  char region_name[region_name_length];
  uint64_t start;
  uint64_t stop;
  char strand;
  uint64_t count;
  double normalized_count;

  bool is_valid(const std::map<std::string, size_t>& chromosome_to_length)
  {
    const auto it = chromosome_to_length.find(chromosome);
    if ( it != chromosome_to_length.end() )
    {
      if ( stop <= it->second )
      {
        return true;
      }
    }

    return false;
  }
};

std::ostream& operator<<(std::ostream& os, const Region& r)
{
  os << "bam file key " << r.bam_file_key << ' ' << r.chromosome << ' '
     << r.region_name << ' ' << r.start << " -> " << r.stop << ' '
     << r.strand << ' ' << r.normalized_count;
  return os;
}

// default_strand: optional argument, default _ indicates to use 
//                 gff strand column or . (both) for .bed region file
std::vector<Region> parse_regions(const std::string& region_file_path,
                                  const std::string& region_format,
                                  const unsigned int bam_file_key,
                                  const std::map<std::string, size_t>& chromosome_to_length, 
                                  const char default_strand = '_') 
{
  int chromosome_column = 0;
  int name_column = 0;
  int start_column = 0;
  int stop_column = 0;
  int strand_column = 0;
  unsigned int min_columns = 0;

  if (region_format == "gff")
  {
    chromosome_column = 0;
    name_column = 1;
    start_column = 3;
    stop_column = 4;
    strand_column = 6;
    min_columns = 7;
  }
  else if (region_format == "bed")
  {
    chromosome_column = 0;
    name_column = 3;
    start_column = 1;
    stop_column = 2;
    strand_column = 5;
    min_columns = 3;
  }
  else
  {
    throw std::runtime_error("unsupported region file format (" + region_format + " for " + region_file_path
                             + "), please supply a .gff or .bed file");
  }

  std::ifstream region_file(region_file_path.c_str());
  if (!region_file.is_open())
  {
    throw std::runtime_error("failed to open region_file " + region_file_path);
  }

  std::vector<Region> regions;
  int line_number = 1;
  for(std::string line; std::getline(region_file, line); ++line_number)
  {
    std::vector<std::string> columns;
    boost::split(columns, line, boost::is_any_of("\t"));
    if (columns.size() < min_columns)
    {
      std::stringstream ss;
      ss << "Not enough columns parsing line " << line_number << " '" << line << "' of " << region_file_path;
      throw std::runtime_error(ss.str());
    }
    Region region;
    region.bam_file_key = bam_file_key; 
    copy(region.chromosome,  columns[chromosome_column], sizeof(Region::chromosome));
    if (columns.size() > name_column)
    {
      copy(region.region_name, columns[name_column], sizeof(Region::region_name));
      if (columns[name_column].size() >= sizeof(Region::region_name))
      {
        Logger::warn() << "Truncated region on line " << line_number << " from '" << columns[name_column] << "' to '" << region.region_name << "'";
      }
    }
    else
    {
      copy(region.region_name, "", sizeof(Region::region_name));
    }
    region.start = boost::lexical_cast<uint64_t>(columns[start_column]);
    region.stop  = boost::lexical_cast<uint64_t>(columns[stop_column]);
    if (region.start > region.stop)
    {
      std::swap(region.start, region.stop);
    }

    if (columns.size() > strand_column)
    {
      if (columns[strand_column].size() != 1)
      {
        std::stringstream ss;
        ss << "error parsing strand: '" << columns[strand_column] << "' on line " << line_number;
        throw std::runtime_error(ss.str());
      }
      region.strand = default_strand == '_'
                    ? columns[strand_column][0]
                    : default_strand;
    }
    else
    {
      region.strand = default_strand == '_'
                    ? '.'
                    : default_strand; 
    }
    region.count = 0;
    region.normalized_count = 0.0;

    if (region.is_valid(chromosome_to_length))
    {
      regions.push_back(region);
    }
    else
    {
      Logger::warn() << "Excluding invalid region on line " << line_number << ": " << region;
    }
  }

  return regions;
}

void write(hid_t& file, std::vector<Region>& regions)
{
  const size_t record_size = sizeof(Region);

  size_t record_offset[] = { HOFFSET(Region, bam_file_key),
                             HOFFSET(Region, chromosome),
                             HOFFSET(Region, region_name),
                             HOFFSET(Region, start),
                             HOFFSET(Region, stop),
                             HOFFSET(Region, strand),
                             HOFFSET(Region, count),
                             HOFFSET(Region, normalized_count) };

  size_t field_sizes[] = { sizeof(Region::bam_file_key),
                           sizeof(Region::chromosome),
                           sizeof(Region::region_name),
                           sizeof(Region::start),
                           sizeof(Region::stop),
                           sizeof(Region::strand),
                           sizeof(Region::count),
                           sizeof(Region::normalized_count) };

  herr_t status = H5TBappend_records(file, "region_counts", regions.size(), record_size, record_offset,
                                     field_sizes, regions.data());
  if (status != 0)
  {
    std::stringstream ss;
    ss << "Error appending record, status = " << status;
    throw std::runtime_error(ss.str());
  }
}

class Liquidator 
{
public:
  Liquidator(const std::string& bam_file_path):
    bam_file_path(bam_file_path),
    fp(nullptr),
    bamidx(nullptr)
  {
    init();
  }

  Liquidator(const Liquidator& other):
    bam_file_path(other.bam_file_path),
    fp(nullptr),
    bamidx(nullptr)
  {
    init();
  }

  Liquidator& operator=(const Liquidator& other) = delete;

  ~Liquidator()
  {
    bam_index_destroy(bamidx);
    samclose(fp);
  }

  double liquidate(const std::string& chromosome, int start, int stop, char strand, unsigned int extension)
  {
    std::vector<double> counts = ::liquidate(fp, bamidx, chromosome, start, stop, strand, 1, extension);
    if (counts.size() != 1)
    {
      throw std::runtime_error("liquidate failed to provide exactly one count (count is " +
        boost::lexical_cast<std::string>(counts.size()) + ")");
    }
    return counts[0];
  }

private:
  std::string bam_file_path;
  samfile_t* fp;
  bam_index_t* bamidx;

  void init()
  {
    fp = samopen(bam_file_path.c_str(),"rb",0);
    if(fp == NULL)
    {
      throw std::runtime_error("samopen() error with " + bam_file_path);
    }

    bamidx = bam_index_load(bam_file_path.c_str());
    if (bamidx == NULL)
    {
      throw std::runtime_error("bam_index_load() error with " + bam_file_path);
    }
  }
};

typedef tbb::enumerable_thread_specific<Liquidator,
                                        tbb::cache_aligned_allocator<Liquidator>,
                                        tbb::ets_key_per_instance>
        Liquidators;

void liquidate_regions(std::vector<Region>& regions, const std::string& bam_file_path,
                       size_t region_begin, size_t region_end, unsigned int extension,
                       Liquidators& liquidators)
{
  Liquidator& liquidator = liquidators.local();

  for (size_t i=region_begin; i < region_end; ++i)
  {
    try
    {
      regions[i].count = liquidator.liquidate(regions[i].chromosome,
                                              regions[i].start, 
                                              regions[i].stop, 
                                              regions[i].strand,
                                              extension);
    } catch(const std::exception& e)
    {
      Logger::error() << "Aborting because failed to parse region " << i+1 << " (" << regions[i] << ") due to error: "
                      << e.what();
      throw;
    }
  }
}

void liquidate_and_write(hid_t& file, std::vector<Region>& regions,
                         unsigned int extension, const std::string& bam_file_path)
{
  Liquidators liquidators((Liquidator(bam_file_path))); 

  tbb::parallel_for(
    tbb::blocked_range<int>(0, regions.size(), 1),
    [&](const tbb::blocked_range<int>& range)
    {
      liquidate_regions(regions, bam_file_path, range.begin(), range.end(), extension, liquidators);
    },
    tbb::auto_partitioner());

  write(file, regions);
}

int main(int argc, char* argv[])
{
  try
  {
    if (argc < 13 || argc % 2 != 1)
    {
      std::cerr << "usage: " << argv[0] << " number_of_threads region_file gff_or_bed_format extension bam_file bam_file_key hdf5_file "
                << "log_file write_warnings_to_stderr strand chr1 length1 ...\n"
        << "\ne.g. " << argv[0] << " /grail/annotations/HG19_SUM159_BRD4_-0_+0.gff gff"
        << "\n      /ifs/labs/bradner/bam/hg18/mm1s/04032013_D1L57ACXX_4.TTAGGC.hg18.bwt.sorted.bam 137 counts.hdf5 "
        << "\n      output/log.txt 1 _ chr1 247249719 chr2 242951149 chr3 199501827\n"
        << "\nstrand value of _ means use strand that is specified in region file (and use . if strand not specified in region file)."
        << "\nnumber of threads <= 0 means use a number of threads equal to the number of logical cpus."
        << "\nnote that this application is intended to be run from bamliquidator_batch.py -- see"
        << "\nhttps://github.com/BradnerLab/pipeline/wiki for more information"
        << std::endl;
      return 1;
    }

    const int number_of_threads = boost::lexical_cast<int>(argv[1]);
    const std::string region_file_path = argv[2];
    const std::string region_format = argv[3];
    const unsigned int extension = boost::lexical_cast<unsigned int>(argv[4]);
    const std::string bam_file_path = argv[5];
    const unsigned int bam_file_key = boost::lexical_cast<unsigned int>(argv[6]);
    const std::string hdf5_file_path = argv[7];
    const std::string log_file_path = argv[8];
    const bool write_warnings_to_stderr = boost::lexical_cast<bool>(argv[9]);
    const char strand = boost::lexical_cast<char>(argv[10]);
    const std::vector<std::pair<std::string, size_t>> chromosome_lengths = extract_chromosome_lengths(argc, argv, 11);

    tbb::task_scheduler_init init( number_of_threads <= 0 
                                 ? tbb::task_scheduler_init::automatic
                                 : number_of_threads); 

    Logger::configure(log_file_path, write_warnings_to_stderr);

    hid_t h5file = H5Fopen(hdf5_file_path.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    if (h5file < 0)
    {
      Logger::error() << "Failed to open H5 file " << hdf5_file_path;
      return 3;
    }

    #ifdef time_region_parsing 
    boost::timer::cpu_timer timer; 
    #endif

    std::map<std::string, size_t> chromosome_to_length;
    for (auto& chr_length : chromosome_lengths)
    {
      chromosome_to_length[chr_length.first] = chr_length.second;
    }

    std::vector<Region> regions = parse_regions(region_file_path,
                                                region_format,
                                                bam_file_key,
                                                chromosome_to_length,
                                                strand);
    #ifdef time_region_parsing 
    timer.stop();
    std::cout << "parsing regions took" << timer.format() << std::endl;
    #endif
    if (regions.size() == 0)
    {
      Logger::warn() << "No valid regions detected in " << region_file_path;
      return 0;
    }

    liquidate_and_write(h5file, regions, extension, bam_file_path);
   
    H5Fclose(h5file);

    return 0;
  }
  catch(const std::exception& e)
  {
    Logger::error() << "Unhandled exception: " << e.what();

    return 4; 
  }
}

/* The MIT License (MIT) 

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
 */
