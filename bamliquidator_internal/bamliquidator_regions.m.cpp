#include "bamliquidator.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <sstream>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include <hdf5.h>
#include <hdf5_hl.h>

#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>
#include <tbb/task_scheduler_init.h>

// this Region must match exactly the structure in HDF5
// -- see bamliquidator_batch.py function create_regions_table
struct Region
{
  char file_name[64];
  char chromosome[16];
  char region_name[64];
  uint64_t start;
  uint64_t stop;
  char strand;
  uint64_t count;
  double normalized_count;
};

std::ostream& operator<<(std::ostream& os, const Region& r)
{
  os << r.file_name << ' ' << r.chromosome << ' ' << r.region_name << ' ' 
     << r.start << " -> " << r.stop << ' ' << r.strand << ' ' << r.normalized_count;
  return os;
}

// default_strand: optional argument, default (space) indicates to use 
//                 gff strand column or . (both) for .bed region file
std::vector<Region> parse_regions(const std::string& region_file_path,
                                  const std::string& file_name,
                                  const char default_strand = ' ') 
{
  const std::string extension = extension_from_file_name(region_file_path);

  int chromosome_column = 0;
  int name_column = 0;
  int start_column = 0;
  int stop_column = 0;
  int strand_column = 0;
  int min_columns = 0;

  if (extension == "gff")
  {
    chromosome_column = 0;
    name_column = 1;
    start_column = 3;
    stop_column = 4;
    strand_column = 6;
    min_columns = 7;
  }
  else if (extension == "bed")
  {
    chromosome_column = 0;
    name_column = 3;
    start_column = 1;
    stop_column = 2;
    strand_column = -1;
    min_columns = 4;
  }
  else
  {
    throw std::runtime_error("unsupported region file format (" + region_file_path
                             + "), please supply a .gff or .bed file");
  }

  std::ifstream region_file(region_file_path.c_str());
  if (!region_file.is_open())
  {
    throw std::runtime_error("failed to open region_file " + region_file_path);
  }

  std::vector<Region> regions;
  for(std::string line; std::getline(region_file, line); )
  {
    std::vector<std::string> columns;
    boost::split(columns, line, boost::is_any_of("\t"));
    if (columns.size() < min_columns)
    {
      throw std::runtime_error("error parsing " + region_file_path);
    }
    Region region;
    strncpy(region.file_name,   file_name.c_str(),                  sizeof(Region::file_name));
    strncpy(region.chromosome,  columns[chromosome_column].c_str(), sizeof(Region::chromosome));
    strncpy(region.region_name, columns[name_column].c_str(),       sizeof(Region::region_name));
    region.start = boost::lexical_cast<uint64_t>(columns[start_column]);
    region.stop  = boost::lexical_cast<uint64_t>(columns[stop_column]);
    if (region.start > region.stop)
    {
      std::swap(region.start, region.stop);
    }

    if (strand_column == -1)
    {
      region.strand = default_strand == ' '
                    ? '.'
                    : default_strand; 
    }
    else
    {
      if (columns[strand_column].size() != 1)
      {
        throw std::runtime_error("error parsing strand: '" + columns[strand_column] + "'");
      }
      region.strand = default_strand == ' '
                    ? columns[strand_column][0]
                    : default_strand;
    }

    regions.push_back(region);
  }

  return regions;
}

void write(hid_t& file, std::vector<Region>& regions)
{
  const size_t record_size = sizeof(Region);

  size_t record_offset[] = { HOFFSET(Region, file_name),
                             HOFFSET(Region, chromosome),
                             HOFFSET(Region, region_name),
                             HOFFSET(Region, start),
                             HOFFSET(Region, stop),
                             HOFFSET(Region, strand),
                             HOFFSET(Region, count),
                             HOFFSET(Region, normalized_count) };

  size_t field_sizes[] = { sizeof(Region::file_name),
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
    std::cerr << "Error appending record, status = " << status << std::endl;
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
      std::cerr << "Skipping region " << i+1 << " (" << regions[i] << ") due to error: "
                << e.what() << std::endl;
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
  tbb::task_scheduler_init init; 

  try
  {
    if (argc != 5 && argc != 6)
    {
      std::cerr << "usage: " << argv[0] << " gff_or_bed_file extension bam_file hdf5_file [strand]\n"
        << "\ne.g. " << argv[0] << " /grail/annotations/HG19_SUM159_BRD4_-0_+0.gff "
        << "\n      /ifs/labs/bradner/bam/hg18/mm1s/04032013_D1L57ACXX_4.TTAGGC.hg18.bwt.sorted.bam\n"
        << "\nnote that this application is intended to be run from bamliquidator_batch.py -- see"
        << "\nhttps://github.com/BradnerLab/pipeline/wiki for more information"
        << std::endl;
      return 1;
    }

    const std::string region_file_path = argv[1];
    const unsigned int extension = boost::lexical_cast<unsigned int>(argv[2]);
    const std::string bam_file_path = argv[3];
    const std::string hdf5_file_path = argv[4];
    const char strand = argc == 6
                      ? boost::lexical_cast<char>(argv[5])
                      : ' ';

    hid_t h5file = H5Fopen(hdf5_file_path.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    if (h5file < 0)
    {
      std::cerr << "Failed to open H5 file " << hdf5_file_path << std::endl;
      return 3;
    }

    std::vector<Region> regions = parse_regions(region_file_path,
                                                file_name_from_path(bam_file_path),
                                                strand);

    liquidate_and_write(h5file, regions, extension, bam_file_path);
   
    H5Fclose(h5file);

    return 0;
  }
  catch(const std::exception& e)
  {
    std::cerr << "Unhandled exception: " << e.what() << std::endl;

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
