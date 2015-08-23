#include "bamliquidator.h"
#include "liquidator_util.h"
#include "bamliquidator_regions.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

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

using namespace liquidator;

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
