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

struct Region
{
  char chromosome[16];
  char name[64];
  uint64_t start;
  uint64_t stop;
  char strand;
  uint64_t count;
};

std::ostream& operator<<(std::ostream& os, const Region& r)
{
    os << r.chromosome << ' ' << r.name << ' ' << r.start << " -> " << r.stop << ' ' << r.strand;
    return os;
}

std::vector<Region> parse_regions(const std::string& region_file_path)
{
  // todo: support both gff and bed, for now just gff
  std::vector<Region> regions;

  const int chromosome_column = 0;
  const int name_column = 1;
  const int start_column = 3;
  const int stop_column = 4;
  const int strand_column = 6;
  const int min_columns = 7;

  std::ifstream region_file(region_file_path);
  if (!region_file.is_open())
  {
    throw std::runtime_error("failed to open region_file " + region_file_path);
  }

  for(std::string line; std::getline(region_file, line); )
  {
    std::vector<std::string> columns;
    boost::split(columns, line, boost::is_any_of("\t"));
    if (columns.size() < min_columns)
    {
      throw std::runtime_error("error parsing " + region_file_path);
    }
    Region region;
    strncpy(region.chromosome, columns[chromosome_column].c_str(), sizeof(Region::chromosome));
    strncpy(region.name,       columns[name_column].c_str(),       sizeof(Region::name));
    region.start = boost::lexical_cast<uint64_t>(columns[start_column]);
    region.stop  = boost::lexical_cast<uint64_t>(columns[stop_column]);

    if (columns[strand_column].size() != 1)
    {
      throw std::runtime_error("error parsing strand: '" + columns[strand_column] + "'");
    }
    region.strand = columns[strand_column][0];

    regions.push_back(region);
  }

  return regions;
}

void liquidate_regions(hid_t& file, std::vector<Region>& regions, const std::string& bam_file_path)
{
  // todo: open bam_file here, instead of repeatedly inside loop
 
  for (size_t i=0; i < regions.size(); ++i)
  {
    try
    {
      std::cout << "about to liquidate " << i << ": " << regions[i] << std::endl;
      std::vector<double> counts = liquidate(bam_file_path,
                                             regions[i].chromosome,
                                             regions[i].start, 
                                             regions[i].stop, 
                                             regions[i].strand, 
                                             1, 0);
      if (counts.size() != 0)
      {
        std::cout << "about to throw..." << std::endl;
        throw std::runtime_error("liquidate failed to provide exactly one count (count is ");// +
        //  boost::lexical_cast<std::string>(counts.size()) + ")");
      }
      std::cout << "done liquidated " << i << std::endl;
    } catch(const std::exception& e)
    {
      std::cerr << "Skipping region " << i+1 << " (" << regions[i] << ") due to error: "
                << e.what() << std::endl;
    }
  }
}

int main(int argc, char* argv[])
{
  try
  {
    if (argc != 4)
    {
      std::cerr << "usage: " << argv[0] << " gff_or_bed_file bam_file hdf5_file\n"
        << "\ne.g. " << argv[0] << " /grail/annotations/HG19_SUM159_BRD4_-0_+0.gff "
        << "\n      /ifs/labs/bradner/bam/hg18/mm1s/04032013_D1L57ACXX_4.TTAGGC.hg18.bwt.sorted.bam\n"
        << "\nnote that this application is intended to be run from bamliquidator_batch.py -- see"
        << "\nhttps://github.com/BradnerLab/pipeline/wiki for more information"
        << std::endl;
      return 1;
    }

    const std::string region_file_path = argv[1];
    const std::string bam_file_path = argv[2];
    const std::string hdf5_file_path = argv[3];

    hid_t h5file = 0; //todo:H5Fopen(hdf5_file_path.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    if (h5file < 0)
    {
      std::cerr << "Failed to open H5 file " << hdf5_file_path << std::endl;
      return 3;
    }

    std::vector<Region> regions = parse_regions(region_file_path);

    liquidate_regions(h5file, regions, bam_file_path);
   
    //todo: H5Fclose(h5file);

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
