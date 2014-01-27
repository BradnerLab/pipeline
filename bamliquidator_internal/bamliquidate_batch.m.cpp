#include "bamliquidator.h"

#include <cmath>
#include <fstream>
#include <future>
#include <iostream>
#include <map>
#include <stdexcept>
#include <sstream>

#include <boost/atomic.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include <hdf5.h>
#include <hdf5_hl.h>

/* todo:
  
   - don't store strings in table records, e.g. have seperate tables for chromosomes, file names,
     and cell types.  in the count table records, just store integers representing the strings.
     this should speed up the writing and cut the file sizes
 */

namespace
{
  const bool logging_enabled = false;
}

#define LOG_INFO(streamable) \
if(logging_enabled) \
{ \
  std::cout << streamable << std::endl; \
}

struct ChromosomeCounts
{
  const std::string chromosome; // e.g. chr1
  const std::string bam_file_name; // e.g. L228_121_DHL6_RNA_POL2.hg18.bwt.sorted.bam 
  const std::string cell_type; // e.g. dhl6
  const unsigned int bin_size; // e.g. 100000 base pairs per bin
  const std::vector<double> bin_counts;
};

class ChromosomeLengths
{
public:
  ChromosomeLengths(const std::string &chrom_size_file)
  {
    const size_t chromosome_column = 0; 
    const size_t type_column = 2;
    const size_t length_column = 4;

    std::ifstream size_file(chrom_size_file);
    if (!size_file.is_open())
    {
      throw std::runtime_error("failed to open chrom_size_file " + chrom_size_file);
    }

    for(std::string line; std::getline(size_file, line); )
    {
      std::vector<std::string> columns;
      boost::split(columns, line, boost::is_any_of("\t"));
      if (columns.size() < length_column + 1)
      {
        throw std::runtime_error("error parsing chrom_size_file");
      }
      m_type_to_chromosome_to_count[columns[type_column]][columns[chromosome_column]] = 
        boost::lexical_cast<size_t>(columns[length_column]);
    }
  }

  size_t operator()(const std::string &bam_file, const std::string& chromosome) const
  {
    for (const auto type : m_type_to_chromosome_to_count)
    {
      if (bam_file.find(type.first) != std::string::npos)
      {
        const auto count = type.second.find(chromosome);
        if (count == type.second.end())
        {
          throw std::runtime_error("failed to find chromosome " + chromosome
                                   + " for " + type.first);
        }
        return count->second;
      }
    }
    throw std::runtime_error("failed to determine type of " + bam_file);
  }

private:
  typedef std::map<std::string, std::map<std::string, size_t>> TypeToChromosomeToCount;
  TypeToChromosomeToCount m_type_to_chromosome_to_count;
};

struct CountH5Record
{
  uint32_t bin_number;
  char cell_type[16];
  char chromosome[16];
  uint64_t count;
  char file_name[64];
};

ChromosomeCounts count(const std::string& chr,
                       const std::string& cell_type,
                       const unsigned int bin_size,
                       const ChromosomeLengths& lengths,
                       const std::string& bam_file)
{
  const size_t last_slash_position = bam_file.find_last_of("/");
  const std::string bam_file_name = last_slash_position == std::string::npos 
                                  ? bam_file
                                  : bam_file.substr(last_slash_position + 1);

  int base_pairs = lengths(bam_file, chr);
  int bins = std::ceil(base_pairs / (double) bin_size);
  int max_base_pair = bins * bin_size;
  LOG_INFO(" - counting chromosome " << chr);
  return ChromosomeCounts {chr, bam_file_name, cell_type, bin_size, 
    liquidate(bam_file, chr, 0, max_base_pair, '.', bins, 0)};
}

void batch(hid_t& file,
           const std::string& cell_type,
           const unsigned int bin_size,
           const ChromosomeLengths& lengths,
           const std::string& bam_file)
{
  std::deque<std::future<ChromosomeCounts>> future_counts;

  const std::vector<std::string> chromosomes {"chr1", "chr2", "chr3", "chr4", "chr5", 
    "chr6", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16",
    "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX"};

  for (auto& chr : chromosomes)
  {
    future_counts.push_back(std::async(count,
                                       std::ref(chr),
                                       std::ref(cell_type),
                                       bin_size,
                                       std::ref(lengths),
                                       std::ref(bam_file)));
  }

  const size_t record_size = sizeof(CountH5Record);

  size_t record_offset[] = { HOFFSET(CountH5Record, bin_number),
                             HOFFSET(CountH5Record, cell_type),
                             HOFFSET(CountH5Record, chromosome),
                             HOFFSET(CountH5Record, count),
                             HOFFSET(CountH5Record, file_name) };

  size_t field_sizes[] = { sizeof(CountH5Record::bin_number),
                           sizeof(CountH5Record::cell_type),
                           sizeof(CountH5Record::chromosome),
                           sizeof(CountH5Record::count),
                           sizeof(CountH5Record::file_name) };

  for (auto& future_count : future_counts)
  {
    ChromosomeCounts counts = future_count.get();

    CountH5Record record;
    record.bin_number = 0;
    strncpy(record.cell_type,  counts.cell_type.c_str(),     sizeof(CountH5Record::cell_type));
    strncpy(record.chromosome, counts.chromosome.c_str(),    sizeof(CountH5Record::chromosome));
    strncpy(record.file_name,  counts.bam_file_name.c_str(), sizeof(CountH5Record::file_name));

    std::vector<CountH5Record> records(counts.bin_counts.size(), record);
    for (size_t bin=0; bin <= counts.bin_counts.size(); ++bin)
    {
      records[bin].bin_number = bin;
      records[bin].count = counts.bin_counts[bin];
    }
   
    LOG_INFO(" - writing " << counts.bin_counts.size() << " from " << counts.chromosome);

    herr_t status = H5TBappend_records(file, "counts", records.size(), record_size, record_offset,
                                       field_sizes, records.data());
    if (status != 0)
    {
      std::cerr << "Error appending record, status = " << status << std::endl;
    }
  }
}

int main(int argc, char* argv[])
{
  try
  {
    if (argc != 6)
    {
      std::cerr << "usage: " << argv[0] << " cell_type bin_size ucsc_chrom_size_path bam_file_path hdf5_file\n"
        << "\ne.g. " << argv[0] << " mm1s 100000 /grail/annotations/ucsc_chromSize.txt"
        << "\n      /ifs/labs/bradner/bam/hg18/mm1s/04032013_D1L57ACXX_4.TTAGGC.hg18.bwt.sorted.bam\n"
        << "\nnote that this application is intended to be run from bamliquidator_batch.py -- see"
        << "\nhttps://github.com/BradnerLab/pipeline/wiki for more information"
        << std::endl;
      return 1;
    }

    const std::string cell_type = argv[1];
    const unsigned int bin_size = boost::lexical_cast<unsigned int>(argv[2]);
    const std::string ucsc_chrom_size_path = argv[3];
    const std::string bam_file_path = argv[4];
    const std::string hdf5_file_path = argv[5];

    if (bin_size == 0)
    {
      std::cerr << "bin size cannot be zero" << std::endl;
      return 2;
    }

    const ChromosomeLengths lengths(ucsc_chrom_size_path);

    hid_t h5file = H5Fopen(hdf5_file_path.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    if (h5file < 0)
    {
      std::cerr << "Failed to open H5 file " << hdf5_file_path << std::endl;
      return 3;
    }

    batch(h5file, cell_type, bin_size, lengths, bam_file_path);
   
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
 */
