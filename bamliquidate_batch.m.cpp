#include "bamliquidator.h"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <sstream>
#include <fstream>

#include <boost/atomic.hpp>
#include <boost/thread/thread.hpp>
#include <boost/lockfree/queue.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

struct ChromosomeCounts
{
  const std::string chromosome; // e.g. chr1
  const std::string bam_file_name; // e.g. L228_121_DHL6_RNA_POL2.hg18.bwt.sorted.bam 
  const std::string cell_type; // e.g. dhl6
  const unsigned int bin_size = 100000; // 100K base pairs per bin
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

void batch()
{
  boost::lockfree::queue<int> computed_chromosome_count_indices(128);

  std::string bam_file_path("../copied_from_tod/");
  std::string bam_file_name("04032013_D1L57ACXX_4.TTAGGC.hg18.bwt.sorted.bam");
  std::string bam_file = bam_file_path+bam_file_name; 
  const std::vector<std::string> chromosomes {"chr1", "chr2", "chr3", "chr4", "chr5", 
    "chr6", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16",
    "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX"};

  const int binSize = 100000; // 100K

  const ChromosomeLengths lengths("../copied_from_tod/ucsc_chromSize.txt");

  for (auto chr : chromosomes)
  {
    int basePairsInChromosome = lengths(bam_file, chr);
    int binsInChromsome = std::ceil(basePairsInChromosome / (double) binSize);
    int maxBasePairForLiquidation = binsInChromsome * binSize;
    // pickup here: call liquidate with proper args
    // first just try a single thread counting from file,
    // and a single thread recording counts in hdf5
    // I'm not sure if samtools is thread safe -- try this after I confirm it working
    // I need to write a script to compare hdf5 and mysql results or something
    // maybe I should just write both to txt and diff them
    std::cout << "liquidate(" << bam_file << ", " << chr << ", 0, " << maxBasePairForLiquidation << ", '+', 1, 0)" << std::endl;
    std::vector<double> counts = liquidate(bam_file, chr, 0, maxBasePairForLiquidation, '+', 1, 0);
  }
}

int main()
{
  batch();
  return 0;
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
