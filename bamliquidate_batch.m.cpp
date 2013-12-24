#include "bamliquidator.h"

#include <boost/thread/thread.hpp>
#include <boost/lockfree/queue.hpp>
#include <iostream>
#include <sstream>
#include <cmath>

#include <boost/atomic.hpp>

struct ChromosomeCounts
{
  const std::string chromosome; // e.g. chr1
  const std::string bam_file_name; // e.g. L228_121_DHL6_RNA_POL2.hg18.bwt.sorted.bam 
  const std::string cell_type; // e.g. dhl6
  const unsigned int bin_size = 100000; // 100K base pairs per bin
  const std::vector<double> bin_counts;
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

  for (auto chr : chromosomes)
  {
    std::stringstream coord;
    int basePairsInChromosome = 247249719;
    int binsInChromsome = std::ceil(basePairsInChromosome / (double) binSize);
    int maxBasePairForLiquidation = binsInChromsome * binSize;
    coord << chr << ":0-" << maxBasePairForLiquidation;
    // pickup here: call liquidate with proper args
    // coordinate is of the form "chr1:0-10000"
    // first just try a single thread counting from file,
    // and a single thread recording counts in hdf5
    // I'm not sure if samtools is thread safe -- try this after I confirm it working
    // I need to write a script to compare hdf5 and mysql results or something
    // maybe I should just write both to txt and diff them
    std::cout << "counting " << coord.str() << std::endl;
    std::vector<double> counts = liquidate(bam_file, coord.str(), 0, maxBasePairForLiquidation, '+', 1, 0);
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
