#include "bamliquidator_util.h"
#include "bamliquidator_regions.h"

#include <samtools/bam.h>

#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <utility>
#include <tuple>

typedef std::vector<std::pair<std::string, size_t>> MotifCounts;

static int bam_fetch_func(const bam1_t* read, void* void_triple)
{
  auto triple = *static_cast<std::tuple<const std::vector<std::string>&, size_t&, MotifCounts&>*>(void_triple);
  const auto& reverse_complements = std::get<0>(triple);
  size_t& read_count = std::get<1>(triple);
  MotifCounts& motif_counts = std::get<2>(triple);

  const bam1_core_t *c = &read->core;
  uint8_t *s = bam1_seq(read);

  static /*thread_local*/ std::string sequence; // todo: why doesn't thread local compile?
  if (sequence.size() != c->l_qseq)
  { 
    // assuming that all reads are uniform length, this will only happen once
    sequence = std::string(c->l_qseq, ' '); 
  }

  bool readHasN = false;
  for (int i = 0; i < c->l_qseq; ++i)
  {
    sequence[i] = bam_nt16_rev_table[bam1_seqi(s, i)];
    if (sequence[i] == 'N') 
    {
      readHasN = true; 
    }
  }
  if (readHasN)
  {
    return 0; // reads with Ns are considered poor and should not be checked
  }
  ++read_count;

  for (size_t i=0; i < motif_counts.size(); ++i)
  {
    motif_counts[i].second += count(sequence, motif_counts[i].first) + count(sequence, reverse_complements[i]);
  }

  return 0;
}

// returns total number of reads in the regions
size_t liquidate(const std::string& input_bam_file,
                 const std::vector<Region>& regions,
                 MotifCounts& motif_counts)
{
  std::vector<std::string> reverse_complements;
  for ( const auto& p : motif_counts )
  {
    std::string reverse_complement = complement(p.first);
    std::reverse(reverse_complement.begin(), reverse_complement.end());
    reverse_complements.push_back(reverse_complement);
  }

  bam1_t* read = bam_init1();
  bamFile input = bam_open(input_bam_file.c_str(), "r");
  bam_header_t* header = bam_header_read(input);
  bam_index_t* index = bam_index_load(input_bam_file.c_str());

  if (read == 0 || input == 0 || header == 0 || index == 0)
  {
    throw std::runtime_error("failed to open " + input_bam_file);
  }

  std::string sequence;
  size_t read_count = 0;

  std::tuple<const std::vector<std::string>&, size_t&, MotifCounts&> data(reverse_complements, read_count, motif_counts);

  for (const auto& r: regions)
  {
    // todo: don't I just need to parse_region once per chromosome to get the tid? perhaps there is a faster way to do this without parsing a whole region string?
    // todo: consider adding a util function to do this and remove duplicate code in bamliquidator.cpp
    std::stringstream coord;
    coord << r.chromosome << ':' << r.start << '-' << r.stop;

    int ref,beg,end;
    const int region_parse_rc = bam_parse_region(header,coord.str().c_str(), &ref, &beg, &end);
    if (region_parse_rc != 0)
    {
      std::stringstream error_msg;
      error_msg << "bam_parse_region failed with return code " << region_parse_rc;
      throw std::runtime_error(error_msg.str());
    }
    if(ref<0)
    {
      // this bam doesn't have this chromosome
      continue;
    }

    const int fetch_rc = bam_fetch(input, index, ref, beg, end, &data, bam_fetch_func);
    if (fetch_rc != 0)
    {
      std::stringstream error_msg;
      error_msg << "bam_fetch failed with return code " << fetch_rc;
      throw std::runtime_error(error_msg.str());
    }
  }

  bam_index_destroy(index);
  bam_header_destroy(header);
  bam_close(input);
  bam_destroy1(read);

  return read_count;
}

int main(int argc, char** argv)
{
  if (argc < 5)
  {
    std::cerr << "Usage: " << argv[0] << " [REGIONS] [BACKGROUND_BAM] [TARGET_BAM] [MOTIF_1] [MOTIF_2] ... [MOTIF_N]" << std::endl;
    std::cerr << "e.g. " << argv[0] << " regions.bed background.bam input.bam TGGGAA AGGG" << std::endl;
    return 1;
  }

  const std::vector<Region> regions(parse_regions(argv[1], "bed", 0)); 
  const std::string background_bam_file = argv[2]; 
  const std::string target_bam_file = argv[3]; 

  MotifCounts background_motif_counts;
  for (int i=4; i < argc; ++i)
  {
    background_motif_counts.push_back(std::make_pair<std::string, size_t>(argv[i], 0));
  }
  MotifCounts target_motif_counts(background_motif_counts);

  try 
  {
    size_t background_count = liquidate(background_bam_file, regions, background_motif_counts);
    size_t target_count     = liquidate(target_bam_file,     regions, target_motif_counts);

    std::cout << "motif\tbackground (normalized)\ttarget (normalized)\n";
    for (size_t i=0; i < background_motif_counts.size(); ++i)
    {
      const std::string& motif = background_motif_counts[i].first;
      if (motif != target_motif_counts[i].first)
      {
        throw std::runtime_error("internal logic error");
      }

      size_t background_matches = background_motif_counts[i].second; 
      size_t target_matches     = target_motif_counts[i].second;
      double background_matches_normalized = background_matches / ( double(background_count) / 1E6);
      double target_matches_normalized     = target_matches     / ( double(target_count)     / 1E6);
      std::cout << motif << "\t" << background_matches << " (" << background_matches_normalized << ")\t" 
                                 << target_matches     << " (" << target_matches_normalized     << ")\n";
    }
    std::cout << std::endl 
              << "background reads: " << background_count << std::endl
              << "target reads: " << target_count << std::endl;
  }
  catch(const std::exception& e)
  {
    std::cerr << e.what() << std::endl;
    return 2;
  }

  return 0;
}

/* The MIT License (MIT) 

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
 */
