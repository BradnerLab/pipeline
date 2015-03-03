#include <samtools/bam.h>

#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <utility>

char complement(char c)
{
    switch(c)
    {
      case 'A': return 'T';
      case 'C': return 'G';
      case 'G': return 'C';
      case 'T': return 'A';
      case 'N': return 'N';
    }
    throw std::runtime_error("no known complement for " + std::string(1, c));
}

std::string complement(const std::string& sequence)
{
  std::string c = sequence;
  for (size_t i=0; i < c.size(); ++i)
  {
    c[i] = complement(sequence[i]);
  }
  return c;
}

bool contains(const std::string& haystack, const std::string& needle)
{
  // naive implementation just to get started
  const char match_any_char = 'N';

  for (size_t i=0; i < haystack.size(); ++i)
  {
    for (size_t j=0; j < needle.size() && (j + i) < haystack.size(); ++j)
    {
      if (needle[j] != match_any_char && haystack[j + i] != match_any_char)
      {
        if (needle[j] != haystack[j + i])
        {
          break;
        }
      }
      if (j == needle.size() - 1)
      {
        return true;
      }
    }
  }
  return false;
}

void liquidate(const std::string& input_bam_file, std::vector<std::pair<std::string, size_t>>& motif_counts)
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

  if (read == 0 || input == 0 || header == 0)
  {
    throw std::runtime_error("failed to open " + input_bam_file);
  }

  std::string sequence;

  while (bam_read1(input, read) >= 0)
  {
    const bam1_core_t *c = &read->core;
    uint8_t *s = bam1_seq(read);

    // [s, s+c->l_qseq) is the sequence, with two bases packed into each byte.
    // I bet we could directly search that instead of first copying into a string
    // but lets get something simple working first. An intermediate step could be
    // to search integers without using bam_nt16_rev_table (and I wouldn't have
    // to worry about the packing complexity).

    if (sequence.size() != c->l_qseq)
    { 
      // assuming that all reads are uniform length, this will only happen once
      sequence = std::string(c->l_qseq, ' '); 
    }

    for (int i = 0; i < c->l_qseq; ++i)
    {
      sequence[i] = bam_nt16_rev_table[bam1_seqi(s, i)];
    }


    for (size_t i=0; i < motif_counts.size(); ++i)
    {
      if (contains(sequence, motif_counts[i].first) || contains(sequence, reverse_complements[i]))
      {
        ++(motif_counts[i].second);
      }
    }
  }

  bam_header_destroy(header);
  bam_close(input);
  bam_destroy1(read);
}

int main(int argc, char** argv)
{
  if (argc < 4)
  {
    std::cerr << "Usage: " << argv[0] << " [BACKGROUND_BAM] [TARGET_BAM] [MOTIF_1] [MOTIF_2] ... [MOTIF_N]" << std::endl;
    std::cerr << "e.g. " << argv[0] << " background.bam input.bam TGGGAA AGGG" << std::endl; 
    return 1;
  }

  const std::string background_bam_file = argv[1]; 
  const std::string target_bam_file = argv[2]; 

  std::vector<std::pair<std::string, size_t>> background_motif_counts;
  for (int i=3; i < argc; ++i)
  {
    background_motif_counts.push_back(std::make_pair<std::string, size_t>(argv[i], 0));
  }
  std::vector<std::pair<std::string, size_t>> target_motif_counts(background_motif_counts);

  try 
  {
    liquidate(background_bam_file, background_motif_counts);
    liquidate(target_bam_file, target_motif_counts);
    std::cout << "motif\tbackground\ttarget\n";
    for (size_t i=0; i < background_motif_counts.size(); ++i)
    {
      const std::string& motif = background_motif_counts[i].first;
      if (motif != target_motif_counts[i].first)
      {
        throw std::runtime_error("internal logic error");
      }

      std::cout << motif << "\t" << background_motif_counts[i].second << "\t" << target_motif_counts[i].second << std::endl;
    }
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
