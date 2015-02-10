#include <iostream>
#include <samtools/bam.h>

char complement(char c)
{
    switch(c)
    {
      case 'A': return 'T';
      case 'C': return 'G';
      case 'G': return 'C';
      case 'T': return 'A';
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

void liquidate(const std::string& input_bam_file, const std::string& target, const std::string& output_bam_file)
{
  const std::string target_complement = complement(target);

  bam1_t* read = bam_init1();
  bamFile input = bam_open(input_bam_file.c_str(), "r");
  bamFile output = bam_open(output_bam_file.c_str(), "w");
  bam_header_t* header = bam_header_read(input);

  if (read == 0 || input == 0 || header == 0)
  {
    throw std::runtime_error("failed to open " + input_bam_file);
  }
  if (output == 0) 
  {
    throw std::runtime_error("failed to open " + output_bam_file);
  }

  bam_header_write(output, header);

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

    if (sequence.find(target) != std::string::npos || sequence.find(target_complement) != std::string::npos)
    {
      bam_write1(output, read);
    }
  }

  bam_header_destroy(header);
  bam_close(input);
  bam_close(output);
  bam_destroy1(read);

  bam_index_build(output_bam_file.c_str());
}

int main(int argc, char** argv)
{
  // todo: support multiple strings in one go?
  if (argc != 4)
  {
    std::cerr << "Usage: " << argv[0] << " [INPUT_BAM] [STRING] [OUTPUT_BAM] " << std::endl;
    std::cerr << "e.g. " << argv[0] << " input.bam TGGGAA output.bam" << std::endl; 
    return 1;
  }

  const std::string input_bam_file = argv[1]; 
  const std::string target = argv[2];
  const std::string output_bam_file = argv[3]; 

  try 
  {
    liquidate(input_bam_file, target, output_bam_file);
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
