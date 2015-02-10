#include <iostream>
#include <samtools/bam.h>

void test(const std::string& bam_file, const std::string& target)
{
  bam1_t* read = bam_init1(); // what is this for?
  bamFile file_handle = bam_open(bam_file.c_str(), "r");
  bam_header_t* header = bam_header_read(file_handle);

  if (file_handle == 0 || read == 0 || header == 0) throw std::runtime_error("failed to open " + bam_file);

  size_t count = 0;
  std::cout << "Hello?" << std::endl;
  while ( bam_read1(file_handle, read) >= 0 )
  {
    ++count;
  }
  std::cout << "count: " << count << std::endl;

  bam_header_destroy(header);
  bam_close(file_handle);
  bam_destroy1(read);
}

int main(int argc, char** argv)
{
  // todo: support multiple strings in one go?
  if (argc != 3)
  {
    std::cout << "Usage: " << argv[0] << " [BAM_FILE] [STRING]" << std::endl;
    std::cout << "e.g. " << argv[0] << " 20110819_580_hg19.sorted.bam TGGGAA" << std::endl; 
    return 1;
  }

  const std::string bam_file = argv[1]; 
  const std::string target = argv[2];

  test(bam_file, target);

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

