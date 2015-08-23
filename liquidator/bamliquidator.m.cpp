#include <stdio.h>
#include <stdlib.h>

#include <vector>
#include <string>

#include "bamliquidator.h"

/* The MIT License (MIT) 

   Copyright (c) 2013 Xin Zhong and Charles Lin

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

int parseArgs(std::string& bamfile, std::string& chromosome, 
              unsigned int& start, unsigned int& stop,
              char& strand, unsigned int& spnum,
              unsigned int& extendlen,
              const int argc, char* argv[])
{
  if(argc!=8)
  {
    printf("[ bamliquidator ] output to stdout\n1. bam file (.bai file has to be at same location)\n2. chromosome\n3. start\n4. stop\n5. strand +/-, use dot (.) for both strands\n6. number of summary points\n7. extension length\n\nNote that each summary point is floor((stop-start)/(number of summary points)) long,\nand if it doesn't divide evenly then the range is truncated.\n");
    return 1;
  }

  bamfile=argv[1];
  chromosome=argv[2];

  char* tail=NULL;
  start=strtol(argv[3],&tail,10);
  if(tail[0]!='\0')
  {
    fprintf(stderr, "wrong start (%s)\n", argv[3]);
    return 1;
  }
  stop=strtol(argv[4],&tail,10);
  if(tail[0]!='\0' || stop<=start)
  {
    fprintf(stderr, "wrong stop (%s)\n", argv[4]);
    return 1;
  }
  strand=argv[5][0];
  if(strand!='+' && strand!='-' && strand!='.')
  {
    fputs("wrong strand, must be +/-/.\n",stderr);
    return 1;
  }
  spnum=strtol(argv[6],&tail,10);
  if(tail[0]!='\0' || spnum<=0)
  {
    fprintf(stderr, "wrong spnum (%s)\n", argv[6]);
    return 1;
  }
  extendlen=(unsigned short)strtol(argv[7],&tail,10);
  if(tail[0]!='\0')
  {
    fprintf(stderr, "wrong extension length (%s)\n", argv[7]);
    return 1;
  }

  return 0;
}

int main(int argc, char* argv[])
{
  std::string bamfile;
  std::string chromosome;
  unsigned int start = 0;
  unsigned int stop  = 0;
  char strand = 0;
  unsigned int spnum = 0;
  unsigned int extendlen = 0;
  if (parseArgs(bamfile, chromosome, start, stop, strand, spnum, extendlen, argc, argv) != 0)
  {
    return 1;
  }

  const std::vector<double> counts = liquidate(bamfile, chromosome, start, stop,
    strand, spnum, extendlen);

  for(double count : counts)
  {
    printf("%d\n", (int) count);
  }

  return 0;
}
