#include "bamliquidator.h"

#include <stdio.h>
#include <samtools/sam.h>

#include <deque>
#include <stdexcept>
#include <sstream>

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

struct ReadItem
{
  unsigned int start;
  /* read stop is start + strlen(seq)
  this *stop* will only be used for computing density
  will not be reported to js for bed plotting
  the actual stop need to be determined by cigar
  */
  unsigned int stop;
  uint32_t flag; // flag from bam
  char strand;
  std::vector<uint32_t> cigar;
};


int intMin(int a, int b)
{
  if(a < b) return a;
  return b;
}

int intMax(int a, int b)
{
  if(a > b) return a;
  return b;
}



struct UserData
{
  std::deque<ReadItem> readItems;
  char strand;
  unsigned int extendlen;
};




static int bam_fetch_func(const bam1_t* b,void* data)
{
  if (b->core.tid < 0) return 0;

  UserData *udata=(UserData *)data;

  const bam1_core_t* c = &b->core;

  char strand= (c->flag&BAM_FREVERSE)?'-':'+';
  if(udata->strand=='+')
  {
    if(strand!='+') return 0;
  }
  else if(udata->strand=='-')
  {
    if(strand!='-') return 0;
  }

  ReadItem r;
  r.strand=strand;

  uint32_t* cigar = bam1_cigar(b);

  // get read length
  int i, readlen;
  if (b->core.tid < 0) return 0;
  for (i = readlen = 0; i < c->n_cigar; ++i)
  {
    int op = cigar[i]&0xf;
    if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CREF_SKIP)
      readlen += cigar[i]>>4;
  }

  r.cigar = std::vector<uint32_t>(cigar, cigar+c->n_cigar);

  r.start=c->pos;
  r.stop=c->pos+readlen;

  //printf("%d\t%d\t%c\t", r.start, r.stop, strand);

  // extend
  if(udata->extendlen>0)
  {
    if(strand=='+')
    {
      r.stop+=udata->extendlen;
    }
    else
    {
      r.start=intMax(0,r.start-udata->extendlen);
    }
  }

  //printf("%d\t%d\n", r.start, r.stop);

  r.flag=c->flag;
  udata->readItems.push_back(r);
  return 0;
}

std::deque<ReadItem> bamQuery_region(const samfile_t* fp, const bam_index_t* idx, const std::string& coord, char strand, unsigned int extendlen)
{
  // will not fill chromidx
  int ref,beg,end;
  int rc = bam_parse_region(fp->header,coord.c_str(),&ref,&beg,&end);
  if (rc != 0)
  {
    std::stringstream error_msg;
    error_msg << "bam_parse_region failed with return code " << rc;
    throw std::runtime_error(error_msg.str());
  }
  if(ref<0)
  {
    return std::deque<ReadItem>();
  }
  UserData d;
  d.strand=strand;
  d.extendlen=extendlen;
  bam_fetch(fp->x.bam,idx,ref,beg,end,&d,bam_fetch_func);
  return d.readItems;
}


std::vector<double> liquidate(const std::string& bamfile, const std::string& chromosome,
                              const unsigned int start, const unsigned int stop,
                              const char strand, const unsigned int spnum,
                              const unsigned int extendlen)
{
	samfile_t* fp=NULL;
	fp=samopen(bamfile.c_str(),"rb",0);
	if(fp == NULL)
	{
		throw std::runtime_error("samopen() error with " + bamfile);
	}

  bam_index_t* bamidx=NULL;
  bamidx=bam_index_load(bamfile.c_str());
	if (bamidx == NULL)
	{
		throw std::runtime_error("bam_index_load() error with " + bamfile);
	}

	std::vector<double> counts = liquidate(fp, bamidx, chromosome, start, stop, strand, spnum, extendlen);

  bam_index_destroy(bamidx);
  samclose(fp);

	return counts;
}

std::vector<double> liquidate(const samfile_t* fp, const bam_index_t* bamidx,
															const std::string& chromosome,
                              const unsigned int start, const unsigned int stop,
                              const char strand, const unsigned int spnum,
                              const unsigned int extendlen)
{
  std::vector<double> data(spnum, 0);
  
  std::string coord;
  {
    std::stringstream ss;
    ss << chromosome << ':' << start << '-' << stop;
    coord = ss.str();
  }

  /* fetch bed items for a region and compute density
  only deal with coord, so use generic item
  */
  if (stop < start) throw std::runtime_error("liquidate called with stop < start");
  const unsigned pieceLength = (stop-start) / spnum;

  unsigned startArr[spnum], stopArr[spnum];
  for(unsigned int i=0; i<spnum; i++)
  {
    startArr[i] = start + pieceLength*i;
    stopArr[i] = start + pieceLength*(i+1);
  }

  std::deque<ReadItem> items = bamQuery_region(fp,bamidx,coord,strand,extendlen);

  for(const ReadItem& item : items)
  {
    // collapse this bed item onto the density counter
    for(unsigned int i=0; i<spnum; i++)
    {
      if(item.start > stopArr[i]) continue;
      if(item.stop < startArr[i]) break;
      int start=intMax(item.start,startArr[i]);
      int stop=intMin(item.stop,stopArr[i]);
      if(start<stop)
      {
        // as Charles suggested, add the fraction of the read (overlapping with the bin)
        // instead of just counting the read
        data[i] += stop-start;
      }
    }
  }

  return data;
}
