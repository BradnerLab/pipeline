#ifndef PIPELINE_BAMLIQUIDATORINTERNAL_BAMLIQUIDATOR_H
#define PIPELINE_BAMLIQUIDATORINTERNAL_BAMLIQUIDATOR_H

#include <samtools/sam.h>

#include <vector>
#include <string>

/** 
 * Count the number of reads in a chromosome between start and stop.  This function 
 * is thread safe.  This function throws if an error is encountered opening or parsing
 * the bamfile with samtools.
 *
 * @param bamfile   the path of the bamfile (which is opened readonly),
 *                  e.g. "../tmp/04032013_D1L57ACXX_4.TTAGGC.hg18.bwt.sorted.bam" 
 * @param chromsome the chromosome to count on, e.g. "chr1" or "chrX" 
 * @param start     the first base pair index to count on (inclusive), e.g. 0 to start 
                    at the beginning   
 * @param stop      the last base pair index to count on (exclusive), e.g. 247249719
 * @param strand    '+' for the forward strand, '-' for the reverse strand, and '.' for both
 * @param spnum     number of summary points, e.g. if 4 with start of 0 and stop of 99,
                    the returned vector will have four counts, the first for the range
                    [0, 24], the second for [25, 49], third for [50, 74], and the last for 
                    [75, 99] -- use 1 for a single summary point
 * @param extenlen  extension length, e.g. 0 is usually used
 * 
 * @return the read counts for the range [start, stop], split into spnum pieces
 */
// todo: why is this a vector of doubles instead of a vector of integers?
std::vector<double> liquidate(const std::string& bamfile, const std::string& chromosome,
                              unsigned int start, unsigned int stop,
                              char strand, unsigned int spnum,
                              unsigned int extendlen);

/** 
 * Same as above function, except this is not thread safe (as bamfile/bamidx cannot be
 * used simultaneously in different threads).  This variant should be preferred when
 * looping over many start/stop values for the same bamfile in a single thread, since
 * opening the file/index can take more time than the liquidation. 
 */
std::vector<double> liquidate(const samfile_t* bamfile, const bam_index_t* bamidx,
															const std::string& chromosome,
                              unsigned int start, unsigned int stop,
                              char strand, unsigned int spnum,
                              unsigned int extendlen);

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

#endif  // PIPELINE_BAMLIQUIDATORINTERNAL_BAMLIQUIDATOR_H
