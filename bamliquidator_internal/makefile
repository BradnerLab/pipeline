###############################################################################
# The MIT License (MIT)
#
# Copyright (c) 2013 Xin Zhong and Charles Lin
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
###############################################################################

# For dependencies, please see https://github.com/BradnerLab/pipeline/wiki/bamliquidator

CPPFLAGS := -std=c++11 -stdlib=libc++ -O -g -Wall -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_GNU_SOURCE -DCOLOR32 
LDFLAGS := -stdlib=libc++ -O -g -Wall 
LDLIBS := -lbam -lz -lpthread 

all: ../bamliquidator bamliquidator_batch bamliquidator_regions 

../bamliquidator: bamliquidator.m.o bamliquidator.o
	clang++ $(LDFLAGS) -o ../bamliquidator bamliquidator.o bamliquidator.m.o $(LDLIBS) 

# note batch has additional dependencies separate from LDLIBS
bamliquidator_batch: bamliquidator_batch.m.o
	clang++ $(LDFLAGS) -o bamliquidator_batch bamliquidator.o bamliquidator_batch.m.o \
					$(LDLIBS) -lhdf5 -lhdf5_hl -ltcmalloc_minimal

# todo: experiment with both tcmalloc and the tbb allocator
# note regions has additional dependencies separate from LDLIBS
bamliquidator_regions: bamliquidator_regions.m.o bamliquidator.o
	clang++ $(LDFLAGS) -o bamliquidator_regions bamliquidator.o bamliquidator_regions.m.o \
					$(LDLIBS) -lhdf5 -lhdf5_hl -ltbb -ltcmalloc_minimal

bamliquidator.m.o: bamliquidator.m.cpp
	clang++ $(CPPFLAGS) -c bamliquidator.m.cpp

bamliquidator_batch.m.o: bamliquidator_batch.m.cpp
	clang++ $(CPPFLAGS) -c bamliquidator_batch.m.cpp

bamliquidator_regions.m.o: bamliquidator_regions.m.cpp
	clang++ $(CPPFLAGS) -c bamliquidator_regions.m.cpp
  
bamliquidator.o: bamliquidator.cpp bamliquidator.h
	clang++ $(CPPFLAGS) -I$(SAM_DIR) -pthread -c bamliquidator.cpp

clean:
	rm -f ../bamliquidator bamliquidator bamliquidator_batch bamliquidator_regions *.o 
