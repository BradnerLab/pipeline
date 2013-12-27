# The MIT License (MIT)

# Copyright (c) 2013 Xin Zhong and Charles Lin

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

###############################################################################

ifndef SAM_DIR
#SAM_DIR='./samtools/'
SAM_DIR:=/home/bradneradmin/samtools
endif

# Please change SAM_DIR to the directory where the samtools program has been
# built (or override it in a condition below).  bamliquidator needs some of the 
# source files in the samtools directory.  note, this is not the location of the
# samtools program.
#
# please change this to the appropriate directory
#
# samtools can be downloaded from here:
# 
#   http://sourceforge.net/projects/samtools/files/latest/download?source=files
#
# samtools has two dependencies:
#   1. ncurses http://www.gnu.org/software/ncurses/ 
#   2. zlib http://zlib.net/
#
# on a Mac, you can also install samtools via homebrew, e.g. homebrew install samtools

CPPFLAGS := -std=c++11 -O -g -Wall -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_GNU_SOURCE -DCOLOR32 
LDFLAGS := -O -g -Wall 
LDLIBS := -L$(SAM_DIR) -lbam -lz -ldl -lpthread -lboost_system 

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S), Darwin)
	SAM_DIR:=/usr/local/Cellar/samtools/0.1.19/include/bam
	CPPFLAGS += -stdlib=libc++ 
	LDFLAGS += -stdlib=libc++ 
endif

# this is a make file, so to build, just run make
# http://www.cprogramming.com/tutorial/makefiles.html


all: bamliquidator bamliquidate_batch 

bamliquidator: bamliquidator.m.o bamliquidator.o
	clang++ $(LDFLAGS) -o bamliquidator bamliquidator.o bamliquidator.m.o $(LDLIBS) 

bamliquidate_batch: bamliquidate_batch.m.o
	clang++ $(LDFLAGS) -o bamliquidate_batch bamliquidator.o bamliquidate_batch.m.o $(LDLIBS) 

bamliquidator.m.o: bamliquidator.m.cpp threadsafe_queue.h
	clang++ $(CPPFLAGS) -c bamliquidator.m.cpp

bamliquidate_batch.m.o: bamliquidate_batch.m.cpp
	clang++ $(CPPFLAGS) -c bamliquidate_batch.m.cpp
  
bamliquidator.o: bamliquidator.cpp
	clang++ $(CPPFLAGS) -I$(SAM_DIR) -pthread -c bamliquidator.cpp

clean:
	rm -f bamliquidator bamliquidate_batch bamliquidate_batch.m.o bamliquidator.o bamliquidator.m.o
