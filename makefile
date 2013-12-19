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

#SAM_DIR='./samtools/'
#SAM_DIR:=/home/bradneradmin/samtools
SAM_DIR:=/usr/local/Cellar/samtools/0.1.19/include/bam

# Please change SAM_DIR to the directory where the samtools program has been
# built.  bamliquidator needs some of the source files in the samtools
# directory.  note, this is not the location of the samtools program.
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


# this is a make file, so to build, just run make
# http://www.cprogramming.com/tutorial/makefiles.html

bamliquidator: bamliquidator.o
	clang++ -std=c++11 -O -g  -Wall -o bamliquidator bamliquidator.o -L$(SAM_DIR) -lbam -lz -ldl 
  
bamliquidator.o: bamliquidator.cpp
	clang++ -std=c++11 -O -g  -Wall -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_GNU_SOURCE -DCOLOR32 -I$(SAM_DIR) -o bamliquidator.o -pthread -c bamliquidator.cpp

clean:
	rm -f bamliquidator bamliquidator.o
