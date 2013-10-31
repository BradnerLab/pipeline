#/home/bradneradmin/samtools points to the directory where the samtools program has been built
#please change this to the appropriate directory
#samtools can be downloaded from here http://sourceforge.net/projects/samtools/files/latest/download?source=files
#samtools has two dependencies: ncurses http://www.gnu.org/software/ncurses/ 
#and zlib http://zlib.net/

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


gcc -O -g  -Wall -Wformat -Wimplicit -Wreturn-type -Wuninitialized -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_GNU_SOURCE -DCOLOR32 -Wall -Wformat -Wimplicit -Wreturn-type -I/home/bradneradmin/samtools -o bamliquidator.o -c bamliquidator.c
gcc -O -g  -Wall -Wformat -Wimplicit -Wreturn-type -Wuninitialized -o bamliquidator bamliquidator.o -pthread -L/home/bradneradmin/samtools -lbam -lz -lpthread -ldl 
