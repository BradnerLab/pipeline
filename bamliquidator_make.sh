# SAM_DIR points to the directory where the samtools program has been built
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

#SAM_DIR=/home/bradneradmin/samtools
SAM_DIR=/usr/local/Cellar/samtools/0.1.19/include/bam

gcc -O -g  -Wall -Wformat -Wimplicit -Wreturn-type -Wuninitialized -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_GNU_SOURCE -DCOLOR32 -Wall -Wformat -Wimplicit -Wreturn-type -I$SAM_DIR -o bamliquidator.o -c bamliquidator.c
gcc -O -g  -Wall -Wformat -Wimplicit -Wreturn-type -Wuninitialized -o bamliquidator bamliquidator.o -pthread -L$SAM_DIR -lbam -lz -lpthread -ldl 
