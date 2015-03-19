###############################################################################
# The MIT License (MIT)
#
# Copyright (c) 2013 John DiMatteo, Xin Zhong, and Charles Lin
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

VERSION := $(shell grep "version =" bamliquidatorbatch/common_util.py | cut -d' ' -f 3 | sed "s/'//g")
export VERSION

# The directory to install in:
prefix = /usr/local
bindir = $(prefix)/bin

# I prefer clang++, but g++ is more easily available, so using that instead
#CC=clang++
CC=g++

# CPPFLAGS used to include march=native, but was removed so executables would 
# work on platforms besides the one used for the build
CPPFLAGS := -std=c++11 -O3 -g -Wall
LDFLAGS := -O3 -g -Wall
LDLIBS := -lbam -lz -lpthread
ADDITIONAL_LDLIBS := -lhdf5 -lhdf5_hl -ltcmalloc_minimal -ltbb

# if someone else does a ppa dput upload, first change the name/email and commit it
UPLOADER := John DiMatteo
UPLOADER_EMAIL := jdimatteo@gmail.com
export UPLOADER
export UPLOADER_EMAIL

define SETUP_PY
from setuptools import setup

setup(
    name='BamLiquidatorBatch',
    version='$(VERSION)',
    description='Python wrapper around bamliquidator for analyzing the density of short DNA sequence read alignments in the BAM file format.',
    maintainer='$(UPLOADER)',
    maintainer_email='$(UPLOADER_EMAIL)',
    packages=['bamliquidatorbatch'],
    url='https://github.com/BradnerLab/pipeline/wiki/bamliquidator',
    license='The MIT License (MIT)',
    entry_points = {
        'console_scripts': [
            'bamliquidator_batch = bamliquidatorbatch.bamliquidator_batch:main',
            'bamliquidator_flattener = bamliquidatorbatch.flattener:main'
        ]
    },
    install_requires=[
        'numpy',
        'pandas',
	'redis',
	'tables',
	'numexpr'
    ]
)

endef
export SETUP_PY

all: bamliquidator bamliquidator_bins bamliquidator_regions motif_liquidator motif_global_liquidator motif_region_liquidator

bamliquidator: bamliquidator.m.o bamliquidator.o
	$(CC) $(LDFLAGS) -o bamliquidator bamliquidator.o bamliquidator.m.o $(LDLIBS) 

bamliquidator_bins: bamliquidator_bins.m.o bamliquidator.o bamliquidator_util.o
	$(CC) $(LDFLAGS) -o bamliquidator_bins bamliquidator.o bamliquidator_bins.m.o bamliquidator_util.o \
					$(LDLIBS) $(ADDITIONAL_LDLIBS)

bamliquidator_regions: bamliquidator_regions.m.o bamliquidator.o bamliquidator_util.o bamliquidator_regions.h
	$(CC) $(LDFLAGS) -o bamliquidator_regions bamliquidator.o bamliquidator_regions.m.o bamliquidator_util.o \
					$(LDLIBS) $(ADDITIONAL_LDLIBS) 

motif_liquidator: motif_liquidator.m.cpp
	$(CC) $(CPPFLAGS) motif_liquidator.m.cpp $(LDFLAGS) -o motif_liquidator $(LDLIBS)

motif_global_liquidator: motif_global_liquidator.m.cpp bamliquidator_util.o
	$(CC) $(CPPFLAGS) motif_global_liquidator.m.cpp bamliquidator_util.o $(LDFLAGS) -o motif_global_liquidator $(LDLIBS)

motif_region_liquidator: motif_region_liquidator.m.cpp bamliquidator_util.o bamliquidator_regions.h
	$(CC) $(CPPFLAGS) motif_region_liquidator.m.cpp bamliquidator_util.o $(LDFLAGS) -o motif_region_liquidator $(LDLIBS)

bamliquidator.m.o: bamliquidator.m.cpp
	$(CC) $(CPPFLAGS) -c bamliquidator.m.cpp

bamliquidator_bins.m.o: bamliquidator_bins.m.cpp
	$(CC) $(CPPFLAGS) -c bamliquidator_bins.m.cpp

bamliquidator_regions.m.o: bamliquidator_regions.m.cpp
	$(CC) $(CPPFLAGS) -c bamliquidator_regions.m.cpp
  
bamliquidator.o: bamliquidator.cpp bamliquidator.h
	$(CC) $(CPPFLAGS) -pthread -c bamliquidator.cpp

bamliquidator_util.o: bamliquidator_util.cpp bamliquidator_util.h
	$(CC) $(CPPFLAGS) -c bamliquidator_util.cpp

EXECUTABLES = bamliquidator bamliquidator_bins bamliquidator_regions motif_liquidator motif_global_liquidator

archive:
	mkdir bamliquidator-$(VERSION)
	cp *.h *.cpp makefile bamliquidator-$(VERSION)
	tar -czf bamliquidator-$(VERSION).tar.gz bamliquidator-$(VERSION)
	rm -rf bamliquidator-$(VERSION)
	mkdir bamliquidatorbatch_$(VERSION)
	cp -r bamliquidatorbatch bamliquidatorbatch_$(VERSION)
	echo "$$SETUP_PY" > bamliquidatorbatch_$(VERSION)/setup.py
	tar -czf bamliquidatorbatch_$(VERSION).orig.tar.gz bamliquidatorbatch_$(VERSION)

deb: clean archive
	./make_deb.sh

dput: clean archive
	debuild_args=-S ./make_deb.sh
	dput ppa:bradner-computation/pipeline deb_dist/bamliquidatorbatch_$(VERSION)-*precise_source.changes
	dput ppa:bradner-computation/pipeline deb_dist/bamliquidatorbatch_$(VERSION)-*trusty_source.changes
	dput ppa:bradner-computation/pipeline bamliquidator_$(VERSION)-*precise_source.changes
	dput ppa:bradner-computation/pipeline bamliquidator_$(VERSION)-*trusty_source.changes

gtest:
	wget https://googletest.googlecode.com/files/gtest-1.7.0.zip
	unzip gtest-1.7.0.zip
	mv gtest-1.7.0 gtest
	mkdir gtest/build
	(cd gtest/build; cmake ..; make)

cpp_test: gtest test.cpp bamliquidator_util.h
	$(CC)	$(CPPFLAGS) -I gtest/include gtest/build/libgtest.a test.cpp -o cpp_test

test: cpp_test 
	./cpp_test
	python bamliquidatorbatch/test.py

clean:
	rm -f $(EXECUTABLES) *.o MANIFEST setup.py bamliquidator*.tar.gz
	rm -rf bamliquidator*precise* bamliquidator*trusty* BamLiquidatorBatch.egg-info dist bamliquidatorbatch_* deb_dist gtest cpp_test

install: all
	install $(EXECUTABLES) $(DESTDIR)$(bindir)
