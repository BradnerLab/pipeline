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

VERSION := $(shell grep "__version__ =" bamliquidatorbatch/bamliquidator_batch.py | cut -d' ' -f 3 | sed "s/'//g")

# The directory to install in:
prefix = /usr/local
bindir = $(prefix)/bin

# I prefer clang++, but g++ is more easily available, so using that instead
#CC=clang++
CC=g++

# CPPFLAGS used to include march=native, but was removed so executables would 
# work on platforms besides the one used for the build
CPPFLAGS := -std=c++0x -O3 -g -Wall -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE
LDFLAGS := -O3 -g -Wall
LDLIBS := -lbam -lz -lpthread
ADDITIONAL_LDLIBS := -lhdf5 -lhdf5_hl -ltcmalloc_minimal -ltbb

GIT_DIFF_RC := $(shell git diff --quiet; echo $$?)
ifeq ($(GIT_DIFF_RC), 0)
CLEAN := clean
else
CLEAN := unclean
endif

# if someone does a ppa and/or pip upload, first change the name/email and commit it
UPLOADER := John DiMatteo
UPLOADER_EMAIL := jdimatteo@gmail.com

define CHANGELOG 
bamliquidator ($(VERSION)-0ppa1~%s) %s; urgency=low

  * Auto generated from makefile
  * $(shell git config --get remote.origin.url)
  * $(shell git rev-parse HEAD) ($(CLEAN)) 

 -- $(UPLOADER) <$(UPLOADER_EMAIL)>  $(shell date +"%a, %d %b %G %H:%M:%S %z")

endef
export CHANGELOG 

define SETUP_PY
from setuptools import setup

setup(
    name='BamLiquidatorBatch',
    version='$(VERSION)',
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
        "openpyxl >=1.6.1, <2.0.0",
        "bokeh == 0.4.4"
    ],
)

endef
export SETUP_PY

all: bamliquidator bamliquidator_bins bamliquidator_regions 

bamliquidator: bamliquidator.m.o bamliquidator.o
	$(CC) $(LDFLAGS) -o bamliquidator bamliquidator.o bamliquidator.m.o $(LDLIBS) 

bamliquidator_bins: bamliquidator_bins.m.o bamliquidator.o bamliquidator_logger.o
	$(CC) $(LDFLAGS) -o bamliquidator_bins bamliquidator.o bamliquidator_bins.m.o bamliquidator_logger.o \
					$(LDLIBS) $(ADDITIONAL_LDLIBS)

bamliquidator_regions: bamliquidator_regions.m.o bamliquidator.o bamliquidator_logger.o
	$(CC) $(LDFLAGS) -o bamliquidator_regions bamliquidator.o bamliquidator_regions.m.o bamliquidator_logger.o \
					$(LDLIBS) $(ADDITIONAL_LDLIBS) 

bamliquidator.m.o: bamliquidator.m.cpp
	$(CC) $(CPPFLAGS) -c bamliquidator.m.cpp

bamliquidator_bins.m.o: bamliquidator_bins.m.cpp
	$(CC) $(CPPFLAGS) -c bamliquidator_bins.m.cpp

bamliquidator_regions.m.o: bamliquidator_regions.m.cpp
	$(CC) $(CPPFLAGS) -c bamliquidator_regions.m.cpp
  
bamliquidator.o: bamliquidator.cpp bamliquidator.h
	$(CC) $(CPPFLAGS) -pthread -c bamliquidator.cpp

bamliquidator_logger.o: bamliquidator_logger.cpp bamliquidator_logger.h
	$(CC) $(CPPFLAGS) -c bamliquidator_logger.cpp

EXECUTABLES = bamliquidator bamliquidator_bins bamliquidator_regions

archive:
	mkdir bamliquidator-$(VERSION)
	cp *.h *.cpp makefile bamliquidator-$(VERSION)
	tar -czf bamliquidator-$(VERSION).tar.gz bamliquidator-$(VERSION)
	rm -rf bamliquidator-$(VERSION)

deb: clean archive
	mv bamliquidator-$(VERSION).tar.gz bamliquidator_$(VERSION).orig.tar.gz 
	tar xf bamliquidator_$(VERSION).orig.tar.gz
	cp -R debian bamliquidator-$(VERSION)
	cp -R bamliquidator-$(VERSION) bamliquidator-$(VERSION)-precise
	mv bamliquidator-$(VERSION) bamliquidator-$(VERSION)-trusty
	printf "$$CHANGELOG" precise precise > bamliquidator-$(VERSION)-precise/debian/changelog 
	printf "$$CHANGELOG" trusty trusty > bamliquidator-$(VERSION)-trusty/debian/changelog 
	cd bamliquidator-$(VERSION)-precise/debian; debuild -S
	cd bamliquidator-$(VERSION)-trusty/debian; debuild -S

dput: deb
	dput ppa:bradner-computation/pipeline bamliquidator_$(VERSION)-*precise_source.changes
	dput ppa:bradner-computation/pipeline bamliquidator_$(VERSION)-*trusty_source.changes

pip:
	echo "$$SETUP_PY" > setup.py
	python setup.py sdist upload

clean:
	rm -f $(EXECUTABLES) *.o MANIFEST setup.py bamliquidator*.tar.gz
	rm -rf bamliquidator*precise* bamliquidator*trusty* BamLiquidatorBatch.egg-info dist

install: all
	install $(EXECUTABLES) $(DESTDIR)$(bindir)
