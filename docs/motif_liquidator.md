## Motif Liquidator: High Performance DNA Motif Analysis

* [Overview](https://github.com/BradnerLab/pipeline/wiki/motif_liquidator#motif-liquidator-high-performance-dna-motif-analysis)
    * [Collaborators](https://github.com/BradnerLab/pipeline/wiki/motif_liquidator#collaborators)
* [Quick Start](https://github.com/BradnerLab/pipeline/wiki/motif_liquidator#quick-start)
    * [Download](https://github.com/BradnerLab/pipeline/wiki/motif_liquidator#1-download)
* [Comparison to FIMO](https://github.com/BradnerLab/pipeline/wiki/motif_liquidator#comparison-to-fimo)
* [Command Line Usage](https://github.com/BradnerLab/pipeline/wiki/motif_liquidator#command-line-usage)
* [Developer Notes](https://github.com/BradnerLab/pipeline/wiki/motif_liquidator#developer-notes)
    * [Compile from Source](https://github.com/BradnerLab/pipeline/wiki/motif_liquidator#compile-from-source)
    * [Creating Binaries For Distribution](https://github.com/BradnerLab/pipeline/wiki/motif_liquidator#creating-binaries-for-distribution)
* [Reproducible Performance Tests](https://github.com/BradnerLab/pipeline/wiki/motif_liquidator#reproducible-performance-tests)

## Overview

* Find/score hundreds of [motifs](https://en.wikipedia.org/wiki/Sequence_motif) across millions of DNA reads in seconds
* Up to 160x faster than [MEME FIMO](http://meme-suite.org/doc/fimo.html)
    * uses same [MEME motif format](http://meme-suite.org/doc/meme-format.html)
    * motif_liquidator and FIMO results are identical
* Multithreaded and optimized C++ to scale
* Supports DNA sequences in FASTA file format, with beta support for the BAM file format

### Collaborators

Dr. Alexander Federation originated the idea to develop motif liquidator.  John DiMatteo is the lead developer.  Dr. Charles Lin and Adam Torgerson contribute and guide development.

* John DiMatteo, Senior Software Engineer at [Boulder Labs](http://www.boulderlabs.com/)
* [Charles Y. Lin, Ph.D.](https://www.charleslinlab.org/)
    * Assistant Professor, [Department of Molecular and Human Genetics](https://www.bcm.edu/departments/molecular-and-human-genetics/), [Dan L. Duncan Cancer Center](https://www.bcm.edu/centers/cancer-center), Baylor College of Medicine
* [Alexander Federation, Ph.D.](http://www.alexfederation.com/)
    * Postdoctoral Fellow, [Altius Institute for Biomedical Sciences](http://www.stamlab.org/), Depts. of Genome Science and Medicine, University of Washington
* Adam Torgerson, Co-owner of [Boulder Labs](http://www.boulderlabs.com/)

## Quick Start

### 1. Download

Download one of the following binaries:
* RedHat/CentOS: [liquidator_1.8_centos_6_8.tar.gz](https://www.dropbox.com/sh/342j6x5cx2r1vr8/AAAdqAB--xjgwgYcYPsBw-QBa/liquidator_1.8_centos_6_8.tar.gz)
* Ubuntu 16.04 Xenial [liquidator_1.8_ubuntu_16_04.tar.gz](https://www.dropbox.com/sh/342j6x5cx2r1vr8/AACWotjJDLb8ya1gSj2aQpMma/liquidator_1.8_ubuntu_16_04.tar.gz)
* Ubuntu 14.04 Trusty [liquidator_1.8_ubuntu_14_04.tar.gz](https://www.dropbox.com/sh/342j6x5cx2r1vr8/AACCLH6HLOZ2DHy3SgBh_lfya/liquidator_1.8_ubuntu_14_04.tar.gz)

If you aren't sure which to download, check `cat /etc/issue`. If you have any problems getting started (e.g. need additional binaries, unit tests fail, etc.), please email John DiMatteo at jdimatteo@boulderlabs.com

This wiki will be updated with links to newer binaries as they become available.

### 2. Verify Tests Pass

For example:

```
$ wget https://www.dropbox.com/sh/342j6x5cx2r1vr8/AACCLH6HLOZ2DHy3SgBh_lfya/liquidator_1.8_ubuntu_14_04.tar.gz
$ tar -xzf liquidator_1.8_ubuntu_14_04.tar.gz
$ cd liquidator
liquidator$ source source_this_from_this_directory.sh
liquidator$ ./cpp_test
...
[  PASSED  ]
```

Verify the last line of cpp_test says "PASSED".

### 3. Run motif_liquidator

For example:

```
$ wget https://www.dropbox.com/s/83gi7y5c2yzbsoz/04032013_D1L57ACXX_4.TTAGGC.hg18.bwt.chr1.fasta
$ wget https://www.dropbox.com/sh/342j6x5cx2r1vr8/AADLE4g39BOcjCkhQWQPquzHa/HOCOMOCOv9_AHR_si.meme
$ motif_liquidator HOCOMOCOv9_AHR_si.meme 04032013_D1L57ACXX_4.TTAGGC.hg18.bwt.chr1.fasta -o motifliquidator.out
$ head -5 motifliquidator.out 
#pattern name	sequence name	start	stop	strand	score	p-value	q-value	matched sequence
AHR_si	WIGTC-HISEQ:4:2106:13465:96775#TTAGGC/1;1	7	15	+	11.4605	8.03e-05		ATTGCGTGT
AHR_si	WIGTC-HISEQ:4:2203:18396:34488#TTAGGC/1;1	18	26	+	11.7763	5.38e-05		CTCGCGTGC
AHR_si	WIGTC-HISEQ:4:1115:19032:44760#TTAGGC/1;0	20	28	-	11.75	5.98e-05		TGTGCGTGC
AHR_si	WIGTC-HISEQ:4:1305:6766:23451#TTAGGC/1;1	20	28	-	11.75	5.98e-05		TGTGCGTGC
```

Note that if you don't otherwise have the prerequisite shared libraries installed, you will need to `source source_this_from_this_directory.sh` as shown in the "Verify Tests Pass" section.

## Comparison to FIMO

In all testing to date, motif_liquidator produces exactly the same results as FIMO, just faster and in different line orders.

For example, following example shows motif_liquidator producing exact same results:

```
$ wget https://www.dropbox.com/s/83gi7y5c2yzbsoz/04032013_D1L57ACXX_4.TTAGGC.hg18.bwt.chr1.fasta
$ wget https://www.dropbox.com/sh/342j6x5cx2r1vr8/AADLE4g39BOcjCkhQWQPquzHa/HOCOMOCOv9_AHR_si.meme
$ motif_liquidator HOCOMOCOv9_AHR_si.meme 04032013_D1L57ACXX_4.TTAGGC.hg18.bwt.chr1.fasta -o motifliquidator_HOCOMOCOv9_AHR_si.out

$ fimo --motif AHR_si -verbosity 1 -text -oc FIMO HOCOMOCOv9_AHR_si.meme 04032013_D1L57ACXX_4.TTAGGC.hg18.bwt.chr1.fasta > fimo_HOCOMOCOv9_AHR_si.out
Warning: text mode turns off computation of q-values

$ diff -s <(sort fimo_HOCOMOCOv9_AHR_si.out) <(sort motifliquidator_HOCOMOCOv9_AHR_si.out)
Files /dev/fd/63 and /dev/fd/62 are identical
```

## Command Line Usage

```
$ motif_liquidator --help
usage: motif_liquidator [options] motif fasta|bam
version 1.8.0-030671e

Score DNA sequences against a motif, with exact same scores as MEME's FIMO.
For more info, see https://github.com/BradnerLab/pipeline/wiki/motif_liquidator

positional arguments:
  motif                   A MEME style position weight matrix file.
  fasta|bam               A fasta file or (sorted and indexed) bam file.

optional arguments:
  -m [ --motif-name ] arg Motif name (default is all motifs in the motif file)
  -b [ --background ] arg MEME style background frequency file.  Note that only
                          0-order background (single nucleotide) frequenceis 
                          are currently used (just like FIMO).  Backgrounds 
                          specified in the motif file are never used (just like
                          default FIMO behavior).
  -h [ --help ]           Display this help and exit.
  -o [ --output ] arg     File to write matches to. Output is fimo style for 
                          fasta input, and output is a (sorted/indexed) .bam 
                          for bam input.
  -r [ --region ] arg     .bed or .gff region file for filtering bam input.
  -u [ --unmapped-only ]  Only scores unmapped reads from bam.
  -p [ --print ] arg      For bams, additionally prints detailed fimo style 
                          output to stdout.  Specify '-p fimo' for fimo style 
                          output or '-p mapped-fimo' for the sequence name to 
                          include the chromosome and the for start/stop values 
                          to be the mapped positions.
```


## Developer Notes

### Compile from Source

For Ubuntu 14.04 and 16.04:

1. `sudo apt-get install git libbam-dev libboost-dev libboost-timer-dev libgoogle-perftools-dev libtbb-dev samtools cmake libz-dev libhdf5-serial-dev libboost-program-options-dev libboost-filesystem-dev libboost-system-dev unzip g++`
2. `git clone -b motif_liquidator_matrix https://github.com/BradnerLab/pipeline.git`
3. `cd pipeline/liquidator`
4. `make cpp_test && make`
5. `./cpp_test`
6. verify last line of cpp_test output says `[  PASSED  ]`

The [bamliquidator build instructions](https://github.com/BradnerLab/pipeline/wiki/bamliquidator#developers) might be helpful for building on RedHat/CentOS and MacOS.

### Creating Binaries For Distribution

First build from source, then copy all the used libraries into a zip, or just update a previously packaged zip.

E.g. 

```
liquidator$ cat /etc/issue
Ubuntu 14.04.5 LTS \n \l
liquidator$ mkdir liquidator
liquidator$ ldd motif_liquidator
	linux-vdso.so.1 =>  (0x00007ffc4c42c000)
	libz.so.1 => /lib/x86_64-linux-gnu/libz.so.1 (0x00007f1bc4803000)
	libpthread.so.0 => /lib/x86_64-linux-gnu/libpthread.so.0 (0x00007f1bc45e5000)
	libtbb.so.2 => /usr/lib/libtbb.so.2 (0x00007f1bc43b0000)
	libboost_program_options.so.1.54.0 => /usr/lib/x86_64-linux-gnu/libboost_program_options.so.1.54.0 (0x00007f1bc4142000)
	libboost_filesystem.so.1.54.0 => /usr/lib/x86_64-linux-gnu/libboost_filesystem.so.1.54.0 (0x00007f1bc3f2c000)
	libboost_system.so.1.54.0 => /usr/lib/x86_64-linux-gnu/libboost_system.so.1.54.0 (0x00007f1bc3d27000)
	libstdc++.so.6 => /usr/lib/x86_64-linux-gnu/libstdc++.so.6 (0x00007f1bc3a23000)
	libm.so.6 => /lib/x86_64-linux-gnu/libm.so.6 (0x00007f1bc371d000)
	libgcc_s.so.1 => /lib/x86_64-linux-gnu/libgcc_s.so.1 (0x00007f1bc3506000)
	libc.so.6 => /lib/x86_64-linux-gnu/libc.so.6 (0x00007f1bc313e000)
	/lib64/ld-linux-x86-64.so.2 (0x000055ca23719000)
jdimatteo@bb8:~/dev/motifliquidator/liquidator$ ldd motif_liquidator
	linux-vdso.so.1 =>  (0x00007fff658e8000)
	libz.so.1 => /lib/x86_64-linux-gnu/libz.so.1 (0x00007fca772c5000)
	libpthread.so.0 => /lib/x86_64-linux-gnu/libpthread.so.0 (0x00007fca770a7000)
	libtbb.so.2 => /usr/lib/libtbb.so.2 (0x00007fca76e72000)
	libboost_program_options.so.1.54.0 => /usr/lib/x86_64-linux-gnu/libboost_program_options.so.1.54.0 (0x00007fca76c04000)
	libboost_filesystem.so.1.54.0 => /usr/lib/x86_64-linux-gnu/libboost_filesystem.so.1.54.0 (0x00007fca769ee000)
	libboost_system.so.1.54.0 => /usr/lib/x86_64-linux-gnu/libboost_system.so.1.54.0 (0x00007fca767e9000)
	libstdc++.so.6 => /usr/lib/x86_64-linux-gnu/libstdc++.so.6 (0x00007fca764e5000)
	libm.so.6 => /lib/x86_64-linux-gnu/libm.so.6 (0x00007fca761df000)
	libgcc_s.so.1 => /lib/x86_64-linux-gnu/libgcc_s.so.1 (0x00007fca75fc8000)
	libc.so.6 => /lib/x86_64-linux-gnu/libc.so.6 (0x00007fca75c00000)
	/lib64/ld-linux-x86-64.so.2 (0x000056160a5d0000)
liquidator$ mkdir liquidator/liquidator_dependencies_20170321
liquidator$ cp /lib/x86_64-linux-gnu/libz.so.1 /lib/x86_64-linux-gnu/libpthread.so.0 /usr/lib/libtbb.so.2 /usr/lib/x86_64-linux-gnu/libboost_program_options.so.1.54.0 /usr/lib/x86_64-linux-gnu/libboost_filesystem.so.1.54.0 /usr/lib/x86_64-linux-gnu/libboost_system.so.1.54.0 /usr/lib/x86_64-linux-gnu/libstdc++.so.6 /lib/x86_64-linux-gnu/libm.so.6 /lib/x86_64-linux-gnu/libgcc_s.so.1 /lib/x86_64-linux-gnu/libc.so.6 /lib64/ld-linux-x86-64.so.2 liquidator/liquidator_dependencies_20170321
liquidator$ cp bamliquidator bamliquidator_bins bamliquidator_regions cpp_test motif_liquidator liquidator
liquidator$ vim liquidator/source_this_from_this_directory.sh # copy from a prior build
liquidator$ tar -czvf liquidator_1.6_ubuntu_14_04.tar.gz liquidator
# share on dropbox or something
```

## Reproducible Performance Tests

Measurements were done on AWS EC2 instances
* [c4.8xlarge](https://aws.amazon.com/ec2/instance-types/): 18 cores (36 hyper threads), high frequency Intel Xeon E5-2666 v3 (Haswell) processors, 60 GB memory
* Ubuntu Server 16.04 LTS (HVM), SSD Volume Type - ami-f4cc1de2

Note that the tests below used Provisioned IOPs SSD (IO1) with 100 GiB and 5000 IOPS, which is not recommended due to the high cost.  AWS General Purpose SSD should be plenty fast, since these tests are not IO bound and warm runs will have all files cached in RAM anyway.  Running motif liquidator first or just doing a `cat` of each input file is sufficient to do a warm FIMO run without waiting an additional 4 hours for FIMO to run a second time.

### setup environment

```
sudo apt-get install git libbam-dev libboost-dev libboost-timer-dev libgoogle-perftools-dev libtbb-dev samtools cmake libz-dev libhdf5-serial-dev libboost-program-options-dev libboost-filesystem-dev libboost-system-dev unzip g++

git clone -b motif_liquidator_matrix https://github.com/BradnerLab/pipeline.git
cd pipeline/liquidator
make cpp_test && make
./cpp_test # passed
./motif_liquidator --version # version 1.8.0-030671e
git log -1 # 3/28/2017 commit 030671e05178180b991bcca2e8abeb59d7c7cded

wget http://meme-suite.org/meme-software/4.11.3/meme_4.11.3_1.tar.gz
tar zxf meme_4.11.3_1.tar.gz
cd meme_4.11.3
./configure --prefix=$HOME/meme --with-url=http://meme-suite.org --enable-build-libxml2 --enable-build-libxslt
make
make test # note that several other meme tools failed but all fimo labeled tests reported passing
make install
export PATH=$HOME/meme/bin:$PATH

wget https://www.dropbox.com/s/83gi7y5c2yzbsoz/04032013_D1L57ACXX_4.TTAGGC.hg18.bwt.chr1.fasta
wget http://meme-suite.org/meme-software/Databases/motifs/motif_databases.12.15.tgz
tar -zxf motif_databases.12.15.tgz
```

### measurements

Cold run (run after rebooting machine):

```
/usr/bin/time -v pipeline/liquidator/motif_liquidator motif_databases/HUMAN/HOCOMOCOv9.meme 04032013_D1L57ACXX_4.TTAGGC.hg18.bwt.chr1.fasta -o motifliquidator.out
	Command being timed: "pipeline/liquidator/motif_liquidator motif_databases/HUMAN/HOCOMOCOv9.meme 04032013_D1L57ACXX_4.TTAGGC.hg18.bwt.chr1.fasta -o motifliquidator.out"
	User time (seconds): 3204.90
	System time (seconds): 25.68
	Percent of CPU this job got: 3420%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 1:34.45
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 1388688
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 1
	Minor (reclaiming a frame) page faults: 654364
	Voluntary context switches: 4506
	Involuntary context switches: 6126
	Swaps: 0
	File system inputs: 436168
	File system outputs: 1869904
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0
```

Warm run 1 (ran right after cold run):

```
/usr/bin/time -v pipeline/liquidator/motif_liquidator motif_databases/HUMAN/HOCOMOCOv9.meme 04032013_D1L57ACXX_4.TTAGGC.hg18.bwt.chr1.fasta -o motifliquidator_warm.out
        Command being timed: "pipeline/liquidator/motif_liquidator motif_databases/HUMAN/HOCOMOCOv9.meme 04032013_D1L57ACXX_4.TTAGGC.hg18.bwt.chr1.fasta -o motifliquidator_warm.out"
	User time (seconds): 3212.70
	System time (seconds): 31.52
	Percent of CPU this job got: 3540%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 1:31.62
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 1404136
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 588163
	Voluntary context switches: 1238
	Involuntary context switches: 6089
	Swaps: 0
	File system inputs: 0
	File system outputs: 1869904
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0
```

Partly warmed run (ran after motif_liquidator so input files already cached in memory):

```
/usr/bin/time -v fimo -verbosity 1 -text -oc FIMO motif_databases/HUMAN/HOCOMOCOv9.meme 04032013_D1L57ACXX_4.TTAGGC.hg18.bwt.chr1.fasta > fimo.out                                 
Warning: text mode turns off computation of q-values
        Command being timed: "fimo -verbosity 1 -text -oc FIMO motif_databases/HUMAN/HOCOMOCOv9.meme 04032013_D1L57ACXX_4.TTAGGC.hg18.bwt.chr1.fasta"
        User time (seconds): 14007.83
        System time (seconds): 811.88
        Percent of CPU this job got: 99%
        Elapsed (wall clock) time (h:mm:ss or m:ss): 4:07:03
        Average shared text size (kbytes): 0
        Average unshared data size (kbytes): 0
        Average stack size (kbytes): 0
        Average total size (kbytes): 0
        Maximum resident set size (kbytes): 641932
        Average resident set size (kbytes): 0
        Major (requiring I/O) page faults: 0
        Minor (reclaiming a frame) page faults: 29749862
        Voluntary context switches: 5074
        Involuntary context switches: 18926
        Swaps: 0
        File system inputs: 0
        File system outputs: 1869896
        Socket messages sent: 0
        Socket messages received: 0
        Signals delivered: 0
        Page size (bytes): 4096
        Exit status: 0
```

Correctness test:
```
sort fimo.out > fimo.out.sorted
sort motifliquidator.out > motifliquidator.out.sorted
diff -s fimo.out.sorted motifliquidator.out.sorted
Files fimo.out.sorted and motifliquidator.out.sorted are identical
```

Speedup Calculation:
```
(fimo wall seconds)/(motif liquidator wall seconds)
(4*60*60+7*60+3)/(1*60+31.62)
161.7878x faster in wall clock time
```