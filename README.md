pipeline
========

bradner lab computation pipeline scripts

[![Build Status](https://travis-ci.org/BradnerLab/pipeline.svg)](https://travis-ci.org/BradnerLab/pipeline)

For more documentation, see the [wiki](https://github.com/BradnerLab/pipeline/wiki/bamliquidator)
associated with this repository.

## Containers

You can use the software here easily using containers, either Docker or Singularity.
Instructions for building and usage are provided below.

### Build

A [Dockerfile](bamliquidator_internal/Dockerfile) is provided for you to build the container locally,
and then (optionally) push to a container registry to pull as Singularity.
Let's walk through steps to build the container:

```bash
cd bamliquidator_internal
docker build -t bioliquidator/bamliquidator .
```

You can then push this to a registry, or use [Singularity](https://sylabs.io/guides/latest/user-guide/) to
build a container from your local docker daemon.

### Option 1. Push to a Registry and Pull with Singularity

```bash
docker push bioliquidator/bamliquidator
singularity pull docker://bioliquidator/bamliquidator
```

## Option 2: Build Singularity from local Docker

```bash
singularity build bamliquidator_latest.sif docker-daemon://bioliquidator/bamliquidator:latest
```

### Docker Usage

With the docker container, the entrypoint provides access to the executable:

```
docker run -it bioliquidator/bamliquidator
usage: bamliquidator_batch.py [-h] [-b BIN_SIZE | -r REGIONS_FILE]
                              [-o OUTPUT_DIRECTORY] [-c COUNTS_FILE] [-f]
                              [-e EXTENSION] [--sense {+,-,.}] [-m]
                              [--region_format {gff,bed}] [--skip_plot]
                              [--black_list BLACK_LIST [BLACK_LIST ...]] [-q]
                              [-n NUMBER_OF_THREADS] [--xml_timings]
                              [--version]
                              bam_file_path
bamliquidator_batch.py: error: the following arguments are required: bam_file_path
```

You'll need to provide input files (and paths) that are bound relative to the container.
For example, if I create folders "inputs" and "output" to write data, I'd need to bind both
to the container as volumes, and reference my input files as if they are in the container.

```bash
$ docker run -it -v $PWD/input:/input $PWD/output:/output bioliquidator/bamliquidator /input/04032013_D1L57ACXX_4.TTAGGC.hg18.bwt.sorted.bam -o /output
```

Singularity usage is more seamless, discussed next.

### Singularity Usage

For this example, let's download some actual data first to test. Note that
we've followed the [build](#build) instructions above and have a container
`bamliquidator_latest.sif` in the present working directory. First, let's download some data.

```bash
# Get input files
mkdir -p input output
cd input
wget https://www.dropbox.com/s/bu75ojqr2ibkf57/04032013_D1L57ACXX_4.TTAGGC.hg18.bwt.sorted.bam
wget https://www.dropbox.com/s/a71ngagu2k8pgiv/04032013_D1L57ACXX_4.TTAGGC.hg18.bwt.sorted.bam.bai
cd ..
```

Now we can run the container and provide the same input and output arguments.

```bash
singularity run bamliquidator_latest.sif ./input/04032013_D1L57ACXX_4.TTAGGC.hg18.bwt.sorted.bam -o ./output
```

This would be functionally equivalent to:

```bash
singularity exec bamliquidator_latest.sif /opt/liquidator/bamliquidatorbatch/bamliquidator_batch.py ./input/04032013_D1L57ACXX_4.TTAGGC.hg18.bwt.sorted.bam -o ./output
Liquidating ./input/04032013_D1L57ACXX_4.TTAGGC.hg18.bwt.sorted.bam (file 1 of 1)
Liquidation completed: 6.141825 seconds, 29059326 reads, 4.731383 millions of reads per second
Cell Types: input
Normalizing and calculating percentiles for cell type input
Indexing normalized counts
ERROR	Skipping plotting because plots require a compatible version of bokeh -- see https://github.com/BradnerLab/pipeline/wiki/bamliquidator#Install . Bokeh module not found; consider running the following command to install:
sudo pip install bokeh==0.9.3 "openpyxl>=1.6.1,<2.0.0"
Summarizing
Post liquidation processing took 0.654222 seconds
```

And note that we can use relative paths because by default, Singularity binds the $PWD. This
is different than Docker, which enforces a totally isolated environment.


Copyright (c) 2013 Charles Lin and collaborators:
  - John DiMatteo
  - Nick Semenkovich
  - Xin Zhong

License: MIT (see LICENSE for details)
