# Bamplot Turbo

* [Overview](#Overview)  
* [Install](#Install)  
* [Usage](#Usage)
* [Examples](#Examples)
    * [Plotting a single locus](#single_locus_example)

# Overview
* bamPlot_turbo is a set of scripts designed to create publication quality vector graphic plots of read density from nextgen sequencing data across specific genomic loci.
    * Read data must be in the form of sorted and indexed [BAM](http://samtools.sourceforge.net/) files
    * Input regions can either be specified individually e.g. `chr1:+:1-1000` or as a batch via a [.gff](http://genome.ucsc.edu/FAQ/FAQformat.html#format3) file.
    * Reference gene annotations are also plotted.  (Currently, only HG18,HG19,MM8,MM9 are supported)
    * Additional reference regions can be plotted when provided in [.bed](http://genome.ucsc.edu/FAQ/FAQformat.html#format1) format

# Install
* bamPlot_turbo is part of the bradnerLab pipeline module.  Install instructions under construction

# Usage
* bamPlot_turbo is called from the command line using python ./bamPlot_turbo.py

```bash
Usage: bamPlot_turbo.py [options] -g [GENOME] -b [SORTED BAMFILE(S)] -i [INPUTFILE] -o [OUTPUTFOLDER]

Options:
  -h, --help            show this help message and exit
  -b BAM, --bam=BAM     Enter a comma separated list of .bam files to be
                        processed.
  -i INPUT, --input=INPUT
                        Enter .gff or genomic region e.g. chr1:+:1-1000.
  -g GENOME, --genome=GENOME
                        specify a genome, HG18,HG19,MM8,MM9 are currently
                        supported
  -o OUTPUT, --output=OUTPUT
                        Enter the output folder.
  -c COLOR, --color=COLOR
                        Enter a colon separated list of colors e.g.
                        255,0,0:255,125,0, default samples the rainbow
  -s SENSE, --sense=SENSE
                        Map to '+','-' or 'both' strands. Default maps to
                        both.
  -e EXTENSION, --extension=EXTENSION
                        Extends reads by n bp. Default value is 200bp
  -r, --rpm             Normalizes density to reads per million (rpm) Default
                        is True
  -y YSCALE, --yScale=YSCALE
                        Choose either relative or uniform y axis scaling.
                        options = 'relative,uniform' Default is relative
                        scaling
  -n NAMES, --names=NAMES
                        Enter a comma separated list of names for your bams
  -p PLOT, --plot=PLOT  Choose either all lines on a single plot or multiple
                        plots. options = 'single,multiple'
  -t TITLE, --title=TITLE
                        Specify a title for the output plot(s), default will
                        be the coordinate region
  --save-temp           If flagged will save temporary files made by bamPlot
  --bed=BED             Add a comma separated list of bam files to plot
```

# Examples


#### Single locus plot
* plotting reads from two datasets at a single genomic region

```bash
python ./bamPlot_turbo.py -g 'hg18' -b $BAM1,$BAM2 -i 'chr1:+:100900000-100980000' -r -y 'UNIFORM' -t 'VCAM1' --bed $BED1 -o './'
```
