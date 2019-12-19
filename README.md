===== BAM2GFF: Genomic Region READ DENSITY distribution =====

#### === Changelog
1) USAGE

Executable anywhere as long as the PATHTO is correctly specified.
> BAM2GFF-call.sh ["GTF file"] ["feature type"] ["BAM file"] ["CHROM SIZES"] ["SAMPLENAME"]


2) DIRECTORY structure

    ├── LICENSE

    ├── README.md

    ├── BAM2GFF-call.sh    : bash wrapper

    ├── lib

        └── utils.py   : utilities method
    │   

    └── src
    
        ├── BAM2GFF_main.py    : calculates density of .bam reads in .gff regions
    
        ├── BAM2GFF_gtftogenes.py    : splits .gff into genomic regions (promoter, upstream, downstream, gene region)
    
        └── BAM2GFF_plots.R    : plots matplots and heatmaps of regions provided

    Total: 2 directories, 7 files


3) DEPENDENCIES

	* samtools
	* R version > 3.4
	* python2


=================================================================
