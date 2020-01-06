#!/usr/bin/env python3
#get region both left and right of position of interest.

import os
import string

pseudogff = open('annotation/genes.gff', 'w')
promotersgff = open('annotation/promoters.gff', 'w')
upstreamgff = open('annotation/upstream.gff', 'w')
downstreamgff = open('annotation/downstream.gff', 'w')
    
def parsingGeneLocations(CHR, results):
  flank = 2000 #flank region
  lines = results.split("\t")
  lines[3] = int(lines[3])
  lines[4] = int(lines[4])
  if lines[6] == "+":
    end = lines[3] + flank
    start = lines[3] - flank
    upend = lines[3] - 1
    upstart = lines[3] - flank
    downstart = lines[4] + 1
    downend = lines[4] + flank
  elif lines[6] == "-":
    end = lines[4] + flank
    start = lines[4] - flank
    upend = lines[4] + flank
    upstart = lines[4] + 1
    downend = lines[3] - 1
    downstart = lines[3] - flank
    
  if downstart < 1:
    downstart = 1
  if upstart < 1:
    upstart = 1
    
  if upend > int(CHR[lines[0]]):
    upend = CHR[lines[0]]
  if downend > int(CHR[lines[0]]):
    downend = CHR[lines[0]]
    
  promotersgff.write ("{0}\t{1}\t{2}\t{3}\n".format("\t".join(lines[0:3]), start, end, "\t".join(lines[5:])))
  upstreamgff.write ("{0}\t{1}\t{2}\t{3}\n".format("\t".join(lines[0:3]), upstart, upend, "\t".join(lines[5:])))
  downstreamgff.write ("{0}\t{1}\t{2}\t{3}\n".format("\t".join(lines[0:3]), downstart, downend, "\t".join(lines[5:])))

def main():
  from optparse import OptionParser
  usage = "usage: %prog -g [GTF/GFF file] -f [FEATURE TYPE (gene/transcript)] -c [CHROMSIZES]"
  parser = OptionParser(usage = usage)
  
  parser.add_option("-g","--gtf", dest="gtf",nargs = 1, default=None,
                      help = "Enter .gtf/gff file to be processed.")
  parser.add_option("-f","--feature", dest="feature",nargs = 1, default=None,
                      help = "Enter feature type [gene/transcript] to be processed.")
  parser.add_option("-c","--chrom", dest="chrom",nargs = 1, default=None,
                      help = "Enter ucsc chrom sizes file to be processed.")
  
  (options, args) = parser.parse_args()
  
  print(options)
  print(args)

  if options.chrom and options.gtf and options.feature:
    chromsizes = open(options.chrom, 'r')
    CHR = {}
    for line in chromsizes:
      line = line.split('\t')
      CHR[line[0]] = line[1]
    
    feature = options.feature
    if options.gtf.split('.')[-1] == 'gff':
      gffFile = open(options.gtf, 'r')
      for line in gffFile:
        if not line.startswith('#'):
          lines = line.split("\t")
          if lines[2] == feature:
            results = ("chr{0}\t{1}".format(lines[0],"\t".join(lines[1:])))
            pseudogff.write (results+"\n")
            parsingGeneLocations(CHR, results) 
    elif options.gtf.split('.')[-1] == 'gtf':
      gffFile = open(options.gtf, 'r')
      for line in gffFile:
        if not line.startswith('#'):
          lines = line.split("\t")
          if lines[2] == feature:
            newline = lines[8].split(' ')
            results = ("chr{0}\t{1}\t{2}={3}".format(lines[0], "\t".join(lines[1:8]), newline[0], newline[1]))
            pseudogff.write (results + "\n")
            parsingGeneLocations(CHR, results)
  else:
    parser.print_help()


if __name__ == "__main__":
    main()


