#dynamicEnhancer_rank.R

#produces a pair of hockey sticks that are red/green labeled

# The MIT License (MIT)

# Copyright (c) 2013 Charles Lin

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

#enhancerFile ='mergeTest/EC_BRD4_CON_ROSE/HG18_EC_MERGED_SUPERS_-0_+0_0KB_STITCHED_ENHANCER_DELTA_ENHANCER_TO_GENE_100KB_RANK.txt'

#nSuper1 = 347
#nSuper2 = 271
#===========================================================
#===============READING IN ARGUMNETS========================
#===========================================================

args <- commandArgs()

print(args[3:7])

enhancerFile = args[3]

name1 = args[4]
name2 = args[5]

nSuper1 = as.numeric(args[6])
nSuper2 = as.numeric(args[7])




#===========================================================
#==================GETTING IN TABLES========================
#===========================================================

enhancerTable = read.delim(enhancerFile)

name1Supers = which(enhancerTable[,15]<=nSuper1)
name2Supers = which(enhancerTable[,16]<=nSuper2)

superRows = union(name1Supers,name2Supers)

conservedSupers = intersect(name1Supers,name2Supers)

plotName = gsub('RANK','RANK_PLOT',enhancerFile)
plotName = gsub('txt','png',plotName)


png(filename=plotName,width = 800,height = 800)

plot(enhancerTable[name1Supers,15],enhancerTable[name1Supers,16],col='green',pch=16,xlim = c(1,max(enhancerTable[,15:16])),ylim = c(1,max(enhancerTable[15:16])),log='xy',ylab=paste('Rank in',name2),xlab= paste('Rank in',name1))
points(enhancerTable[name2Supers,15],enhancerTable[name2Supers,16],col='red',pch=16)

points(enhancerTable[conservedSupers,15],enhancerTable[conservedSupers,16],col='grey',pch=16)
abline(h=nSuper2)
abline(v=nSuper1)
text(1,1,paste(length(conservedSupers),'conserved supers'),pos=4)

text(1,nSuper1+100,paste(nSuper1,name1,'\nonly supers'),pos=4)
text(nSuper2+100,1,paste(nSuper2,name2,'\nonly supers'),pos=4)

dev.off()