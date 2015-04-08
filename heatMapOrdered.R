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

library(graphics)


args <- commandArgs()

print(args[3:9])

referenceGFF = args[3]
mappedGFF = args[4]
colorVector = as.numeric(unlist(strsplit(args[5],',')))
color = rgb(colorVector[1],colorVector[2],colorVector[3],maxColorValue=255)
output = args[6]
geneListFile = args[7]
relative = as.numeric(args[8])
backgroundGFF = args[9]


#getting the reference order and color spectrum
referenceData <- read.delim(file=referenceGFF,sep="\t",header=TRUE)


if(nchar(geneListFile) > 5){
geneListTable = read.delim(geneListFile,header=FALSE)
geneList = as.vector(geneListTable[,1])}else{geneList = as.vector(seq(1,nrow(referenceData),1))}

#loading mappedGFF
mappedData <- read.delim(file=mappedGFF,sep="\t",header=TRUE)
mappedData <- as.matrix(mappedData[geneList,3:ncol(mappedData)])  #remove GENE_ID  & locusLine and force to matrix
colnames(mappedData) <- NULL

#performing background correction
if(backgroundGFF != 'NONE' | nchar(backgroundGFF) > 5){

  print('PERFORMING BACKGROUND CORRECTION')
  backgroundData <- read.delim(file=backgroundGFF,sep="\t",header=TRUE)	
  backgroundData <- as.matrix(backgroundData[geneList,3:ncol(backgroundData)])  #remove GENE_ID  & locusLine and force to matrix
  colnames(backgroundData) <- NULL
  mappedData = mappedData - backgroundData
  mappedData[which(mappedData < 0)] <- 0
}




referenceData <- as.matrix(referenceData[geneList,3:ncol(referenceData)])  #remove GENE_ID  & locusLine and force to matrix
colnames(referenceData) <- NULL
referenceOrder = order(apply(referenceData,1,mean,na.rm=TRUE))



#if scaling by reference
colorSpectrum <- colorRampPalette(c("white",color))(100)
minValue <- quantile(referenceData,na.rm=TRUE,prob=0.6,names=FALSE)
print('min value is')
print(minValue)
#maxValue <- quantile(referenceData,na.rm=TRUE,prob=0.95,names=FALSE)
maxValue = 3
print('max value is')
print(maxValue)
color_cuts <- seq(minValue,maxValue,length=100)

#add extreme points and one extra min color
#color_cuts <- c(min(referenceData,na.rm=TRUE), color_cuts,max(referenceData,na.rm=TRUE))
color_cuts <- c(min(referenceData,na.rm=TRUE), color_cuts,max(5,max(referenceData,na.rm=TRUE)))
colorSpectrum <- c(colorSpectrum[1],colorSpectrum)

#if relative scaling
if(relative==1){
  colorSpectrum <- colorRampPalette(c("white",color))(100)
  minValue <- quantile(mappedData,na.rm=TRUE,prob=0.6,names=FALSE)
  print('min value is')
  print(minValue)
  maxValue <- quantile(mappedData,na.rm=TRUE,prob=0.95,names=FALSE)
  print('max value is')
  print(maxValue)
  color_cuts <- seq(minValue,maxValue,length=100)

  #add extreme points and one extra min color
  color_cuts <- c(min(mappedData,na.rm=TRUE), color_cuts,max(mappedData,na.rm=TRUE))
  colorSpectrum <- c(colorSpectrum[1],colorSpectrum)
}


#reorder by reference
mappedData <- mappedData[referenceOrder,]

png(filename = output,width = 480,height = 1600)
layout(matrix(data=c(1,1,1,1,1,2),nrow=1))
image(1:ncol(mappedData),1:nrow(mappedData),t(mappedData),breaks=color_cuts,col=colorSpectrum,xaxt="n",yaxt="n",xlab="",ylab="")
image(1:2,color_cuts[2:101],t(matrix(data=color_cuts[2:101],ncol=2,nrow=100)),breaks=color_cuts,col=colorSpectrum,xaxt="n",xlab="",ylab="")
dev.off()
