#setwd('/Volumes/young_ata4/myc_111311/bamPlot/')
library(graphics)

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



args <- commandArgs()

print(args[3:7])


nameTable = read.delim(args[3],header=FALSE)
diagramTable = read.delim(args[4],header=FALSE)
plotTable = read.delim(args[5])
yScale = args[6]
plotStyle = args[7]
fileName = args[8]

plotHeight = (nrow(plotTable)+1)*3

nameFile = unlist(strsplit(args[3],'/'))
outFolder = paste(as.vector(nameFile[2:(length(nameFile)-1)]),collapse='/')
inputName = nameFile[length(nameFile)]
start = unlist(strsplit(inputName,'_'))[3]
sense = unlist(strsplit(inputName,'_'))[2]
end = unlist(strsplit(inputName,'_'))[4]
chrom = unlist(strsplit(inputName,'_'))[1]

print('fileName is')
print(fileName)
print(length(fileName))
if(nchar(fileName) > 1){
  outputName = fileName}else{outputName = paste(chrom,sense,start,end,'plot.pdf',sep='_')}
outFile = paste('',outFolder,outputName,sep='/')

print(outFile)






nBins = length(plotTable[1,])-7
yMinDiagram = min(diagramTable[,2]-3)

#first bring in the colors
colorVector = c()
for(i in 1:nrow(plotTable)){
	color = rgb(plotTable[i,5],plotTable[i,6],plotTable[i,7],maxColorValue=255)
	colorVector = c(colorVector,color)
}


#do the plotting
pdf(file=outFile,width = 8.5,height =plotHeight)

if(plotStyle == 'SINGLE'){
	m = matrix(c(2,2,2,2,2,2,2,2,1,1,1),nrow=11,ncol=8)	
	layout(m)
	#plotting the diagram
	plot(0,0,xlim = c(0,nBins),ylim = c(yMinDiagram,2),col=rgb(1,1,1),xaxt='n',yaxt='n',ylab='',xlab='',main ='')
	for(i in 2:nrow(diagramTable)){
		rect(diagramTable[i,1],diagramTable[i,2],diagramTable[i,3],diagramTable[i,4],col='black')
	
	
	}
	
	#plotting the names
	for(i in 2:nrow(nameTable)){
		text(nameTable[i,2],nameTable[i,3],nameTable[i,1],cex=1)
	}

	#for all on the same plot
	

	yMax = 	1.2*max(plotTable[1,(8:(nBins+7))])
		
	if(yScale =='RELATIVE'){
		color = colorVector[1]
		plot(spline(1:nBins,as.numeric(plotTable[1,(8:(nBins+7))]),n=nBins),ylim = c(0.1,yMax),type='l',col= color,lwd=2.5,xaxt='n',yaxt='n',xlab='',ylab='Relative peak heights')
		if(sense =='-'){
			axis(1,at = c(0,nBins),labels= c(paste(chrom,end,sep=':'),paste(chrom,start,sep=':')))
			}else{
			axis(1,at = c(0,nBins),labels= c(paste(chrom,start,sep=':'),paste(chrom,end,sep=':')))
			}
		legend(0,yMax,as.vector(plotTable[,3]),col=colorVector,lwd=2.5,cex=1.2)

		for(i in 2:nrow(plotTable)){
			scaleFactor = max(plotTable[1,(8:(nBins+7))])/(1.2*max(plotTable[i,(8:(nBins+7))]))

			print(scaleFactor)
			#scaleFactor = 1
			color = colorVector[i]
			lines(spline(1:nBins,scaleFactor*as.numeric(plotTable[i,(8:(nBins+7))]),n=3*nBins),lwd=2,col = color)
			}
		
	}else{
		color = colorVector[1]  ## BJA tweaked to style of RELATIVE
		plot(spline(1:nBins,as.numeric(plotTable[1,(8:(nBins+7))]),n=nBins),ylim = c(0.1,yMax),type='l',col= color,lwd=2.5,xaxt='n',xlab='',ylab='ChIP-Seq reads')
		if(sense =='-'){
			axis(1,at = c(0,nBins),labels= c(paste(chrom,end,sep=':'),paste(chrom,start,sep=':')))
			}else{
				axis(1,at = c(0,nBins),labels= c(paste(chrom,start,sep=':'),paste(chrom,end,sep=':')))
			}
		legend(0,yMax,as.vector(plotTable[,3]),col=colorVector,lwd=2.5,cex=1.2)

		for(i in 2:nrow(plotTable)){
			color = colorVector[i] ## BJA tweaked to style of RELATIVE
# 			color = rgb(plotTable[i,5],plotTable[i,6],plotTable[i,7],maxColorValue=255)
			lines(spline(1:nBins,as.numeric(plotTable[i,(8:(nBins+7))]),n=3*nBins),lwd=2,col = color)
			}
		}
		
	}





#for different plots
if(plotStyle == 'MULTIPLE'){
	par(mfrow = c(nrow(plotTable)+1,1))
	if(yScale == 'UNIFORM'){
		yMax = 1.2*max(plotTable[,(8:(nBins+7))])
	}
	for(i in 1:nrow(plotTable)){
		if(yScale == 'RELATIVE'){
			yMax = 1.2*max(plotTable[i,(8:(nBins+7))])
		}
		color = colorVector[i]
		plot(spline(1:nBins,as.numeric(plotTable[i,(8:(nBins+7))]),n=2*nBins),ylim = c(0.05*yMax,yMax),type='l',col= color,lwd=2,xlab='',ylab='ChIP-Seq Reads',xaxt = 'n')
		legend(0,yMax,as.vector(plotTable[i,3]),col=colorVector[i],lwd=2.5,cex=1.2)
		if(sense =='-'){
			axis(1,at = c(0,nBins),labels= c(paste(chrom,end,sep=':'),paste(chrom,start,sep=':')))
			}else{
			axis(1,at = c(0,nBins),labels= c(paste(chrom,start,sep=':'),paste(chrom,end,sep=':')))
			}
		}

	plot(0,0,xlim = c(0,nBins),ylim = c(yMinDiagram,2),col=rgb(1,1,1),xaxt='n',yaxt='n',ylab='',xlab='',main ='')
	for(i in 2:nrow(diagramTable)){
		rect(diagramTable[i,1],diagramTable[i,2],diagramTable[i,3],diagramTable[i,4],col='black')
		}
	for(i in 2:nrow(nameTable)){
		text(nameTable[i,2],nameTable[i,3],nameTable[i,1],cex=1)
	}
	
	
}

dev.off()
