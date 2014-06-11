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


print(args[3:6])

summaryFile = args[3]
outFile = args[4]
yScale = args[5]
plotStyle = args[6]

#==========================================================
#==================DEBUG SECTION===========================
#==========================================================
#setwd('/Users/charles/Dropbox/src/temp/')
#summaryFile = 'ACTB/ACTB_summary.txt'
#outFile = '/Users/charles/Dropbox/src/temp/ACTB_plots_test.pdf'
#yScale = 'UNIFORM'
#plotStyle = 'MULTIPLE'

#==========================================================
#==========================================================
#==========================================================


#obtain the outfile name
print('fileName is')
print(outFile)

#load in the summary file
summaryTable = read.delim(summaryFile)


#establish the plot height from the first entry
plotFile = as.character(summaryTable$PLOT_TABLE[1])
plotTable = read.delim(plotFile)
plotHeight = (nrow(plotTable)+1)*3

#now open up the pdf
pdf(file=outFile,width = 8.5,height =plotHeight)

#now loop through the summary table
for(i in 1:nrow(summaryTable)){
	
	diagramFile = as.character(summaryTable$DIAGRAM_TABLE[i])
	nameFile = as.character(summaryTable$NAME_TABLE[i])
	bedDiagramFile = as.character(summaryTable$BED_DIAGRAM_TABLE[i])
	bedNameFile = as.character(summaryTable$BED_NAME_TABLE[i])
	plotFile = as.character(summaryTable$PLOT_TABLE[i])
	
	diagramTable = read.delim(diagramFile,header=FALSE)
	nameTable = read.delim(nameFile,header=FALSE)
	bedDiagramTable = read.delim(bedDiagramFile,header=FALSE)
	bedNameTable = read.delim(bedNameFile,header=FALSE)
	plotTable = read.delim(plotFile)
	
	chrom = as.character(summaryTable$CHROM[i])
	name = as.character(summaryTable$ID[i])
	sense = as.character(summaryTable$SENSE[i])
	start = as.numeric(summaryTable$START[i])
	end = as.numeric(summaryTable$END[i])
	
	#check if the beds have any data
	if(nrow(bedNameTable) >1){
		hasBed=TRUE
		print("Plotting with a bed file")
	}else{hasBed=FALSE}

	#don't attempt to plot regions w/o data
	if(is.na(plotTable[1,8])){
		next
		}

	nBins = length(plotTable[1,])-7
	yMinDiagram = min(diagramTable[,2]-3)
	yMinBedDiagram = min(bedDiagramTable[,2])

	#first bring in the colors
	colorVector = c()
	for(i in 1:nrow(plotTable)){
		color = rgb(plotTable[i,5],plotTable[i,6],plotTable[i,7],maxColorValue=255)
		colorVector = c(colorVector,color)
	}
	
	#now the actual plotting
	if(plotStyle == 'SINGLE'){
		if(hasBed){
			m = matrix(c(3,3,3,3,3,2,2,1,1),nrow=9,ncol=8)	
		}else{
			m = matrix(c(2,2,2,2,2,2,2,2,1,1,1),nrow=11,ncol=8)	
		}
		layout(m)
		#plotting the diagram
		par(mai=c(0,1.5,0.2,0.2772))
		plot(0,0,xlim = c(0,nBins),ylim = c(yMinDiagram,2),col=rgb(1,1,1),xaxt='n',yaxt='n',ylab='',xlab='',main ='')
		for(i in 2:nrow(diagramTable)){
			rect(diagramTable[i,1],diagramTable[i,2],diagramTable[i,3],diagramTable[i,4],col='black')
		}
		#plotting the names
		for(i in 2:nrow(nameTable)){
			text(nameTable[i,2],nameTable[i,3],nameTable[i,1],cex=1)
		}
		if(hasBed){
			#plotting the beds
			par(mai=c(0,1.5,.2,0.2772))
			plot(0,0,xlim = c(0,nBins),ylim = c(yMinBedDiagram,.5),col=rgb(1,1,1),xaxt='n',yaxt='n',ylab='',xlab='',main ='')
			for(i in 2:nrow(bedDiagramTable)){
				rect(bedDiagramTable[i,1],bedDiagramTable[i,2],bedDiagramTable[i,3],bedDiagramTable[i,4],col='black')
			}
			
			#the bed names		
			axis(2,bedNameTable[2:nrow(bedNameTable),3],labels=bedNameTable[2:nrow(bedNameTable),1],las=1)
		}
		
		#for all on the same plot
		yMax = 	1.2*max(plotTable[1,(8:(nBins+7))])
		par(mai=c(0.1,1.5,0.1,0.2772))	
		if(yScale =='RELATIVE'){
			#establish a blank plot			
			plotSpline = spline(1:nBins,scaleFactor*as.numeric(plotTable[1,(8:(nBins+7))]),n=2*nBins)			
			xVector = c(1,plotSpline$x,max(plotSpline$x))

			plot(0,0,ylim = c(0.05*yMax,yMax),cex=0,xlim = range(xVector),xlab='',ylab='Relative peak heights',xaxt = 'n',main=name)

			if(sense =='-'){
				axis(1,at = c(0,nBins),labels= c(paste(chrom,end,sep=':'),paste(chrom,start,sep=':')))
				}else{
				axis(1,at = c(0,nBins),labels= c(paste(chrom,start,sep=':'),paste(chrom,end,sep=':')))
				}
			legend(0,yMax,as.vector(plotTable[,3]),col=colorVector,lwd=2.5,cex=1.2)
	
			for(i in 1:nrow(plotTable)){
				scaleFactor = max(plotTable[1,(8:(nBins+7))])/(1.2*max(plotTable[i,(8:(nBins+7))]))
				color = colorVector[i]
				plotSpline = spline(1:nBins,scaleFactor*as.numeric(plotTable[i,(8:(nBins+7))]),n=2*nBins)			
				xVector = c(1,plotSpline$x,max(plotSpline$x))
				yVector = c(0,plotSpline$y,0)			
				polygon(xVector,yVector,col= color,lty=0)
				}
			
		}else{
			color = colorVector[1]  ## BJA tweaked to style of RELATIVE

			#establish a blank plot
			plotSpline = spline(1:nBins,as.numeric(plotTable[1,(8:(nBins+7))]),n=2*nBins)
			xVector = c(1,plotSpline$x,max(plotSpline$x))
			yVector = c(0,plotSpline$y,0)			

			plot(0,0,ylim = c(0.05*yMax,yMax),cex=0,xlim = range(xVector),xlab='',ylab='Relative peak heights',xaxt = 'n',main=name)

			if(sense =='-'){
				axis(1,at = c(0,nBins),labels= c(paste(chrom,end,sep=':'),paste(chrom,start,sep=':')))
				}else{
					axis(1,at = c(0,nBins),labels= c(paste(chrom,start,sep=':'),paste(chrom,end,sep=':')))
				}
			legend(0,yMax,as.vector(plotTable[,3]),col=colorVector,lwd=2.5,cex=1.2)
	
			for(i in 1:nrow(plotTable)){
				color = colorVector[i] ## BJA tweaked to style of RELATIVE
	# 			color = rgb(plotTable[i,5],plotTable[i,6],plotTable[i,7],maxColorValue=255)
				plotSpline = spline(1:nBins,as.numeric(plotTable[i,(8:(nBins+7))]),n=2*nBins)
				xVector = c(1,plotSpline$x,max(plotSpline$x))
				yVector = c(0,plotSpline$y,0)			
				polygon(xVector,yVector,col= color,lty=0)

				}
			}
			
		}
	
	
	#for different plots
	if(plotStyle == 'MULTIPLE'){
		if(hasBed){
			par(mfrow = c(nrow(plotTable)+2,1))
		}else{
			par(mfrow = c(nrow(plotTable)+1,1))

		}
		par(mai=c(0.2,1.5,0.2,0.2772))	

		if(yScale == 'UNIFORM'){
			yMax = 1.2*max(plotTable[,(8:(nBins+7))])
		}
		for(i in 1:nrow(plotTable)){
			if(yScale == 'RELATIVE'){
				yMax = 1.2*max(plotTable[i,(8:(nBins+7))])
			}
			color = colorVector[i]
			plotSpline = spline(1:nBins,as.numeric(plotTable[i,(8:(nBins+7))]),n=2*nBins)
	
			xVector = c(1,plotSpline$x,max(plotSpline$x))
			yVector = c(0,plotSpline$y,0)			
			plot(0,0,ylim = c(0.05*yMax,yMax),cex=0,xlim = range(xVector),xlab='',ylab='ChIP-Seq Reads',xaxt = 'n',main=name)
			polygon(xVector,yVector,col= color,lty=0)
			legend(0,yMax,as.vector(plotTable[i,3]),col=colorVector[i],lwd=2.5,cex=1.2)
			if(sense =='-'){
				axis(1,at = c(0,nBins),labels= c(paste(chrom,end,sep=':'),paste(chrom,start,sep=':')))
				}else{
				axis(1,at = c(0,nBins),labels= c(paste(chrom,start,sep=':'),paste(chrom,end,sep=':')))
				}
			}
		#plotting the beds
		if(hasBed){
			par(mai=c(0,1.5,0.2,0.2772))
			plot(0,0,xlim = c(0,nBins),ylim = c(yMinBedDiagram,.5),col=rgb(1,1,1),xaxt='n',yaxt='n',ylab='',xlab='',main ='')
			for(i in 2:nrow(bedDiagramTable)){
				rect(bedDiagramTable[i,1],bedDiagramTable[i,2],bedDiagramTable[i,3],bedDiagramTable[i,4],col='black')
			}
			
			#the bed names		
			axis(2,bedNameTable[2:nrow(bedNameTable),3],labels=bedNameTable[2:nrow(bedNameTable),1],las=1)
		}
		#the gene diagram
		par(mai=c(0,1.5,0.2,0.2772))

		plot(0,0,xlim = c(0,nBins),ylim = c(yMinDiagram,2),col=rgb(1,1,1),xaxt='n',yaxt='n',ylab='',xlab='',main ='')
		for(i in 2:nrow(diagramTable)){
			rect(diagramTable[i,1],diagramTable[i,2],diagramTable[i,3],diagramTable[i,4],col='black')
			}
		for(i in 2:nrow(nameTable)){
			text(nameTable[i,2],nameTable[i,3],nameTable[i,1],cex=1)
		}


		
	}

	
}


dev.off()
