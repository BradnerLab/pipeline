#dynamicEnhancer_plot.R

#produces a ranked table from the region map table

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


#===========================================================
#======================DEBUG PARAM==========================
#===========================================================
#name1 = 'EC_CON_BRD4'
#name2 = 'EC_TNF_BRD4'
#stitchedFile = 'mergeTest/EC_BRD4_CON_ROSE/HG18_EC_MERGED_SUPERS_-0_+0_0KB_STITCHED_ENHANCER_REGION_MAP.txt'
#outFile = gsub('REGION_MAP','DELTA',stitchedFile)

#setwd('/Volumes/bradnerlab/projects/athero')


#stitched_regions=  read.delim('mergeTest/EC_BRD4_CON_ROSE/HG18_EC_MERGED_SUPERS_-0_+0_0KB_STITCHED_ENHANCER_REGION_MAP.txt',header=TRUE,sep='\t')

#===========================================================
#===============READING IN ARGUMNETS========================
#===========================================================

args <- commandArgs()

print(args[3:5])

stitchedFile = args[3]
name1 = args[4]
name2 = args[5]

#===========================================================
#=================MAKING DELTA TABLE========================
#===========================================================
outFile = gsub('REGION_MAP','DELTA',stitchedFile)
stitched_regions = read.delim(stitchedFile,header=TRUE)

factor1 = stitched_regions[,7] - stitched_regions[,8]
factor2 = stitched_regions[,9] - stitched_regions[,10]

factor1[factor1 <= 1] <- 1 
factor2[factor2 <= 1] <- 1 

deltaFactor = log2(factor2/factor1)
deltaFactor[is.infinite(deltaFactor)] <- -6
deltaFactor[deltaFactor < -8] <- -8
deltaFactor[deltaFactor > 8] <- 8

deltaOrder = order(deltaFactor,decreasing=TRUE)


#write the new table
newTable = cbind(stitched_regions[deltaOrder,1:6],deltaFactor[deltaOrder],factor1[deltaOrder],factor2[deltaOrder],1:length(deltaFactor),1)
colnames(newTable)[7:11] = c(paste('LOG2_',name2,'_VS_',name1,'_CHANGE',sep=''),name1,name2,'RANK','IS_SUPER')

write.table(newTable,outFile,quote=FALSE,sep='\t',row.names=FALSE)


#===========================================================
#===============PLOTTING WATERFALL==========================
#===========================================================


colorSpectrum <- colorRampPalette(c("green","black","black","red"))(100)

#setting a color data range
minValue <- -1.5
maxValue <- 1.5
color_cuts <- seq(minValue,maxValue,length=100)
color_cuts <- c(min(deltaFactor,na.rm=TRUE), color_cuts,max(deltaFactor,na.rm=TRUE))


#add one extra min color to even out sampling
colorSpectrum <- c(colorSpectrum[1],colorSpectrum[1],colorSpectrum)

colorVector = c()
for(i in deltaOrder){
	delta = deltaFactor[i]
	color = colorSpectrum[max(which(color_cuts <= delta))]
	colorVector =c(colorVector,color)
	
	
	
}

outPlot = gsub('txt','pdf',outFile)
pdf(file=outPlot,width =10,height=5)
plot(1:length(deltaFactor),deltaFactor[deltaOrder],ylim =c(-1.2*max(abs(deltaFactor)),1.2*max(abs(deltaFactor))),type='l',lwd=2,ylab=paste('Log2 change in signal ',name2,' over ',name1,sep=''),xlab=paste('Super-enhancer regions in either ',name1,' or ',name2,sep=''))
lines(1:length(deltaFactor),deltaFactor[deltaOrder],type='h',lwd=3,col=colorVector)
dev.off()
