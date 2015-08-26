library(affy)
library(graphics)



# The MIT License (MIT)

# Copyright (c) 2015 Charles Lin

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



#========================================================================
#===========================JOB PARAMETERS===============================
#========================================================================

args = commandArgs()
print(args)
geneFPKMFile = args[3]
outputFolder = args[4]
name = args[5]
groupString = args[6]

if(args[7] == 'TRUE'){
	useERCC = TRUE
	
}else{
	useERCC = FALSE
}

print(geneFPKMFile)
print(name)


groupVector=unlist(strsplit(groupString,','))

#========================================================================
#============================DEBUGGING===================================
#========================================================================

#========================================================================
#=========================HARD CODED STUFF===============================
#========================================================================
#set as path of the ERCC_Controls_Analysis.txt
erccTable = read.delim("/grail/genomes/ERCC_Technical_Data/ERCC_Controls_Analysis.txt")


#========================================================================
#=============================FUNCTIONS==================================
#========================================================================

#simple moving average
ma <- function(x,n=5){filter(x,rep(1/n,n), sides=2)}


## Add an alpha value to a colour
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
                     function(x) 
                       rgb(x[1], x[2], x[3], alpha=alpha))  
}

#panel function to do a scatter with a red diagonal line
panel.awesome <- function(x, y, col = par("col"), bg = NA, pch = par("pch"), 
    cex = 1, col.smooth = "red", span = 2/3, iter = 3, ...) 
{
    points(x, y, pch = pch, col = col, bg = bg, cex = cex,ylab='log2 expression (a.u.)',xlab='log2 expression (a.u.)')
    ok <- is.finite(x) & is.finite(y)
    if (any(ok)) 
        #lines(stats::lowess(x[ok], y[ok], f = span, iter = iter),col = 'red', ...)
        abline(a=0,b=1,lwd=2,col='red')
}

#panel function to do correlation
#adapted from http://www.r-bloggers.com/five-ways-to-visualize-your-pairwise-comparisons/
panel.cor <- function(x,y,digits=2,prefix="",...){
	usr <- par("usr"); on.exit(par(usr))
	par(usr=c(0,1,0,1))
	r <- abs(cor(x,y,method='spearman',use='complete'))
	txt <- round(r,4)
	txt <- paste(prefix,txt,sep="")
	cex <- 2
	test <- cor.test(x,y,method='spearman',use='complete')
	Signif <- symnum(test$p.value,corr=FALSE,na=FALSE,cutpoints = c(0,0.001,0.01,0.05,0.1,1),symbols = c("***","**","*",".","N.S"))
	text(0.5,0.5,txt,cex=cex*r)
	text(.8,.8,Signif,cex=cex,col=2)
	
}



#returns a vector of the concentrations for each ercc probe
plot_ercc <- function(erccTable,all_fpkm_exprs,tag){
	#first get the erccRows
	erccRows = grep("ERCC-",rownames(all_fpkm_exprs))
	erccList = rownames(all_fpkm_exprs)[erccRows]
	exprsRowVector = c()
	concVector = c()
	for(i in 1:length(erccList)){
		erccProbe = erccList[i]
		erccName = substr(erccProbe,1,10)
		#print(erccName)
		row = which(erccTable[,2] ==erccName)
		#print(row)
		if(length(row) >0){
			concentration = as.numeric(erccTable[row,4])
			#now check to see if probe is detected
			if(min(all_fpkm_exprs[erccRows[i],]) > 0){
				concVector = c(concVector,concentration)
				exprsRowVector = c(exprsRowVector,erccRows[i])
			}
		}
	}
	
	#now let's do some cute plotting
	plot(log10(concVector),log2(all_fpkm_exprs[exprsRowVector,1]),cex=0,xlab='log10 attomoles/ul',ylab='log2 expression (fpkm)',main=paste(tag,' spike-in expression',sep=""))
	palette = rainbow(ncol(all_fpkm_exprs),alpha=0.3)
	for(i in 1:ncol(all_fpkm_exprs)){
		#color = add.alpha(i,0.2)
		points(log10(concVector),log2(all_fpkm_exprs[exprsRowVector,i]),pch=19,col =add.alpha(i,0.2),cex=0.4)
		lines(loess.smooth(log10(concVector),log2(all_fpkm_exprs[exprsRowVector,i])),lwd=2,col=i)	
		
		
	}				
	legend(-1.5,.95*max(log2(all_fpkm_exprs[exprsRowVector,])),colnames(all_fpkm_exprs),col=1:ncol(all_fpkm_exprs),lwd=2)

}


magnitude <- function(x){
	
	
	return(sqrt((x[1])^2 + (log10(x[2]))^2))

}


makeMagnitudeMatrix <- function(m){
	
	magnitudeMatrix = cbind(m,apply(m[,3:4],1,magnitude))
	magnitudeOrder = order(magnitudeMatrix[,5],decreasing=TRUE)
	magnitudeMatrix = magnitudeMatrix[magnitudeOrder,]
	colnames(magnitudeMatrix)[5] = 'MAGNITUDE_POP_POP'
	return(magnitudeMatrix)
	
	
}


#makes a color vector to accompany a numeric vector x
makeColorVector <- function(x){
	
	
	colorSpectrum <- colorRampPalette(c("blue","grey","grey","red"))(100)

	#setting a color data range
	minValue <- -2
	maxValue <- 2
	color_cuts <- seq(minValue,maxValue,length=100)
	color_cuts <- c(min(x,na.rm=TRUE), color_cuts,max(x,na.rm=TRUE))
	
	
	#add one extra min color to even out sampling
	colorSpectrum <- c(colorSpectrum[1],colorSpectrum[1],colorSpectrum)
	
	colorVector = c()
	for(i in x){
		color = colorSpectrum[max(which(color_cuts <= i))]
		colorVector =c(colorVector,color)
	}
	return(colorVector)
}

#========================================================================
#============================DATA PROCESSING=============================
#========================================================================


#formatting the genes.fpkm file

all_fpkm_exprs = read.table(geneFPKMFile,header=TRUE)

print(all_fpkm_exprs[1:5,])

#gene row names must be unique.
#this finds all uniquely named rows (takes first instance)
usedGenes = c()
uniqueRows = c()
for(i in 1:nrow(all_fpkm_exprs)){
	geneName = as.character(all_fpkm_exprs[i,1])

	if(!(geneName %in% usedGenes)){
		uniqueRows = c(uniqueRows,i)
		usedGenes = c(usedGenes,geneName)
		
	}
	
}



if(length(uniqueRows) != nrow(all_fpkm_exprs)){
	print("WARNING: GENE ROW NAMES NOT UNIQUE. USING FIRST INSTANCE OF EACH GENE")
	}


#now get the unique gene row names

geneRowNames = as.character(all_fpkm_exprs[uniqueRows,1])

#now we need to remove any NAs
geneRowNames[which(is.na(geneRowNames))] <- 'GENE_NA'

#now subset the initial expression table
all_fpkm_exprs = all_fpkm_exprs[uniqueRows,2:ncol(all_fpkm_exprs)]
rownames(all_fpkm_exprs)= geneRowNames

#set a sane lower limit on expression
all_fpkm_exprs = apply(all_fpkm_exprs,c(1,2),function(x){x[intersect(which(x < 0.01),which(x >0))] = 0.01;x})

#write probe level expression raw
filename_raw = paste(outputFolder,name,'_all_fpkm_exprs_raw.txt',sep='')
write.table(all_fpkm_exprs,file=filename_raw,quote=FALSE,sep='\t')


if(useERCC == TRUE){
	subset=grep("ERCC-",rownames(all_fpkm_exprs))
	
	#epsilon adjustment to allow loess normalization to work
	#all_fpkm_exprs[subset,] = all_fpkm_exprs[subset,] +.1
	all_fpkm_exprs = all_fpkm_exprs+.1

	all_fpkm_exprs_norm <- loess.normalize(all_fpkm_exprs,subset=grep("ERCC-",rownames(all_fpkm_exprs)),log.it=TRUE,family.loess='gaussian')
	all_fpkm_exprs_norm[is.na(all_fpkm_exprs_norm)] <- 0
	#get rid of any negative values and set to 0
	all_fpkm_exprs_norm = apply(all_fpkm_exprs_norm,c(1,2),function(x){x[x<0] = 0;x})
	#set a postive expression floor of 0.01
	all_fpkm_exprs_norm = apply(all_fpkm_exprs_norm,c(1,2),function(x){x[intersect(which(x < 0.01),which(x >0))] = 0.01;x})
	


	#write probe level expression spikey normy
	filename_norm = paste(outputFolder,name,'_all_fpkm_exprs_norm.txt',sep='')
	write.table(all_fpkm_exprs_norm,file=filename_norm,quote=FALSE,sep='\t')
}


# 
#========================================================================
#====================BASIC ANALYSIS WITH ERCC============================
#========================================================================

if(useERCC == TRUE){
	

	#plotting spike-ins raw
	filename_spike = paste(outputFolder,name,'_spike_raw.pdf',sep='')
	pdf(file=filename_spike,width = 8,height =8)
	plot_ercc(erccTable,all_fpkm_exprs,'Raw')
	dev.off()
	
	
	#plotting spike-ins raw
	filename_spike = paste(outputFolder,name,'_spike_norm.pdf',sep='')
	pdf(file=filename_spike,width = 8,height =8)
	plot_ercc(erccTable,all_fpkm_exprs_norm,'Normalized')
	dev.off()
	
	#require an fpkm of at least 1 in at least 1 sample
	expressedProbesNorm = which(apply(all_fpkm_exprs_norm,1,max)>1)
	expressedProbesRaw = which(apply(all_fpkm_exprs,1,max)>1)

	png_size = 200 * ncol(all_fpkm_exprs_norm)

	#now do a pairwise scatter plot either raw or norm
	axisMinRaw = 0
	axisMaxRaw = max(log2(all_fpkm_exprs[expressedProbesRaw,]))
	filename_raw = paste(outputFolder,name,'_all_fpkm_exprs_raw_scatter.png',sep='')

	png(filename=filename_raw,width =png_size,height =png_size,pointsize=24)
	
	sampleRows = sample(expressedProbesRaw,1000)
	pairs(log2(all_fpkm_exprs[sampleRows,]),lower.panel=panel.awesome,upper.panel=panel.cor,cex.labels=0.8,xlim =c(axisMinRaw,axisMaxRaw),ylim = c(axisMinRaw,axisMaxRaw),pch=19,col=rgb(0.5,0.5,0.5,0.4),cex=1,main='Unnormalized log2 expression (fpkm)')
	dev.off()
	
	axisMinNorm = 0
	axisMaxNorm = max(log2(all_fpkm_exprs_norm[expressedProbesNorm,]))
	filename_norm = paste(outputFolder,name,'_all_fpkm_exprs_norm_scatter.png',sep='')
	png(filename=filename_norm,width =png_size,height =png_size,pointsize=24)
	
	sampleRows = sample(expressedProbesNorm,1000)
	pairs(log2(all_fpkm_exprs_norm[sampleRows,]),lower.panel=panel.awesome,upper.panel=panel.cor,cex.labels=0.8,xlim =c(axisMinNorm,axisMaxNorm),ylim = c(axisMinRaw,axisMaxRaw),pch=19,col=rgb(0.5,0.5,0.5,0.4),cex=1,main='Spike-in normalized log2 expression (fpkm)')
	dev.off()
	

	#now make some boxplots	
	filename_box = paste(outputFolder,name,'_exprs_boxplot.pdf',sep='')
	pdf(file=filename_box,width = 10,height = 8)
	par(mfrow=c(1,2))
	par(mar=c(12,6,3,1))
	axisMinBox = -4
	axisMaxBox = max(axisMaxRaw,axisMaxNorm)
	boxplot(log2(all_fpkm_exprs[expressedProbesRaw[1:1000],]),cex=0,main='Unnormalized expression',ylab='log2 expression (fpkm)',las=3,ylim = c(axisMinBox,axisMaxBox))
	
	boxplot(log2(all_fpkm_exprs_norm[expressedProbesNorm[1:1000],]),cex=0,main='Spike-in normalized expression',ylab='log2 expression (fpkm)',las=3,ylim = c(axisMinBox,axisMaxBox))
	dev.off()
}

#========================================================================
#=====================BASIC ANALYSIS WITHOUT ERCC========================
#========================================================================

if(useERCC == FALSE){

	#identify expressed probes
	#at least 1 probe above 1 fpkm
	expressedProbesRaw = which(apply(all_fpkm_exprs,1,max)>1)
	
	#provide a size scaling factor for the pngs
	png_size = 200 * ncol(all_fpkm_exprs)
	
	#now do a pairwise scatter plot either raw or norm
	axisMinRaw = 0
	axisMaxRaw = max(log2(all_fpkm_exprs[expressedProbesRaw,]))
	filename_raw = paste(outputFolder,name,'_all_fpkm_exprs_raw_scatter.png',sep='')
	
	png(filename=filename_raw,width =png_size,height =png_size,pointsize=24)
	pairs(log2(all_fpkm_exprs[expressedProbesRaw[1:1000],]),lower.panel=panel.awesome,upper.panel=panel.cor,cex.labels=0.8,xlim =c(axisMinRaw,axisMaxRaw),ylim = c(axisMinRaw,axisMaxRaw),pch=19,col=rgb(0.5,0.5,0.5,0.4),cex=1,main='Unnormalized log2 expression (fpkm)')
	dev.off()
	
	#now make some boxplots
	filename_box = paste(outputFolder,name,'_exprs_boxplot.pdf',sep='')
	pdf(file=filename_box,width = 10,height = 8)

	par(mar=c(12,6,3,1))
	axisMinBox = 0
	axisMaxBox = max(axisMaxRaw,axisMaxRaw)
	boxplot(log2(all_fpkm_exprs[expressedProbesRaw[1:1000],]),cex=0,main='Unnormalized expression',ylab='log2 expression (fpkm)',las=3,ylim = c(axisMinBox,axisMaxBox))	
	dev.off()
}



#========================================================================
#=================PICK A DATASET FOR DOWNSTREAM ANALYSIS=================
#========================================================================

if(useERCC){
	
	expressionTable = all_fpkm_exprs_norm
	
}else{
	
	expressionTable = all_fpkm_exprs
}

#focus only on expressed genes
#genes with an fpkm of at least 1 in 1 sample
#and genes with detectable expression in all samples
expressedGeneRows = intersect(which(apply(expressionTable,1,max)>1),which(apply(expressionTable,1,min)>0))

#setting axis limits for future plots to cover 99% of the data
axisMax=quantile(expressionTable[expressedGeneRows,1],0.995)

#set a floor of 0.5 fpkm for anything we want to use
axisMin=max(quantile(expressionTable[expressedGeneRows,1],0.005),0.5)

#get the sample names
sampleNames = colnames(expressionTable)

#========================================================================
#============================REPLICATE ANALYSIS==========================
#========================================================================

#make a multipage pdf with pairiwise scatters for replicates
#put the average correlation in the title

filename_replicates = paste(outputFolder,name,'_replicate_correlations.pdf',sep='')
pdf(file=filename_replicates,width =8,height =9)
for(group in groupVector){
	groupColumns = grep(group,sampleNames)
	groupColumns = sort(groupColumns)
	corVector = c()
	if(length(groupColumns) > 1){
		for(i in 1:length(groupColumns)){
			j=i+1
			while(j <= length(groupColumns)){			
				corVector = c(corVector,cor(expressionTable[expressedGeneRows,groupColumns[i]],expressionTable[expressedGeneRows,groupColumns[j]],use='complete',method='spearman'))
				j=j+1
			}	
		}
		
		avgCor=round(mean(corVector),4)
		figureTitle = paste('Replicates for ',group,' Avg. pairwise correlation of ',avgCor,sep='')
		
		
		par(mar=c(5.1,4.1,8.1,2.1))
	
		sampleRows = sample(expressedGeneRows,1000)#to keep scatter plots from bogging down
		pairs(log2(expressionTable[sampleRows,groupColumns]),lower.panel=panel.awesome,upper.panel=panel.cor,cex.labels=0.8,xlim =c(log2(axisMin),log2(axisMax)),ylim = c(log2(axisMin),log2(axisMax)),pch=19,col=rgb(0.5,0.5,0.5,0.4),cex=1,main=figureTitle)
	}
}
dev.off()


#========================================================================
#==========================MAKING MEAN MATRIX============================
#========================================================================

meanMatrixAll= matrix(nrow =nrow(expressionTable),ncol=length(groupVector))

for(i in 1:length(groupVector)){
	group = groupVector[i]
	groupColumns = grep(group,sampleNames)
	if(length(groupColumns) == 1){
		meanMatrixAll[,i] = expressionTable[,groupColumns]
	}else{
	meanMatrixAll[,i] = apply(expressionTable[,groupColumns],1,mean)
	}
}

rownames(meanMatrixAll) = rownames(expressionTable)
colnames(meanMatrixAll) = groupVector

filename_means = paste(outputFolder,name,'_all_fpkm_means.txt',sep='')
write.table(meanMatrixAll,file=filename_means,quote=FALSE,sep='\t')



meanMatrix= matrix(nrow =length(expressedGeneRows),ncol=length(groupVector))

for(i in 1:length(groupVector)){
	group = groupVector[i]
	groupColumns = grep(group,sampleNames)
	if(length(groupColumns) == 1){
		meanMatrix[,i] = expressionTable[expressedGeneRows,groupColumns]
	}else{
	meanMatrix[,i] = apply(expressionTable[expressedGeneRows,groupColumns],1,mean)
	}
}

rownames(meanMatrix) = rownames(expressionTable)[expressedGeneRows]
colnames(meanMatrix) = groupVector

filename_means = paste(outputFolder,name,'_exprs_fpkm_means.txt',sep='')
write.table(meanMatrix,file=filename_means,quote=FALSE,sep='\t')



#========================================================================
#=================CROSS COMPAIRSONS BETWEEN GROUPS=======================
#========================================================================

#for each group versus group
#make several kinds of figures
#1. volcano plot (w/ outliers highlighted)
#2 scatter w/ outliers annotated
#2. pairwise amplifier plots, genes ranked by expression in A with expression in B plotted and vice versa
#if a group has only 1 member, skip volcano plot and put NA in p-value column

#then spit out a GSEA style ranked table
#with all genes ranked by log2 change in expression w/ p-value annotated

#loop through pairwise comparisons
for(i in 1:length(groupVector)){
	j=i+1
	while(j <= length(groupVector)){			
		group1 = groupVector[i]
		group2 = groupVector[j]
		
		print(paste("Running analysis on ",group1," vs ",group2))
		print(paste("Running analysis on ",i," vs ",j))

		
		group1Columns = grep(group1,sampleNames)
		group2Columns = grep(group2,sampleNames)
		
		#set up the report pdf
		filename_pair = paste(outputFolder,name,'_',group1,'_vs_',group2,'.pdf',sep='')
		pdf(file = filename_pair,width = 11,height = 8.5)		
		
		expMatrix = matrix(nrow=length(expressedGeneRows),ncol=4)
		rownames(expMatrix) = rownames(expressionTable)[expressedGeneRows]
		colnames(expMatrix) = c(group1,group2,'LOG2_FOLD_CHANGE','P_VALUE')
		
		if(length(group1Columns)==1){
			expMatrix[,1] = expressionTable[expressedGeneRows,group1Columns]
		}else{
		expMatrix[,1] = apply(expressionTable[expressedGeneRows,group1Columns],1,mean)
		}
		
		if(length(group2Columns)==1){
			expMatrix[,2] = expressionTable[expressedGeneRows,group2Columns]
		}else{		
		expMatrix[,2] = apply(expressionTable[expressedGeneRows,group2Columns],1,mean)
		}
		expMatrix[,3] = log2(expMatrix[,2]/expMatrix[,1])
		
		#check to see if we have enough data to caluclate p-values
		#if so make the volcano plot
		if(min(length(group1Columns),length(group2Columns)) >1){
			pValueVector = c()
			for(n in expressedGeneRows){
			      	 
				expVector = c(expressionTable[n,group1Columns],expressionTable[n,group2Columns])
				if(min(expVector) == max(expVector)){
					pValue = 1
				}else{
					pValue = t.test(expressionTable[n,group1Columns],expressionTable[n,group2Columns])$p.value
				}
				pValueVector = c(pValueVector,pValue)	
			}
			expMatrix[,4] = pValueVector		
			
			#make volcano plot only w/ pvalue vector
			#find the min
			minPvalue = log10(expMatrix[which.min(log10(expMatrix[,4])),4])
			if(is.infinite(minPvalue)){
				minPvalue = -10
			}
			xTitle = paste('Log2 fold change ',group2,' vs. ',group1,sep='')
			xMin = quantile(expMatrix[,3],0.005) -2
			xMax = quantile(expMatrix[,3],0.995) +2

			plot(expMatrix[,3],log10(expMatrix[,4]),xlim = c(xMin,xMax),ylim =c(0,minPvalue),cex=0.5,pch=19,col=rgb(0.5,0.5,0.5,0.2),xlab=xTitle,ylab='Log10 p-value')
			abline(h=log10(0.05),lty=2)
			abline(v=1,lty=2)
			abline(v=-1,lty=2)
			downRows = intersect(which(expMatrix[,3] < -1),which(log10(expMatrix[,4]) < log10(0.05)))
			upRows = intersect(which(expMatrix[,3] > 1),which(log10(expMatrix[,4]) < log10(0.05)))
			points(expMatrix[downRows,3],log10(expMatrix[downRows,4]),pch=19,col=rgb(0,0,1,.2),cex=0.8)
			points(expMatrix[upRows,3],log10(expMatrix[upRows,4]),pch=19,col=rgb(1,0,0,.2),cex=0.8)

			
			#plot top 10 genes w/ max divergence
			#circle
			downMagnitudeMatrix = makeMagnitudeMatrix(expMatrix[downRows,])
			upMagnitudeMatrix = makeMagnitudeMatrix(expMatrix[upRows,])
			
			if(nrow(downMagnitudeMatrix) > 10){
				points(downMagnitudeMatrix[1:10,3],log10(downMagnitudeMatrix[1:10,4]),pch=19,col=rgb(0,0,1),cex=1)			
				for(n in 1:10){
				text(downMagnitudeMatrix[n,3],log10(downMagnitudeMatrix[n,4]),rownames(downMagnitudeMatrix)[n],pos=2,col='blue')
				}
			}

			if(nrow(upMagnitudeMatrix) > 10){
				points(upMagnitudeMatrix[1:10,3],log10(upMagnitudeMatrix[1:10,4]),pch=19,col=rgb(1,0,0),cex=1)
				for(n in 1:10){
				text(upMagnitudeMatrix[n,3],log10(upMagnitudeMatrix[n,4]),rownames(upMagnitudeMatrix)[n],pos=4,col='red')
				}								
			}
			#label w/ text top5 of each
		

			
		}
		
		axisLimits = c(log2(axisMin),log2(axisMax))
		#next is the scatter
		plot(log2(expMatrix[,1]),log2(expMatrix[,2]),xlim=axisLimits,ylim=axisLimits,xlab=paste(group1,'log2 fpkm'),ylab=paste(group2,'log2 fpkm'),pch=19,cex=0.5,col=rgb(0.5,0.5,0.5,0.2))
		abline(a=0,b=1,lty=2)
		abline(a=1,b=1,col='red')
		abline(a=-1,b=1,col='blue')
		#if we have pvalues
		if(min(length(group1Columns),length(group2Columns)) >1){
			upSig = intersect(which(expMatrix[,3]>1),which(expMatrix[,4]<0.05))
			points(log2(expMatrix[upSig,1]),log2(expMatrix[upSig,2]),pch=19,cex=0.8,col=rgb(1,0,0,.2))
			downSig = intersect(which(expMatrix[,3]< -1),which(expMatrix[,4]<0.05))
			points(log2(expMatrix[downSig,1]),log2(expMatrix[downSig,2]),pch=19,cex=0.8,col=rgb(0,0,1,.2))		
		}
		
		#next is the waterfall
		changeOrder = order(expMatrix[,3])
		yTitle = paste('Log2 fold change ',group2,' vs. ',group1,sep='')
		xTitle = paste('Genes ranked by increasing Log2 fold change ',group2,' vs. ',group1,sep='')
		yLimits=c(min(quantile(expMatrix[,3],c(0.0005)),-1),max(quantile(expMatrix[,3],c(0.9995)),1))
		
		plot(1:length(changeOrder),expMatrix[changeOrder,3],ylim= yLimits,type='l',ylab=yTitle,xlab=xTitle)
		abline(h=0)
		colorVector = makeColorVector(expMatrix[changeOrder,3])
		lines(1:length(changeOrder),expMatrix[changeOrder,3],type='h',lwd=3,col=colorVector)
		
		upCount =length(which(expMatrix[,3] > 1))
		downCount =length(which(expMatrix[,3]< -1))
		abline(v = length(changeOrder)-upCount,col='red')
		abline(v = downCount,col='blue')
		text(length(changeOrder)-upCount,-0.5,col='red',paste(upCount,'\nup'),pos=4)
		text(downCount,0.5,col='blue',paste(downCount,'\ndown'),pos=2)

		#last is amplifier plot
		par(mfrow=c(1,2))
		#do linear and log?
		group1RankOrder = order(expMatrix[,1])
		group2RankOrder = order(expMatrix[,2])
		binSize = length(group1RankOrder)/50
		logString=''
		
		#ranked by group1
		plot(1:length(group1RankOrder),log2(expMatrix[group1RankOrder,1]),ylim =c(log2(axisMin),log2(max(expMatrix,na.rm=TRUE))),xlab=paste('Genes ranked by expression in',group1),ylab='Log2 expression (fpkm)',lwd=2,log=logString,type='l',col = 'grey',xlim = c(1,length(group1RankOrder)))
		points(1:length(group1RankOrder), log2(expMatrix[group1RankOrder,2]),col=add.alpha('red',.1),pch=16,cex=0.5)
		legend(0,.7*log2(max(expMatrix,na.rm=TRUE)),c(group1,group2),lwd=2,col=c('grey','red'))
		x= ma(1:length(group1RankOrder),n=binSize)	
		x = x[is.finite(x)]
		y= ma(log2(expMatrix[group1RankOrder,2]),n=binSize)
		y = y[is.finite(y)]
		lines(x,y,col='red',lwd=2)
		
		#ranked by group2
		plot(1:length(group2RankOrder),log2(expMatrix[group2RankOrder,2]),ylim =c(log2(axisMin),log2(max(expMatrix,na.rm=TRUE))),xlab=paste('Genes ranked by expression in',group2),ylab='Log2 expression (fpkm)',lwd=2,log=logString,type='l',col = 'grey')
		points(1:length(group2RankOrder),log2(expMatrix[group2RankOrder,1]),col=add.alpha('red',.1),pch=16,cex=0.5)
		legend(0,.7*log2(max(expMatrix,na.rm=TRUE)),c(group2,group1),lwd=2,col=c('grey','red'))
		x= ma(1:length(group2RankOrder),n=binSize)	
		x = x[is.finite(x)]
		y= ma(log2(expMatrix[group2RankOrder,1]),n=binSize)
		y = y[is.finite(y)]
		lines(x,y,col='red',lwd=2)
	
		#close the plot
		dev.off()
		
		#now write out the exp table
		
		filename_exp= paste(outputFolder,name,'_',group1,'_vs_',group2,'_exprs_matrix.txt',sep='')
		write.table(expMatrix,file=filename_exp,quote=FALSE,sep='\t')
		
		
		#making the gct
		filename_gct= paste(outputFolder,name,'_',group1,'_vs_',group2,'.gct',sep='')
		gctMatrix =matrix(ncol=4,nrow=nrow(expMatrix))
		colnames(gctMatrix) = c('NAME','DESCRIPTION',group1,group2)
		gctMatrix[,1]= rownames(expMatrix)
		gctMatrix[,3]= expMatrix[,1]
		gctMatrix[,4]=expMatrix[,2]
		
		gctHeader = matrix(data='',ncol=4,nrow=3)
		gctHeader[1,1]='#1.2'
		gctHeader[2,1]=nrow(expMatrix)
		gctHeader[2,2]='2'
		gctHeader[3,]=c('NAME','DESCRIPTION',group1,group2)
		gctCombined = rbind(gctHeader,gctMatrix)
		write.table(gctCombined,file=filename_gct,quote=FALSE,sep='\t',row.names=FALSE,col.names=FALSE)		

		#making the cls
		filename_cls= paste(outputFolder,name,'_',group1,'_vs_',group2,'.cls',sep='')
		clsTable = matrix(data='',ncol=3,nrow=3)
		clsTable[1,] =c(2,2,1)
		clsTable[2,1]=paste('#',group1,sep='')
		clsTable[2,2]=group2
		clsTable[3,1]=group1
		clsTable[3,2]=group2
		write.table(clsTable,file=filename_cls,quote=FALSE,sep='\t',row.names=FALSE,col.names=FALSE)		

		
		
				
		j=j+1
		}	
	}
