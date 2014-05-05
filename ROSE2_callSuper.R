#============================================================================
#==============SUPER-ENHANCER CALLING AND PLOTTING FUNCTIONS=================
#============================================================================

 #This function calculates the cutoff by sliding a diagonal line and finding where it is tangential (or as close as possible)
 calculate_cutoff <- function(inputVector, drawPlot=TRUE,...){
 	inputVector <- sort(inputVector)
	inputVector[inputVector<0]<-0 #set those regions with more control than ranking equal to zero
	slope <- (max(inputVector)-min(inputVector))/length(inputVector) #This is the slope of the line we want to slide. This is the diagonal.
	xPt <- floor(optimize(numPts_below_line,lower=1,upper=length(inputVector),myVector= inputVector,slope=slope)$minimum) #Find the x-axis point where a line passing through that point has the minimum number of points below it. (ie. tangent)
	y_cutoff <- inputVector[xPt] #The y-value at this x point. This is our cutoff.
	
	if(drawPlot){  #if TRUE, draw the plot
		plot(1:length(inputVector), inputVector,type="l",...)
		b <- y_cutoff-(slope* xPt)
		abline(v= xPt,h= y_cutoff,lty=2,col=8)
		points(xPt,y_cutoff,pch=16,cex=0.9,col=2)
		abline(coef=c(b,slope),col=2)
		title(paste("x=",xPt,"\ny=",signif(y_cutoff,3),"\nFold over Median=",signif(y_cutoff/median(inputVector),3),"x\nFold over Mean=",signif(y_cutoff/mean(inputVector),3),"x",sep=""))
		axis(1,sum(inputVector==0),sum(inputVector==0),col.axis="pink",col="pink") #Number of regions with zero signal
	}
	return(list(absolute=y_cutoff,overMedian=y_cutoff/median(inputVector),overMean=y_cutoff/mean(inputVector)))
}

#this is an accessory function, that determines the number of points below a diagnoal passing through [x,yPt]
numPts_below_line <- function(myVector,slope,x){
	yPt <- myVector[x]
	b <- yPt-(slope*x)
	xPts <- 1:length(myVector)
	return(sum(myVector<=(xPts*slope+b)))
}


convert_stitched_to_bed <- function(inputStitched,trackName,trackDescription,outputFile,splitSuper=TRUE,score=c(),superRows=c(),baseColor="0,0,0",superColor="255,0,0"){
	outMatrix <- matrix(data="",ncol=4+ifelse(length(score)==nrow(inputStitched),1,0),nrow=nrow(inputStitched))
	
	outMatrix[,1] <- as.character(inputStitched$CHROM)
	outMatrix[,2] <- as.character(inputStitched$START)
	outMatrix[,3] <- as.character(inputStitched$STOP)
	outMatrix[,4] <- as.character(inputStitched$REGION_ID)
	if(length(score)==nrow(inputStitched)){
		score <- rank(score,ties.method="first")
		score <- length(score)-score+1  #Stupid rank only does smallest to largest. 
		outMatrix[,5] <- as.character(score)
	}
	trackDescription <- paste(trackDescription,"\nCreated on ",format(Sys.time(), "%b %d %Y"),collapse="",sep="")
	trackDescription <- gsub("\n","\t", trackDescription)
	tName <- gsub(" ","_",trackName)
	cat('track name="', tName,'" description="', trackDescription,'" itemRGB=On color=',baseColor,"\n",sep="",file=outputFile)
	write.table(file= outputFile,outMatrix,sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE)
	if(splitSuper==TRUE){
		cat("\ntrack name=\"Super_", tName,'" description="Super ', trackDescription,'" itemRGB=On color=', superColor,"\n",sep="",file=outputFile,append=TRUE)
		write.table(file= outputFile,outMatrix[superRows,],sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE)
	}
}







writeSuperEnhancer_table <- function(superEnhancer,description,outputFile,additionalData=NA){
	description <- paste("#",description,"\nCreated on ",format(Sys.time(), "%b %d %Y"),collapse="",sep="")
	description <- gsub("\n","\n#",description)
	cat(description,"\n",file=outputFile)
	if(is.matrix(additionalData)){
		if(nrow(additionalData)!=nrow(superEnhancer)){
			warning("Additional data does not have the same number of rows as the number of super enhancers.\n--->>> ADDITIONAL DATA NOT INCLUDED <<<---\n")
		}else{
			superEnhancer <- cbind(superEnhancer,additionalData)
			superEnhancer = superEnhancer[order(superEnhancer$enhancerRank),]
			
		}
	}
	write.table(file=outputFile,superEnhancer,sep="\t",quote=FALSE,row.names=FALSE,append=TRUE)
}



#============================================================================
#============================HELPER FUNCTIONS================================
#============================================================================

#http://stackoverflow.com/questions/9837766/r-plot-circle-with-radius-1-and-angle-0-2pi-in-polar-coordinates
circle <- function(x, y, rad = 1, nvert = 500, ...){
    rads <- seq(0,2*pi,length.out = nvert)
    xcoords <- cos(rads) * rad + x
    ycoords <- sin(rads) * rad + y
    polygon(xcoords, ycoords, ...)
}


magnitude <- function(x,y){
	magnitudeVector=c()
	for(i in 1:length(x)){
		
		magnitudeVector = c(magnitudeVector,sqrt((x[i])^2 + (y[i])^2))

	}
	return(magnitudeVector)

}



geneToRefseq <- function(geneName,transcribedTable){
	refseqIDs = c()
	rowID = which(transcribedTable[,3] == geneName)
	for(row in rowID){
		refseqIDs = c(refseqIDs,as.character(transcribedTable[row,2]))
	}
	
	return(refseqIDs)

	
}


#get the row by enhancer ID for enhancer tables that are sorted uniquely
enhancerIDToRow <- function(enhancerID,targetTable){
	return(which(targetTable[,1]==enhancerID))
	
}


#gets genes associated w/ an enhancer by ID
getEnhancerGenes <- function(enhancerID,enhancerTable){
	
	enhancerGenes = c()
	row = enhancerIDToRow(enhancerID,enhancerTable)
		
	foo = as.character(enhancerTable[row,7])
	enhancerGenes = c(enhancerGenes,unlist(strsplit(foo,',')))
		foo = as.character(enhancerTable[row,8])
		enhancerGenes = c(enhancerGenes,unlist(strsplit(foo,',')))
		foo = as.character(enhancerTable[row,9])
		enhancerGenes = c(enhancerGenes,unlist(strsplit(foo,',')))
	
	enhancerGenes = unique(enhancerGenes)
	return(enhancerGenes)
	
}


getRefseqIDs <- function(enhancerIDList,enhancerTable,transcribedTable){
	
	refIDs = c()
	for(enhancerID in enhancerIDList){
		
		enhancerGenes = getEnhancerGenes(enhancerID,enhancerTable)
		for(geneName in enhancerGenes){
			
			refIDs = c(refIDs,geneToRefseq(geneName,transcribedTable))
			
		}
	}
	#print(refIDs)
	return(refIDs)
}



#============================================================================
#===================SUPER-ENHANCER CALLING AND PLOTTING======================
#============================================================================


#============================================================================
#==============================INPUT ARGUMENTS===============================
#============================================================================



args <- commandArgs()

print('THESE ARE THE ARGUMENTS')
print(args)

#ARGS
outFolder = args[3]
enhancerFile = args[4]
enhancerName = args[5]
wceName = args[6]


#============================================================================
#================================WORKSHOPPING================================
#============================================================================

# setwd('/Volumes/grail/projects/newRose/')



# enhancerFile = 'roseTest/MM1S_H3K27AC_DMSO_peaks_6KB_STITCHED_TSS_DISTAL_ENHANCER_REGION_MAP.txt'

# outFolder = 'roseTest/'

# #wceName = 'MM1S_WCE_DMSO_MERGED.hg18.bwt.sorted.bam'
# wceName = 'NONE'
# enhancerName = 'MM1S_H3K27AC_DMSO_peaks'

#============================================================================
#================================DATA INPUT==================================
#============================================================================


#Read enhancer regions with closestGene columns
stitched_regions <- read.delim(file= enhancerFile,sep="\t")


#perform WCE subtraction. Using pipeline table to match samples to proper background. 
rankBy_factor = colnames(stitched_regions)[7]

prefix = unlist(strsplit(rankBy_factor,'_'))[1]

if(wceName == 'NONE'){
	
	
	rankBy_vector = as.numeric(stitched_regions[,7])
	
}else{
	wceName = colnames(stitched_regions)[8]
	print('HERE IS THE WCE NAME')
	print(wceName)
	
	rankBy_vector = as.numeric(stitched_regions[,7])-as.numeric(stitched_regions[,8])
}	
	

#SETTING NEGATIVE VALUES IN THE rankBy_vector to 0

rankBy_vector[rankBy_vector < 0] <- 0

#============================================================================
#======================SETTING ORIGINAL ROSE CUTOFFS=========================
#============================================================================



#FIGURING OUT THE CUTOFF

cutoff_options <- calculate_cutoff(rankBy_vector, drawPlot=FALSE,xlab=paste(rankBy_factor,'_enhancers'),ylab=paste(rankBy_factor,' Signal','- ',wceName),lwd=2,col=4)


#These are the super-enhancers
superEnhancerRows <- which(rankBy_vector> cutoff_options$absolute)
typicalEnhancers = setdiff(1:nrow(stitched_regions),superEnhancerRows)
enhancerDescription <- paste(enhancerName," Enhancers\nCreated from ", enhancerFile,"\nRanked by ",rankBy_factor,"\nUsing cutoff of ",cutoff_options$absolute," for Super-Enhancers",sep="",collapse="")


#============================================================================
#========================MAKING SUPER HOCKEY STICK===========================
#============================================================================


#MAKING HOCKEY STICK PLOT
plotFileName = paste(outFolder,enhancerName,'_Plot_points.png',sep='')
png(filename=plotFileName,height=600,width=600)
signalOrder = order(rankBy_vector,decreasing=TRUE)
if(wceName == 'NONE'){
	plot(length(rankBy_vector):1,rankBy_vector[signalOrder], col='red',xlab=paste(rankBy_factor,'_enhancers'),ylab=paste(rankBy_factor,' Signal'),pch=19,cex=2)	
	
}else{
	plot(length(rankBy_vector):1,rankBy_vector[signalOrder], col='red',xlab=paste(rankBy_factor,'_enhancers'),ylab=paste(rankBy_factor,' Signal','- ',wceName),pch=19,cex=2)
}
abline(h=cutoff_options$absolute,col='grey',lty=2)
abline(v=length(rankBy_vector)-length(superEnhancerRows),col='grey',lty=2)
lines(length(rankBy_vector):1,rankBy_vector[signalOrder],lwd=4, col='red')
text(0,0.8*max(rankBy_vector),paste(' Cutoff used: ',cutoff_options$absolute,'\n','Super-Enhancers identified: ',length(superEnhancerRows)),pos=4)

dev.off()



#============================================================================
#======================SETTING STRETCH ROSE CUTOFFS==========================
#============================================================================



#FIGURING OUT THE CUTOFF

stretch_vector = abs(as.numeric(stitched_regions[,6]))

stretch_cutoff_options <- calculate_cutoff(stretch_vector, drawPlot=FALSE,xlab=paste(rankBy_factor,'_enhancers'),ylab=paste(rankBy_factor,' enhancer lengths'),lwd=2,col=4)


#These are the stretch-enhancers
stretchEnhancerRows <- which(stretch_vector > stretch_cutoff_options$absolute)
typicalStretchEnhancers = setdiff(1:nrow(stitched_regions), stretchEnhancerRows)
stretchEnhancerDescription <- paste(enhancerName," Enhancers\nCreated from ", enhancerFile,"\nRanked by ",rankBy_factor," lengths\nUsing cutoff of ", stretch_cutoff_options$absolute," for Stretch-Enhancers",sep="",collapse="")


#============================================================================
#=========================MAKING STRETCH ROSE PLOTS==========================
#============================================================================

#MAKING HOCKEY STICK PLOT
plotFileName = paste(outFolder,enhancerName,'_Plot_points_stretch.png',sep='')
png(filename=plotFileName,height=600,width=600)
signalOrder = order(stretch_vector,decreasing=TRUE)

plot(length(stretch_vector):1, stretch_vector[signalOrder], col='red',xlab=paste(rankBy_factor,'_enhancers'),ylab=paste(rankBy_factor,' lengths (bp)'),pch=19,cex=2)	
	

abline(h=stretch_cutoff_options$absolute,col='grey',lty=2)
abline(v=length(stretch_vector)-length(stretchEnhancerRows),col='grey',lty=2)
lines(length(stretch_vector):1, stretch_vector[signalOrder],lwd=4, col='red')
text(0,0.8*max(stretch_vector),paste(' Cutoff used: ',stretch_cutoff_options$absolute,'\n','Stretch-Enhancers identified: ',length(stretchEnhancerRows)),pos=4)

dev.off()



#============================================================================
#================================MAKING PANEL PLOTS==========================
#============================================================================

#MAKING NEW HOCKEY STICK PLOT
plotFileName = paste(outFolder,enhancerName,'_Plot_panel.png',sep='')
png(filename=plotFileName,height=600,width=1200)
par(mfrow= c(1,3))


#FIRST THE HOCKEY
signalOrder = order(rankBy_vector,decreasing=TRUE)
enhancerOrder = signalOrder
plot(length(rankBy_vector):1, rankBy_vector[enhancerOrder], col='red',xlab='Enhancers ranked by increasing signal',ylab='Enhancer signal (total rpm)',lwd=2,type='l')

points(length(rankBy_vector):(length(rankBy_vector)-length(superEnhancerRows)+1),rankBy_vector[enhancerOrder[1:length(superEnhancerRows)]],pch=19,cex=1,col='red')
points((length(rankBy_vector)-length(superEnhancerRows)):1,rankBy_vector[enhancerOrder[(length(superEnhancerRows)+1):length(enhancerOrder)]],pch=19,cex=0.75,col='grey')

abline(h=cutoff_options$absolute,col=rgb(0.3,0.3,0.3),lty=2)
abline(v=length(rankBy_vector)-length(superEnhancerRows),col=rgb(0.3,0.3,0.3),lty=2)
text(0,0.8*max(rankBy_vector),paste(' Cutoff used: ',cutoff_options$absolute,'\n','Super-Enhancers identified: ',length(superEnhancerRows)),pos=4)



#THEN THE SCATTER

allSEs = union(superEnhancerRows,stretchEnhancerRows)
superStretch = intersect(superEnhancerRows,stretchEnhancerRows)


enhMagnitude = magnitude(stretch_vector[allSEs]/max(stretch_vector),rankBy_vector[allSEs]/max(rankBy_vector))

m =  as.matrix(cbind(stretch_vector[allSEs]/max(stretch_vector),rankBy_vector[allSEs]/max(rankBy_vector)))
mDiag = apply(m,1,sum)/2
mDist = sqrt(2*(m[,1]-mDiag)^2)
mDist[which(m[,2] > m[,1])] <- mDist[which(m[,2] > m[,1])]*-1



plot(mDist, enhMagnitude,cex=0.75,col='grey',ylim =c(-.05,1),xlim = c(-0.5,0.5),xlab='Enhancer skew',ylab='Enhancer combined magnitude')

ssSubset = c()
for(x in 1:length(allSEs)){
	if(length(which(superStretch == allSEs[x])) > 0){
		ssSubset = c(ssSubset,x)
		
	}
	
}

points(mDist[ssSubset],enhMagnitude[ssSubset],pch=19,cex=1,col='red')
abline(h=0)
abline(v=0)

text(0,-.05,"MORE SUPER",pos=2)
text(0,-.05,"MORE STRETCH",pos=4)
legend(-.5,.95,c(paste(length(superStretch),'SUPER AND STRETCH')),pch=19,col='red')



#THEN STRETCH
signalOrder = order(stretch_vector,decreasing=FALSE)
enhancerOrder = signalOrder


plot(1:length(stretch_vector), stretch_vector[rev(enhancerOrder)], col='red',xlab='Enhancers ranked by decreasing length',ylab='Enhancer length (bm)',lwd=2,type='l')

points(1:length(stretchEnhancerRows), stretch_vector[enhancerOrder[length(stretch_vector):(length(stretch_vector)-length(stretchEnhancerRows)+1)]],pch=19,cex=1,col='red')
points(length(stretchEnhancerRows):length(stretch_vector), stretch_vector[enhancerOrder[(length(typicalStretchEnhancers)+1):1]],pch=19,cex=0.75,col='grey')

abline(h=stretch_cutoff_options$absolute,col=rgb(0.3,0.3,0.3),lty=2)
abline(v=length(stretchEnhancerRows),col=rgb(0.3,0.3,0.3),lty=2)
text(length(stretch_vector),0.8*max(stretch_vector),paste(' Cutoff used: ',stretch_cutoff_options$absolute,'\n','Stretch-Enhancers identified: ',length(stretchEnhancerRows)),pos=2)

dev.off()


#============================================================================
#============================WRITING SUPER OUTPUT============================
#============================================================================


#Writing a bed file
bedFileName = paste(outFolder,enhancerName,'_Enhancers_withSuper.bed',sep='')
convert_stitched_to_bed(stitched_regions,paste(rankBy_factor,"Enhancers"), enhancerDescription,bedFileName,score=rankBy_vector,splitSuper=TRUE,superRows= superEnhancerRows,baseColor="0,0,0",superColor="255,0,0")

#This matrix is just the super_enhancers
true_super_enhancers <- stitched_regions[superEnhancerRows,]

additionalTableData <- matrix(data=NA,ncol=2,nrow=nrow(stitched_regions))
colnames(additionalTableData) <- c("enhancerRank","isSuper")
additionalTableData[,1] <- nrow(stitched_regions)-rank(rankBy_vector,ties.method="first")+1
additionalTableData[,2] <- 0
additionalTableData[superEnhancerRows,2] <- 1


#Writing enhancer and super-enhancer tables with enhancers ranked and super status annotated
enhancerTableFile = paste(outFolder,enhancerName,'_AllEnhancers.table.txt',sep='')
writeSuperEnhancer_table(stitched_regions, enhancerDescription,enhancerTableFile, additionalData= additionalTableData)

superTableFile = paste(outFolder,enhancerName,'_SuperEnhancers.table.txt',sep='')
writeSuperEnhancer_table(true_super_enhancers, enhancerDescription,superTableFile, additionalData= additionalTableData[superEnhancerRows,])



#============================================================================
#============================WRITING STRETCH ROSE============================
#============================================================================

#Writing a bed file
bedFileName = paste(outFolder,enhancerName,'_Enhancers_withStretch.bed',sep='')
convert_stitched_to_bed(stitched_regions,paste(rankBy_factor,"Enhancers"), enhancerDescription,bedFileName,score= stretch_vector,splitSuper=TRUE,superRows= stretchEnhancerRows,baseColor="0,0,0",superColor="255,0,0")



#This matrix is just the super_enhancers
true_stretch_enhancers <- stitched_regions[stretchEnhancerRows,]

additionalTableData <- matrix(data=NA,ncol=2,nrow=nrow(stitched_regions))
colnames(additionalTableData) <- c("enhancerRank","isStretch")
additionalTableData[,1] <- nrow(stitched_regions)-rank(stretch_vector,ties.method="first")+1
additionalTableData[,2] <- 0
additionalTableData[stretchEnhancerRows,2] <- 1


#Writing enhancer and stretch-enhancer tables with enhancers ranked and stretch status annotated
enhancerTableFile = paste(outFolder,enhancerName,'_AllEnhancers_Length.table.txt',sep='')
writeSuperEnhancer_table(stitched_regions, enhancerDescription,enhancerTableFile, additionalData= additionalTableData)

stretchTableFile = paste(outFolder,enhancerName,'_StretchEnhancers.table.txt',sep='')
writeSuperEnhancer_table(true_stretch_enhancers, enhancerDescription,stretchTableFile, additionalData= additionalTableData[stretchEnhancerRows,])



#============================================================================
#================================WRITING 2D ROSE=============================
#============================================================================

#Writing a bed file
bedFileName = paste(outFolder,enhancerName,'_Enhancers_withSuperStretch.bed',sep='')
convert_stitched_to_bed(stitched_regions,paste(rankBy_factor,"Enhancers"), enhancerDescription,bedFileName,score= stretch_vector,splitSuper=TRUE,superRows= superStretch,baseColor="0,0,0",superColor="255,0,0")





#This matrix is just the super_enhancers
true_superStretch_enhancers <- stitched_regions[superStretch,]
print(length(superStretch))
print(dim(true_superStretch_enhancers))
additionalTableData <- matrix(data=NA,ncol=2,nrow=nrow(stitched_regions))
colnames(additionalTableData) <- c("enhancerRank","isSuperStretch")


enhMagnitude = magnitude(stretch_vector/max(stretch_vector),rankBy_vector/max(rankBy_vector))
additionalTableData[,1] <- nrow(stitched_regions)-rank(enhMagnitude,ties.method="first")+1
additionalTableData[,2] <- 0
additionalTableData[superStretch,2] <- 1


#Writing enhancer and superStretch-enhancer tables with enhancers ranked and superStretch status annotated
enhancerTableFile = paste(outFolder,enhancerName,'_AllEnhancers_SuperStretch.table.txt',sep='')
writeSuperEnhancer_table(stitched_regions, enhancerDescription,enhancerTableFile, additionalData= additionalTableData)

superStretchTableFile = paste(outFolder,enhancerName,'_SuperStretchEnhancers.table.txt',sep='')
writeSuperEnhancer_table(true_superStretch_enhancers, enhancerDescription,superStretchTableFile, additionalData= additionalTableData[superStretch,])

