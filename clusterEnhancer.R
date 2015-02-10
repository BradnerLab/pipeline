#clusterEnhancer.R


#140121
#Charles Lin

# '''                                                                                                                                                    
# The MIT License (MIT)                                                                                                                                  
                                                                                                                                                       
# Copyright (c) 2014 Charles Lin                                                                                                                         
                                                                                                                                                       
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
# '''


# '''                                                                                                                                                    
# program to perform 2D clustering by enhancer signal                                                                                                    
# can perform initial enhancer mapping or draw from a set of finished rose outputs                                                                       
# '''

#===================================================================
#========================INPUT PARAMETERS===========================
#===================================================================

args <- commandArgs()
print(args[3:6])


genome = args[3]
outputFolder = args[4]
analysisName = args[5]
enhancerFile = args[6]




#===================================================================
#====================MAKING ENHANCER MATRIX=========================
#===================================================================



#enhancerTable = read.delim('/Volumes/raider/temp/enhancerClusterTesting/output/HG18_dlbcl_enhancers_signalTable.txt')
enhancerTable = read.delim(enhancerFile)


enhancerMatrix = as.matrix(enhancerTable[,7:ncol(enhancerTable)])
#rownames(enhancerMatrix) = enhancerTable[,1]
#===================================================================
#==================SAMPLE DISTANCE MATRICIES========================
#===================================================================

#distance by samples
sampleDist = as.dist(1-cor(enhancerMatrix))
sampleHC = hclust(sampleDist)


#===================================================================
#=====================PLOTTING SAMPLE TREE==========================
#===================================================================


#PLOTTING THE SAMPLE TREE

treePDFFile = paste(outputFolder,genome,'_',analysisName,"_treePlot.pdf",sep='')
pdf(treePDFFile,width=8,height=8)


plot(sampleHC,xlab='',main = paste(genome,'_',analysisName,sep=''))
dev.off()



#===================================================================
#=======================SAMPLE CLUSTER ORDERING=====================
#===================================================================

sampleOrder = sampleHC$order



#===================================================================
#=================MAKING SAMPLE PAIRWISE HEATMAP====================
#===================================================================



#distance by samples
sampleSimMatrix = cor(enhancerMatrix)


#Set the color spectrum
colorSpectrum <- colorRampPalette(c("white","red"))(100)

#setting a color data range
minValue <- .1
maxValue <- .9
color_cuts <- seq(minValue,maxValue,length=100)
color_cuts <- c(0, color_cuts,1)

#add one extra min color to even out sampling
colorSpectrum <- c(colorSpectrum[1],colorSpectrum)


clusterSampleFile = paste(outputFolder,genome,'_',analysisName,"_clusterSamples.pdf",sep='')

pdf(file = clusterSampleFile,width = 10,height =10)

layout(matrix(data=c(1,2,2,2,1,3,3,3,1,3,3,3,1,3,3,3,1,3,3,3,1,3,3,3,1,4,4,4),ncol= 7))
par(mar=c(8,9,5,9))
plot(as.dendrogram(sampleHC),ylab='Distance')
par(mar=c(4,2,2,2))

plot(as.dendrogram(sampleHC),horiz=TRUE,xlab='Distance',leaflab='none')
par(mar=c(6,4,4,2))

image(1:ncol(sampleSimMatrix),1:nrow(sampleSimMatrix),t(sampleSimMatrix[sampleOrder,sampleOrder]),breaks=color_cuts,col=colorSpectrum,xaxt="n",yaxt="n",ylab='',xlab='')

par(mar=c(6,5,4,2))

image(1:2,color_cuts[2:101],t(matrix(data=color_cuts[2:101],ncol=2,nrow=100)),breaks=color_cuts,col=colorSpectrum,xaxt="n",xlab="",ylab="Similarity")
dev.off()


#===================================================================
#===================ENHANCER DISTANCE MATRICIES=====================
#===================================================================

#distances by SEs
seDist = as.dist(1-cor(t(enhancerMatrix)))
seHC = hclust(seDist)
#plot(seHC)

#===================================================================
#======================ENHANCER CLUSTER ORDERING====================
#===================================================================

seOrder = seHC$order

#===================================================================
#======================INDIVIDUAL CLUSTERS==========================
#===================================================================

#establishing individual clusters
seClusterOrder = cutree(seHC,k=10)

# #try to guess the distance cutoff
# seDistRange = range(seDist)
# n = (seDistRange[2]-seDistRange[1])/20

# cutVector = c()
# nClustVector = c()
# for(i in seq(seDistRange[1]+n,seDistRange[2],n)){
	# n = max(cutree(seHC,h=i))
	
	# cutVector = c(cutVector,i)
	# nClustVector = c(nClustVector,n)
	
	
	
	
	
# }



#===================================================================
#======================WRITING CLUSTER TABLE========================
#===================================================================

clusterOrderTable = enhancerTable[which(seClusterOrder==1),]
for(i in 2:max(seClusterOrder)){

      clusterOrderTable =rbind(clusterOrderTable,enhancerTable[which(seClusterOrder==i),])

}

clusterColVector = c()
for(i in 1:max(seClusterOrder)){
	clusterColVector = c(clusterColVector,rep(i,length(which(seClusterOrder==i))))
	
	
}

clusterOrderTable = cbind(clusterOrderTable,clusterColVector)
colnames(clusterOrderTable)[ncol(clusterOrderTable)] = 'CLUSTER'

clusterOrderTable = cbind(clusterOrderTable,clusterColVector)
colnames(clusterOrderTable)[ncol(clusterOrderTable)] = 'CLUSTER'


clusterTableFile = paste(outputFolder,genome,'_',analysisName,"_clusterTable.txt",sep='')
write.table(clusterOrderTable,file= clusterTableFile,quote=FALSE,sep='\t',col.names=TRUE,row.names=FALSE)



#===================================================================
#======================PLOTTING EXEMPLARS===========================
#===================================================================
exemplarPDFFile = paste(outputFolder,genome,'_',analysisName,"_exemplarPlots.pdf",sep='')
pdf(exemplarPDFFile,width=5,height=5)

exemplarMatrix = matrix(nrow=max(seClusterOrder),ncol = ncol(enhancerMatrix))
#Plot the examplars
for(i in 1:max(seClusterOrder)){
	
	cluster = matrix(enhancerMatrix[which(seClusterOrder==i),],ncol=ncol(enhancerMatrix))
	plot(1:ncol(enhancerMatrix),rep(0,ncol(enhancerMatrix)),ylim =quantile(cluster,c(0.025,0.975)),cex=0,xlab='Samples',xaxt='n',ylab='Enhancer signal (fold/median)',main = paste('Cluster',i))
	axis(1,1:ncol(enhancerMatrix),colnames(enhancerMatrix)[sampleOrder])
	for(j in 1:nrow(cluster)){
		lines(1:ncol(enhancerMatrix),cluster[j,sampleOrder],col = rgb(0.2,0.2,0.2,0.1),lwd=1)
		
	}
	lines(1:ncol(enhancerMatrix),apply(cluster,2,median)[sampleOrder],col='red',lwd=4)
	legend(1,quantile(cluster,.95),paste("n =",nrow(cluster)))
	if(nrow(cluster) == 1){
		exemplarMatrix[i,] = cluster[1,sampleOrder]
	}else{
		exemplarMatrix[i,] = apply(cluster[,sampleOrder],2,median)
	}
}

colnames(exemplarMatrix) = colnames(enhancerMatrix)[sampleOrder]
rownames(exemplarMatrix) = 1:max(seClusterOrder)
dev.off()


#===================================================================
#======================WRITING EXEMPLAR TAB=========================
#===================================================================


exemplarTableFile = paste(outputFolder,genome,'_',analysisName,"_exmplarTable.txt",sep='')
write.table(signif(exemplarMatrix,4),file= exemplarTableFile,quote=FALSE,sep='\t')


#===================================================================
#=======================MAKING ENHANCER HEATMAPS====================
#===================================================================

#Set the color spectrum
colorSpectrum <- colorRampPalette(c("white","red"))(100)

#setting a color data range
minValue <- quantile(enhancerMatrix,na.rm=TRUE,prob=0.1,names=FALSE)
maxValue <- quantile(enhancerMatrix,na.rm=TRUE,prob=0.9,names=FALSE)
color_cuts <- seq(minValue,maxValue,length=100)
color_cuts <- c(0, color_cuts,max(enhancerMatrix))

#add one extra min color to even out sampling
colorSpectrum <- c(colorSpectrum[1],colorSpectrum)

#Making png

clusterPNGFile = paste(outputFolder,genome,'_',analysisName,'_2d_cluster.png',sep='')
png(filename = clusterPNGFile,width = 800,height =800)
layout(matrix(data=c(1,2,2,2,1,3,3,3,1,3,3,3,1,3,3,3,1,3,3,3,1,3,3,3,1,4,4,4),ncol= 7))
par(mar=c(7,9,2,9))
plot(as.dendrogram(sampleHC),ylab='Distance')
par(mar=c(4,2,2,0))

plot(as.dendrogram(seHC),horiz=TRUE,xlab='Distance',leaflab='none')
par(mar=c(6,2,4,2))

image(1:ncol(enhancerMatrix),1:nrow(enhancerMatrix),t(enhancerMatrix[seOrder,sampleOrder]),breaks=color_cuts,col=colorSpectrum,xaxt="n",yaxt="n",ylab='',xlab='')

par(mar=c(6,5,4,2))

image(1:2,color_cuts[2:101],t(matrix(data=color_cuts[2:101],ncol=2,nrow=100)),breaks=color_cuts,col=colorSpectrum,xaxt="n",xlab="",ylab="Fold vs. median")
dev.off()


#making pdf
clusterPDFFile = paste(outputFolder,genome,'_',analysisName,'_2d_cluster.pdf',sep='')

pdf(file = clusterPDFFile,width = 8,height =8)
layout(matrix(data=c(1,2,2,2,1,3,3,3,1,3,3,3,1,3,3,3,1,3,3,3,1,3,3,3,1,4,4,4),ncol= 7))
par(mar=c(7,9,2,9))
plot(as.dendrogram(sampleHC),ylab='Distance')
par(mar=c(4,2,2,0))

plot(as.dendrogram(seHC),horiz=TRUE,xlab='Distance',leaflab='none')
par(mar=c(6,2,4,2))

image(1:ncol(enhancerMatrix),1:nrow(enhancerMatrix),t(enhancerMatrix[seOrder,sampleOrder]),breaks=color_cuts,col=colorSpectrum,xaxt="n",yaxt="n",ylab='',xlab='')

par(mar=c(6,5,4,2))

image(1:2,color_cuts[2:101],t(matrix(data=color_cuts[2:101],ncol=2,nrow=100)),breaks=color_cuts,col=colorSpectrum,xaxt="n",xlab="",ylab="Fold vs. median")
dev.off()




#===================================================================
#====================MAKING EXEMPLAR HEATMAP========================
#===================================================================






#Set the color spectrum
colorSpectrum <- colorRampPalette(c("white","red"))(100)

#setting a color data range
minValue <- quantile(exemplarMatrix,na.rm=TRUE,prob=0.1,names=FALSE)
maxValue <- quantile(exemplarMatrix,na.rm=TRUE,prob=0.9,names=FALSE)
color_cuts <- seq(minValue,maxValue,length=100)
color_cuts <- c(0, color_cuts,max(exemplarMatrix))

#add one extra min color to even out sampling
colorSpectrum <- c(colorSpectrum[1],colorSpectrum)


clusterExemplarFile = paste(outputFolder,genome,'_',analysisName,"_clusterExemplars.pdf",sep='')
pdf(file = clusterExemplarFile,width = 8,height =8)
layout(matrix(data=c(1,2,2,2,1,2,2,2,1,2,2,2,1,2,2,2,1,3,3,3),ncol= 5))
par(mar=c(7,9,2,9))
plot(as.dendrogram(sampleHC),ylab='Distance')
par(mar=c(4,12,2,0))


image(1:ncol(exemplarMatrix),1:nrow(exemplarMatrix),t(exemplarMatrix),breaks=color_cuts,col=colorSpectrum,xaxt="n",yaxt="n",ylab='Clusters',xlab='')
axis(2,1:10,1:10)




par(mar=c(4,5,2,2))

image(1:2,color_cuts[2:101],t(matrix(data=color_cuts[2:101],ncol=2,nrow=100)),breaks=color_cuts,col=colorSpectrum,xaxt="n",xlab="",ylab="Fold vs. median")


dev.off()


