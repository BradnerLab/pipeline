stitchFile = commandArgs(TRUE)[1]

outFolder = commandArgs(TRUE)[2]

name = commandArgs(TRUE)[3]

#print(stitchFile)
#print(outFolder)
#print(name)


stitchTable = read.table(stitchFile,header=TRUE)





xTitle = 'Number of enhancer regions'
xVector =stitchTable$NUM_REGIONS

yVector=stitchTable$MEAN_CONSTIT/stitchTable$MEAN_REGION
yTitle = 'Average fraction of \nenhancer sequence in each region'
yDerivTitle = 'Decrease in fraction of \nenhancer sequence per region'
yDeriv = diff(yVector)

# print('x vector has a length of')
# print(length(xVector))
# print(xVector)
# print(xVector[2:length(xVector)])
# print(length(xVector[2:length(xVector)]))
# print('y spline vector has a length of')
# print(length(yDeriv))
# print(yVector)
# print(yDeriv)

yDerviFit = smooth.spline(xVector[2:length(xVector)],yDeriv)
optX= yDerviFit$x[which.min(yDerviFit$y)]

#get the optimal stitching parameter
optStitch = stitchTable[which(stitchTable[,2]==optX),1]


stitchPDFName = paste(outFolder,name,'_stitch_parameter.pdf',sep='')
#print(stitchPDFName)

pdf(file=stitchPDFName,width=8.5,height=11)
par(mai=c(1,1.5,.2,.5))
par(mfrow=c(2,1))
plot(xVector,yVector,type='l',xlab=xTitle,ylab=yTitle)
abline(v= optX)
stitchText = paste('OPTIMUM STITCHING AT: ',optStitch,'bp',sep='')
text(min(xVector),.8*max(yVector),stitchText,pos=4)

plot(xVector[2:length(xVector)],yDeriv,xlab=xTitle,ylab=yDerivTitle)
lines(yDerviFit,col='blue')
abline(v= optX)

dev.off()
write(optStitch,stdout())