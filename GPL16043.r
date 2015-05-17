library(affy)
library(makecdfenv)
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



#========================================================================
#===========================JOB PARAMETERS===============================
#========================================================================

args = commandArgs()
cel_file_directory = args[3]
name = args[4]

print(cel_file_directory)
print(name)
#quit()
#name = "JB20130926st"
#cel_file_directory <- "/ark/home/cl512/ressrv19/raw/expression/JB20130926st/"

#========================================================================
#=========================HARD CODED STUFF===============================
#========================================================================
erccTable = read.delim("/grail/genomes/ERCC_Technical_Data/ERCC_Controls_Analysis.txt")

#primeviewcdf_env <- make.cdf.env(filename="PrimeView_withERCC_binary.cdf", cdf.path="/ark/home/cl512/ressrv19/annotations/platforms/GPL16043/annotation/",compress=FALSE)

primeviewcdf_env <- make.cdf.env(filename="PrimeView_withERCC_binary.cdf", cdf.path="/grail/annotations/platforms/GPL16043/annotation/",compress=FALSE)

#========================================================================
#=============================FUNCTIONS==================================
#========================================================================
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
	r <- abs(cor(x,y,method='spearman'))
	txt <- round(r,2)
	txt <- paste(prefix,txt,sep="")
	cex <- 2
	test <- cor.test(x,y,method='spearman')
	Signif <- symnum(test$p.value,corr=FALSE,na=FALSE,cutpoints = c(0,0.001,0.01,0.05,0.1,1),symbols = c("***","**","*",".","N.S"))
	text(0.5,0.5,txt,cex=cex*r)
	text(.8,.8,Signif,cex=cex,col=2)
	
}



#returns a vector of the concentrations for each ercc probe
plot_ercc <- function(erccTable,all_mas5_exprs,tag){
	#first get the erccRows
	erccRows = grep("ERCC-",rownames(all_mas5_exprs))
	erccList = rownames(all_mas5_exprs)[erccRows]
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
			concVector = c(concVector,concentration)
			exprsRowVector = c(exprsRowVector,erccRows[i])
		}
	}
	
	#now let's do some cute plotting
	plot(log10(concVector),log2(all_mas5_exprs[exprsRowVector]),cex=0,xlab='log10 attomoles/ul',ylab='log2 expression (a.u.)',main=paste(tag,' spike-in expression',sep=""))
	palette = rainbow(ncol(all_mas5_exprs),alpha=0.3)
	for(i in 1:ncol(all_mas5_exprs)){
		#color = add.alpha(i,0.2)
		points(log10(concVector),log2(all_mas5_exprs[exprsRowVector,i]),pch=19,col =add.alpha(i,0.2),cex=0.4)
		lines(loess.smooth(log10(concVector),log2(all_mas5_exprs[exprsRowVector,i])),lwd=2,col=i)	
		
		
	}				
	legend(-1.5,.95*max(log2(all_mas5_exprs[exprsRowVector])),colnames(all_mas5_exprs),col=1:ncol(all_mas5_exprs),lwd=2)

}






#========================================================================
#============================DATA PROCESSING=============================
#========================================================================


cel_files <- list.celfiles(path= cel_file_directory,full.names=TRUE)
raw_data <- read.affybatch(cel_files, cdfname="primeviewcdf_env")

mas5Result <- expresso(raw_data,bgcorrect.method="mas",normalize=TRUE,pmcorrect.method="pmonly",summary.method="mas")


all_mas5_exprs <- exprs(mas5Result)
all_mas5_exprs_norm <- loess.normalize(all_mas5_exprs,subset=grep("ERCC-",rownames(all_mas5_exprs)))

#write probe level expression raw
filename_raw = paste(cel_file_directory,'output/',name,'_all_mas5_probe_exprs_raw.txt',sep='')
write.table(all_mas5_exprs,file=filename_raw,quote=FALSE,sep='\t')

#write probe level expression spikey normy
filename_raw = paste(cel_file_directory,'output/',name,'_all_mas5_probe_exprs_norm.txt',sep='')
write.table(all_mas5_exprs_norm,file=filename_raw,quote=FALSE,sep='\t')

#========================================================================
#================================ANALYSIS================================
#========================================================================

#plotting spike-ins raw
filename_spike = paste(cel_file_directory,'output/',name,'_spike_raw.pdf',sep='')
pdf(file=filename_spike,width = 8,height =8)
plot_ercc(erccTable,all_mas5_exprs,'Raw')
dev.off()


#plotting spike-ins raw
filename_spike = paste(cel_file_directory,'output/',name,'_spike_norm.pdf',sep='')
pdf(file=filename_spike,width = 8,height =8)
plot_ercc(erccTable,all_mas5_exprs_norm,'Normalized')
dev.off()


#identify expressed probes
#at least 1 probe above 50
expressedProbesRaw = which(apply(all_mas5_exprs,1,max)>100)
expressedProbesNorm = which(apply(all_mas5_exprs_norm,1,max)>100)

#provide a size scaling factor for the pngs
png_size = 200 * ncol(all_mas5_exprs)

#now do a pairwise scatter plot either raw or norm
axisMinRaw = min(log2(all_mas5_exprs[expressedProbesRaw,]))
axisMaxRaw = max(log2(all_mas5_exprs[expressedProbesRaw,]))
filename_raw = paste(cel_file_directory,'output/',name,'_all_mas5_probe_exprs_raw_scatter.png',sep='')

png(filename=filename_raw,width =png_size,height =png_size,pointsize=24)
pairs(log2(all_mas5_exprs[expressedProbesRaw[1:1000],]),lower.panel=panel.awesome,upper.panel=panel.cor,cex.labels=0.8,xlim =c(axisMinRaw,axisMaxRaw),ylim = c(axisMinRaw,axisMaxRaw),pch=19,col=rgb(0.5,0.5,0.5,0.4),cex=1,main='Unnormalized log2 expression (a.u.)')
dev.off()

#now do a pairwise scatter plot either raw or norm
axisMinNorm = min(log2(all_mas5_exprs_norm[expressedProbesNorm,]))
axisMaxNorm = max(log2(all_mas5_exprs_norm[expressedProbesNorm,]))
filename_norm = paste(cel_file_directory,'output/',name,'_all_mas5_probe_exprs_norm_scatter.png',sep='')
png(filename=filename_norm,width =png_size,height =png_size,pointsize=24)
pairs(log2(all_mas5_exprs_norm[expressedProbesNorm[1:1000],]),lower.panel=panel.awesome,upper.panel=panel.cor,cex.labels=0.8,xlim =c(axisMinNorm,axisMaxNorm),ylim = c(axisMinRaw,axisMaxRaw),pch=19,col=rgb(0.5,0.5,0.5,0.4),cex=1,main='Spike-in normalized log2 expression (a.u.)')
dev.off()



#now make some boxplots

filename_box = paste(cel_file_directory,'output/',name,'_probe_exprs_boxplot.pdf',sep='')
pdf(file=filename_box,width = 10,height = 8)
par(mfrow=c(1,2))
par(mar=c(12,6,3,1))
axisMinBox = min(axisMinRaw,axisMinNorm)
axisMaxBox = max(axisMaxRaw,axisMaxNorm)
boxplot(log2(all_mas5_exprs[expressedProbesRaw[1:1000],]),cex=0,main='Unnormalized expression',ylab='log2 expression (a.u.)',las=3,ylim = c(axisMinBox,axisMaxBox))

boxplot(log2(all_mas5_exprs_norm[expressedProbesNorm[1:1000],]),cex=0,main='Spike-in normalized expression',ylab='log2 expression (a.u.)',las=3,ylim = c(axisMinBox,axisMaxBox))
dev.off()

