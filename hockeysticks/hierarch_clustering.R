matrix_a <- read.table(file.choose(), sep='\t', header=T);

t <- cor(matrix_a, method="pearson")

t <- 1-t

matrix_d <- dist(t);

hc <- hclust(matrix_d,"average");

plot(hc)
hmcols<-colorRampPalette(c("red","white"))(256)

hclust.ave <- function(x) hclust(x, method="complete")

heatmap(t, Colv=T,Rowv=T, scale='none', col=hmcols, hclustfun=hclust.ave)

write.table(t, file="corr.csv")

svg("mymap.svg")
heatmap(...)
dev.off()