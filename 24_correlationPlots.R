# File: 24_correlationPlots.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: create correlation matrix plots for selected datasets
# Date: 21/5/2021

source('header.R')

lFiles = list.files('results/', pattern='DEAnalysis*', full.names = T, ignore.case = T)

ldfData = lapply(lFiles, function(x) as.data.frame(read.csv(x, header=T, row.names=1, stringsAsFactors = F)))
names(ldfData) = lFiles
sapply(ldfData, nrow)

# put everything in one order by row names
rn = rownames(ldfData[[1]])
head(rn)

ldfData = lapply(ldfData, function(df){
  df = df[rn,]
})

sapply(ldfData, function(df) identical(rownames(df), rn))

cvTitle = gsub('results//DEAnalysis_', '', names(ldfData))
cvTitle = gsub('.xls', '', cvTitle)

names(ldfData)
## select significant genes
dfContrast1.sub = ldfData[[1]][ldfData[[1]]$adj.P.Val < 0.01 & abs(ldfData[[1]]$logFC) > 0.5,]
dfContrast2.sub = ldfData[[2]][ldfData[[2]]$adj.P.Val < 0.01 & abs(ldfData[[2]]$logFC) > 0.5,]
dfContrast3.sub = ldfData[[3]][ldfData[[3]]$adj.P.Val < 0.01 & abs(ldfData[[3]]$logFC) > 0.5,]
dfContrast4.sub = ldfData[[4]][ldfData[[4]]$adj.P.Val < 0.01 & abs(ldfData[[4]]$logFC) > 0.5,]
dfContrast5.sub = ldfData[[5]][ldfData[[5]]$adj.P.Val < 0.01 & abs(ldfData[[5]]$logFC) > 0.5,]
dfContrast6.sub = ldfData[[6]][ldfData[[6]]$adj.P.Val < 0.01 & abs(ldfData[[6]]$logFC) > 0.5,]
dfContrast7.sub = ldfData[[7]][ldfData[[7]]$adj.P.Val < 0.01 & abs(ldfData[[7]]$logFC) > 0.5,]
dfContrast8.sub = ldfData[[8]][ldfData[[8]]$adj.P.Val < 0.01 & abs(ldfData[[8]]$logFC) > 0.5,]
dfContrast9.sub = ldfData[[9]][ldfData[[9]]$adj.P.Val < 0.01 & abs(ldfData[[9]]$logFC) > 0.5,]
dfContrast10.sub = ldfData[[10]][ldfData[[10]]$adj.P.Val < 0.01 & abs(ldfData[[10]]$logFC) > 0.5,]

# create a list for overlaps
lVenn = list(rownames(dfContrast1.sub), rownames(dfContrast2.sub), rownames(dfContrast3.sub),
             rownames(dfContrast4.sub), rownames(dfContrast5.sub), rownames(dfContrast6.sub),
             rownames(dfContrast7.sub), rownames(dfContrast8.sub), rownames(dfContrast9.sub),
             rownames(dfContrast10.sub)
)
names(ldfData)
names(lVenn) = cvTitle
cvCommonGenes = Reduce(intersect, lVenn)

## create a binary matrix
cvCommonGenes = unique(do.call(c, lVenn))
mCommonGenes = matrix(NA, nrow=length(cvCommonGenes), ncol=length(lVenn))
colnames(mCommonGenes) = cvTitle

for (i in 1:ncol(mCommonGenes)){
  mCommonGenes[,i] = ldfData[[i]][cvCommonGenes, 'logFC'] 
}
rownames(mCommonGenes) = gsub('X', '', cvCommonGenes)
dim(mCommonGenes)

## create correlation plot for all DEGs
library(NMF)
library(RColorBrewer)

m = (cor(mCommonGenes))
aheatmap(m, annRow = NA, scale = 'none', Rowv = NA, breaks=0.4,
         Colv=NA, cexRow=1, cexCol = 1,  
         col=(brewer.pal(9, 'YlGnBu')))

############
## repeat the figure with WT induced
source('header.R')

lFiles = list.files('results/', pattern='DEAnalysis*', full.names = T, ignore.case = T)
lFiles = lFiles[c(1, 3, 5, 7)]

ldfData = lapply(lFiles, function(x) as.data.frame(read.csv(x, header=T, row.names=1, stringsAsFactors = F)))
names(ldfData) = lFiles
sapply(ldfData, nrow)

# put everything in one order by row names
rn = rownames(ldfData[[1]])
head(rn)

ldfData = lapply(ldfData, function(df){
  df = df[rn,]
})

sapply(ldfData, function(df) identical(rownames(df), rn))

cvTitle = gsub('results//DEAnalysis_', '', names(ldfData))
cvTitle = gsub('.xls', '', cvTitle)

names(ldfData)
## select significant genes
dfContrast1.sub = ldfData[[1]][ldfData[[1]]$adj.P.Val < 0.01 & abs(ldfData[[1]]$logFC) > 0.5,]
dfContrast2.sub = ldfData[[2]][ldfData[[2]]$adj.P.Val < 0.01 & abs(ldfData[[2]]$logFC) > 0.5,]
dfContrast3.sub = ldfData[[3]][ldfData[[3]]$adj.P.Val < 0.01 & abs(ldfData[[3]]$logFC) > 0.5,]
dfContrast4.sub = ldfData[[4]][ldfData[[4]]$adj.P.Val < 0.01 & abs(ldfData[[4]]$logFC) > 0.5,]

# create a list for overlaps
lVenn = list(rownames(dfContrast1.sub), rownames(dfContrast2.sub), rownames(dfContrast3.sub),
             rownames(dfContrast4.sub))

names(ldfData)
names(lVenn) = cvTitle
cvCommonGenes = Reduce(intersect, lVenn)

## create a binary matrix
cvCommonGenes = unique(do.call(c, lVenn))
mCommonGenes = matrix(NA, nrow=length(cvCommonGenes), ncol=length(lVenn))
colnames(mCommonGenes) = cvTitle

for (i in 1:ncol(mCommonGenes)){
  mCommonGenes[,i] = ldfData[[i]][cvCommonGenes, 'logFC'] 
}
rownames(mCommonGenes) = gsub('X', '', cvCommonGenes)
dim(mCommonGenes)

## create correlation plot for all DEGs
library(NMF)
library(RColorBrewer)

m = (cor(mCommonGenes))
aheatmap(m, annRow = NA, scale = 'none', Rowv = NA, breaks=0.5,
         Colv=NA, cexRow=1, cexCol = 1,  
         col=(brewer.pal(9, 'YlGnBu')))

library(amap)
range(mCommonGenes)
quantile(as.vector(mCommonGenes), 0:20/20)
m = mCommonGenes
m[m < -2] = -2
m[m > 2] = 2
d = Dist(t(m), method='correlation')
hc = hclust(d)
plot(hc)

d = Dist(m, method='correlation')
hcr = hclust(d)
plot(hcr)

aheatmap(m, annRow = NA, scale = 'row', Rowv = hcr, breaks=0,
         Colv=hc, cexRow=1, cexCol = 1, 
         col=rev(brewer.pal(9, 'RdBu')))
