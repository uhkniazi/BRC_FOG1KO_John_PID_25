# File: 11_vennDiagram.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: using the results for each contrast, create venn diagrams and some plots
# Date: 14/8/2020

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

cvTitle[c(1, 3, 5)]
cvTitle[c(2, 4, 6)]
cvTitle[c(8, 9, 10)]

## create a binary matrix
cvCommonGenes = unique(do.call(c, lVenn))
mCommonGenes = matrix(NA, nrow=length(cvCommonGenes), ncol=length(lVenn))
for (i in 1:ncol(mCommonGenes)){
  mCommonGenes[,i] = cvCommonGenes %in% lVenn[[i]]
}
rownames(mCommonGenes) = cvCommonGenes
colnames(mCommonGenes) = names(lVenn)

dfCommonGenes = read.csv('results/commonDEGenes.xls', header=T, stringsAsFactors = F)
colnames(dfCommonGenes)[1] = 'ENTREZID'
data.frame(colnames(dfCommonGenes))

## choose the appropriate combination of contrasts
m = as.matrix(dfCommonGenes[,c(9, 10, 11)])
head(m)
i = which(rowSums(m) >= 2)
length(i)

dfCommonGenes.ind = dfCommonGenes[i,]
dim(dfCommonGenes.ind)

names(ldfData)
names(ldfData)[c(8, 9, 10)]
## extract the relevant information
pv = sapply(ldfData[c(8, 9, 10)], function(x) return(x$adj.P.Val))
rownames(pv) = ldfData[[1]]$ind
pv = pv[as.character(dfCommonGenes.ind$ENTREZID), ]
colnames(pv) = colnames(m)
identical(rownames(pv), as.character(dfCommonGenes.ind$ENTREZID))

fc = sapply(ldfData[c(8, 9, 10)], function(x) return(x$logFC))
rownames(fc) = ldfData[[1]]$ind
fc = fc[as.character(dfCommonGenes.ind$ENTREZID), ]
colnames(fc) = colnames(m)
identical(rownames(fc), rownames(pv))

### apply the filters
f_findNearest2 = function(v){
  d = as.matrix(dist(v))
  diag(d) = NA
  d = arrayInd(which.min(d), dim(d))
  return(as.vector(d))
}

mN = t(apply(fc, 1, f_findNearest2))

## take average of the nearest 2
iAv = sapply(1:nrow(fc), function(x){
  return(mean(fc[x, mN[x,]]))
})

## is the absolute difference between the 2 fc average
## and the 3rd farthest fc greater than 1
i = sapply(1:nrow(fc), function(x){
  return(fc[x,-mN[x,]])
})
iAvd = abs(iAv - i)
fDrop = iAvd > 1
table(fDrop)

iAv3 = rowMeans(fc)
### replace the specific averages that are from 2 fc
iAv3[fDrop] = iAv[fDrop]

## calculate fisher p-values
fp3 = apply(pv, 1, f_fishersMethod)
fp2 = sapply(1:nrow(pv), function(x){
  return(f_fishersMethod(pv[x, mN[x,]]))
})

fpC = fp3
fpC[fDrop] = fp2[fDrop]

identical(as.character(dfCommonGenes.ind$ENTREZID), names(iAv3))
identical(as.character(dfCommonGenes.ind$ENTREZID), names(fpC))

dfCommonGenes.ind$logFC = iAv3
dfCommonGenes.ind$f.pvalue = fpC
dfCommonGenes.ind$Outliers = fDrop
data.frame(colnames(dfCommonGenes.ind))
colnames(fc)
temp = dfCommonGenes.ind[,c(1, 9, 10, 11, 14:17)]

write.csv(temp, file='results/commonDEGenes_not_induced_2of3_logFC_average.xls')
