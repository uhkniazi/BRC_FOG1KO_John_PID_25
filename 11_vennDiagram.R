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

cvTitle[c(1, 3, 5)]
cvTitle[c(2, 4, 6)]
cvTitle[c(8, 9, 10)]

dfCommonGenes.ind = ldfData$`results//DEAnalysis_ind:C8VSind:MELWT.xls`
dfCommonGenes.ind = dfCommonGenes.ind[,c(1, 10)]
colnames(dfCommonGenes.ind)[1] = 'ENTREZID'
rownames(dfCommonGenes.ind) = dfCommonGenes.ind$ENTREZID
names(ldfData)
names(ldfData)[c(1, 3, 5)]
## extract the relevant information
pv = sapply(ldfData[c(1, 3, 5)], function(x) return(x$adj.P.Val))
rownames(pv) = ldfData[[1]]$ind
pv = pv[as.character(dfCommonGenes.ind$ENTREZID), ]
colnames(pv) = cvTitle[c(1, 3, 5)]
identical(rownames(pv), as.character(dfCommonGenes.ind$ENTREZID))

fc = sapply(ldfData[c(1, 3, 5)], function(x) return(x$logFC))
rownames(fc) = ldfData[[1]]$ind
fc = fc[as.character(dfCommonGenes.ind$ENTREZID), ]
colnames(fc) = cvTitle[c(1, 3, 5)]
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

## get the log fold change and adjusted p-value for wt
temp = ldfData$`results//DEAnalysis_ind:MELWTVSn.i:MELWT.xls`
dfCommonGenes.ind$IndWtVsNiWt_logFC = temp$logFC
dfCommonGenes.ind$IndWtVsNiWt_adj.P.Val = temp$adj.P.Val

write.csv(temp, file='results/crossSectional_T2_KOmerged.xls')
