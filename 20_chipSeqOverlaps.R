# File: 20_chipSeqOverlaps.R
# Auth: umar.niazi@kcl.ac.uk
# Desc: integrating the merged results with legacy chip-seq data
# Date: 23/3/2021

source('header.R')

lFiles = c('results/NotInduced_KOvsWT_merged.xls',
           'results/induced_KOvsWT_merged.xls',
           'results/KO_indVSni_merged.xls')

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

ldfData$`results/WT_indVSni.xls` = read.csv('results/DEAnalysis_ind:MELWTVSn.i:MELWT.xls',
                                            header=T, row.names=1, stringsAsFactors = F)

cvTitle = gsub('results/', '', names(ldfData))
cvTitle = gsub('.xls', '', cvTitle)

######## select genes based on fisher pvalues and fold changes
dfContrast1.sub = ldfData[[1]][ldfData[[1]]$f.pvalue < 0.01 & abs(ldfData[[1]]$logFC) > 0.5,]
dfContrast2.sub = ldfData[[2]][ldfData[[2]]$f.pvalue < 0.01 & abs(ldfData[[2]]$logFC) > 0.5,]
dfContrast3.sub = ldfData[[3]][ldfData[[3]]$f.pvalue < 0.01 & abs(ldfData[[3]]$logFC) > 0.5,]
dfContrast4.sub = ldfData[[4]][ldfData[[4]]$adj.P.Val < 0.01 & abs(ldfData[[4]]$logFC) > 0.5,]

####### create a common genes binary matrix using the selected
####### chip-seq dataset
lVenn = list(rownames(dfContrast1.sub),
             rownames(dfContrast2.sub),
             rownames(dfContrast3.sub),
             as.character(dfContrast4.sub$ind))
names(ldfData)
cvTitle
names(lVenn) = cvTitle

## create a binary matrix
# read in the appropriate data file
dfChip = read.csv('dataExternal/GATA-1_ChIP-seq_in_Fetal_liver.csv', header=T, stringsAsFactors = F)
head(dfChip)
cvSymbol = unique(dfChip$GATA.1_targets_T119.ve.1)
# convert symbols to enterez id
library(org.Mm.eg.db)
df = AnnotationDbi::select(org.Mm.eg.db, keys = cvSymbol, keytype = 'SYMBOL', columns = 'ENTREZID')
dim(df)
dim(na.omit(df))
df = na.omit(df)
table(cvSymbol %in% df$SYMBOL)
cvSymbol = cvSymbol[cvSymbol %in% df$SYMBOL]
i = match(cvSymbol, df$SYMBOL)
df = df[i,]
identical(cvSymbol, df$SYMBOL)
cvEnterez = df$ENTREZID

cvCommonGenes = unique(cvEnterez)
mCommonGenes = matrix(NA, nrow=length(cvCommonGenes), ncol=length(lVenn))
for (i in 1:ncol(mCommonGenes)){
  mCommonGenes[,i] = cvCommonGenes %in% lVenn[[i]]
}
rownames(mCommonGenes) = cvCommonGenes
colnames(mCommonGenes) = names(lVenn)

# create groups in the data based on permutation with repetition: 2 ^ ncol 
# https://www.mathsisfun.com/combinatorics/combinations-permutations.html
mCommonGenes.grp = mCommonGenes
set.seed(123)
dm = dist(mCommonGenes.grp, method='binary')
hc = hclust(dm)
plot(hc)
# cut the tree at the bottom to create groups
cp = cutree(hc, h = 0.1)
# sanity checks
table(cp)
length(cp)
length(unique(cp))

mCommonGenes.grp = cbind(mCommonGenes.grp, cp)

### print and observe this table and select the groups you are interested in
temp = mCommonGenes.grp
temp = (temp[!duplicated(cp),])
temp2 = cbind(temp, table(cp))
rownames(temp2) = NULL
print(temp2)

## extract symbols 
df = AnnotationDbi::select(org.Mm.eg.db, keys = rownames(mCommonGenes), columns = 'SYMBOL', keytype = 'ENTREZID')
dim(df)
dim(na.omit(df))
identical(rownames(mCommonGenes), df$ENTREZID)

## put these results together
dfCommonGenes = data.frame(mCommonGenes, Common=rowSums(mCommonGenes), groups=cp, Symbol=df$SYMBOL)

write.csv(dfCommonGenes, file='results/GATA-1_ChIP-seq_in_Fetal_liver_T119_positive_commonDEGs.xls')
