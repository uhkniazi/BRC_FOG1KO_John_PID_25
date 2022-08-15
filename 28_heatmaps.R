# File: 28_heatmaps.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: create heatmaps from selected gene lists
# Date: 15/8/2022

source('header.R')

## load and filter the data
dfData.niKOvsWT = read.csv('results/paper/NotInduced_KOvsWT_merged.csv', header = T, stringsAsFactors = F)
dfData.inKOvsWT = read.csv('results/paper/induced_KOvsWT_merged.csv', header = T, stringsAsFactors = F)
dfData.KOIndvsNi = read.csv('results/paper/KO_indVSni_merged.csv', header = T, stringsAsFactors = F)
identical(dfData.inKOvsWT$ENTREZID, dfData.niKOvsWT$ENTREZID)
identical(dfData.KOIndvsNi$ENTREZID, dfData.niKOvsWT$ENTREZID)

## select matching genes
cvSelSymbol = scan(what=character(), sep = '\n')
table(cvSelSymbol %in% dfData.inKOvsWT$SYMBOL)
cvSelSymbol = cvSelSymbol[cvSelSymbol %in% dfData.inKOvsWT$SYMBOL]
i = match(cvSelSymbol, dfData.inKOvsWT$SYMBOL)
dfData.inKOvsWT = dfData.inKOvsWT[i,]
dfData.niKOvsWT = dfData.niKOvsWT[i,]
dfData.KOIndvsNi = dfData.KOIndvsNi[i,]
identical(cvSelSymbol, dfData.KOIndvsNi$SYMBOL)

## create a heatmap for the logFC
library(NMF)
library(RColorBrewer)

dfData = data.frame(inKOvsWT = dfData.inKOvsWT$logFC,
                    niKOvsWT = dfData.niKOvsWT$logFC,
                    KOindVSni = dfData.KOIndvsNi$logFC,
                    WTindVSni = dfData.inKOvsWT$WT_indVSni_logFC, row.names = dfData.inKOvsWT$SYMBOL)
mCommonGenes = as.matrix(dfData)

m = (cor(mCommonGenes))
aheatmap(m, annRow = NA, scale = 'none', Rowv = NA,# breaks=0.4,
         Colv=NA, cexRow=1, cexCol = 1,  
         col=(brewer.pal(9, 'YlGnBu')))

## create a heatmap of the logFC instead of correlation
library(amap)
range(mCommonGenes)
quantile(as.vector(mCommonGenes), 0:20/20)
m = mCommonGenes[,c(1:2)]
m[m < -2.5] = -2.5
m[m > 3] = 3
d = Dist(t(m), method='correlation')
hc = hclust(d)
plot(hc)

d = Dist(m, method='correlation')
hcr = hclust(d)
plot(hcr)

aheatmap(m, annRow = NA, scale = 'row', Rowv = hcr, breaks=0,
         Colv=hc, cexRow=1, cexCol = 1, 
         col=rev(brewer.pal(9, 'RdBu')))

aheatmap(m, annRow = NA, scale = 'none', Rowv = hcr, #breaks=0,
         Colv=hc, cexRow=1, cexCol = 1, 
         col=rev(brewer.pal(9, 'RdBu')))



