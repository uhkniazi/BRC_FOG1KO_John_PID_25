# File: 18_gseaFilteredMerged.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: gene set enrichment analysis for the filtered and merged results
# Date: 8/3/2021


## set variables and source libraries
source('header.R')
# libraries to load
library(gage)

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

## get the results for the WT induced vs not induced control
ldfData$`results/WT_indVSni.xls` = data.frame(ENTREZID=ldfData[[1]]$ENTREZID,
                       SYMBOL=ldfData[[1]]$SYMBOL,
                       logFC=ldfData[[1]]$WT_indVSni_logFC, 
                       row.names = rownames(ldfData[[1]]))

cvTitle = gsub('results/', '', names(ldfData))
cvTitle = gsub('.xls', '', cvTitle)

## load msig db data
oMsigGS.c2 = readList('~/Data/MetaData/msigdb.v5.2.symbols_mouse.gmt')

## choose a contrast to work with loop through
for (i in 1:length(ldfData)){
  dfContrast = ldfData[[i]]
  # for a contrats of choice create the list
  iContFc = dfContrast$logFC
  ## add enterez ids
  names(iContFc) = as.character(dfContrast$SYMBOL)
  head(iContFc)
  head(dfContrast)
  oGage = gage(iContFc, oMsigGS.c2)
  
  dfGreater = data.frame(oGage$greater)
  #str(dfGreater)
  #i = which(dfGreater$p.val < 0.01)
  #rownames(dfGreater[i,])
  
  dfLess = data.frame(oGage$less)
  #str(dfLess)
  #i = which(dfLess$p.val < 0.01)
  #rownames(dfLess[i,])
  
  write.csv(dfGreater[,c('p.val', 'q.val', 'set.size')], file=paste('results/mergedDownstream/', cvTitle[i], '_upregulated_pathways_mSigDb_c2_curated.xls', sep=''))
  write.csv(dfLess[,c('p.val', 'q.val', 'set.size')], file=paste('results/mergedDownstream/', cvTitle[i], '_downregulated_pathways_mSigDb_c2_curated.xls', sep=''))
  # 
  # ## c5
  # oGage = gage(iContFc, oMsigGS.c5)
  # 
  # dfGreater = data.frame(oGage$greater)
  # dfLess = data.frame(oGage$less)
  # 
  # write.csv(dfGreater[,c('p.val', 'q.val', 'set.size')], file=paste('results/', cvTitle[i], '_upregulated_pathways_mSigDb_c5.xls', sep=''))
  # write.csv(dfLess[,c('p.val', 'q.val', 'set.size')], file=paste('results/', cvTitle[i], '_downregulated_pathways_mSigDb_c5.xls', sep=''))
  # 
  # ## c7
  # oGage = gage(iContFc, oMsigGS.c7)
  # 
  # dfGreater = data.frame(oGage$greater)
  # dfLess = data.frame(oGage$less)
  # 
  # write.csv(dfGreater[,c('p.val', 'q.val', 'set.size')], file=paste('results/', cvTitle[i], '_upregulated_pathways_mSigDb_c7.xls', sep=''))
  # write.csv(dfLess[,c('p.val', 'q.val', 'set.size')], file=paste('results/', cvTitle[i], '_downregulated_pathways_mSigDb_c7.xls', sep=''))
}