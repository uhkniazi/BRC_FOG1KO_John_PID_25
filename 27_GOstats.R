# File: 27_GOstats.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: over-representation analysis for contrasts of interest
# Date: 19/7/2022

source('header.R')

## load and filter the data
dfData = read.csv('results/paper/NotInduced_KOvsWT_merged.csv', header = T, stringsAsFactors = F)
hist(dfData$f.pvalue)
hist(dfData$logFC)
table(dfData$f.pvalue < 0.01 & abs(dfData$logFC) > 0.5)
fSel = (dfData$f.pvalue < 0.01 & abs(dfData$logFC) > 0.5)
table(fSel)
cv.ni_KOvsWT = dfData$ENTREZID[fSel]

dfData = read.csv('results/paper/induced_KOvsWT_merged.csv', header = T, stringsAsFactors = F)
hist(dfData$f.pvalue)
hist(dfData$logFC)
table(dfData$f.pvalue < 0.01 & abs(dfData$logFC) > 0.5)
fSel = (dfData$f.pvalue < 0.01 & abs(dfData$logFC) > 0.5)
table(fSel)
cv.in_KOvsWT = dfData$ENTREZID[fSel]

library(GOstats)
library(org.Mm.eg.db)

goTest = function(cvSeed, univ = keys(org.Mm.eg.db, 'ENTREZID')){
  ## set up universe background
  dfUniv = AnnotationDbi::select(org.Mm.eg.db, keys = univ, columns = c('GO'), keytype = 'ENTREZID')
  dfUniv = na.omit(dfUniv)
  univ = unique(dfUniv$ENTREZID)
  
  ## make hypergeometric test object for each type, CC, BP and MF
  params = new('GOHyperGParams', geneIds=unique(cvSeed),
               annotation='org.Mm.eg.db',
               universeGeneIds=univ,
               ontology='BP',
               pvalueCutoff= 0.01,
               conditional=FALSE,
               testDirection='over')
  
  oGOStat = tryCatch(hyperGTest(params), error=function(e) NULL)
  return(oGOStat)
}

iSelGenes = as.character(c(cv.ni_KOvsWT, cv.in_KOvsWT))
cvGroups = c(rep('ni', times=length(cv.ni_KOvsWT)),
             rep('in', times=length(cv.in_KOvsWT)))

## perform analysis on most selected genes
lGO.results = lapply(unique(cvGroups), function(group){
  return(goTest(iSelGenes[cvGroups == group]))
})

names(lGO.results) = unique(cvGroups) 

df = summary(lGO.results[['ni']])
grep('cholesterol', df$Term)
write.csv(df, file='results/GO_ni_KOvsWT.csv')

df = summary(lGO.results[['in']])
grep('cholesterol', df$Term)
write.csv(df, file='results/GO_in_KOvsWT.csv')

# oFile.go = file('results/GO_groups.csv', 'wt')
# temp = sapply(as.character(iGroups), function(group){
#   p1 = paste('Contrast Comparison group ', group)
#   df = summary(lGO.results[[group]])
#   p2 = paste(colnames(df), collapse = ',')
#   writeLines(p1, oFile.go)
#   writeLines(p2, oFile.go)
#   sapply(1:10, function(x){
#     p3 = gsub(',', replacement = '-', df[x,])
#     p3 = paste(p3, collapse=',')
#     writeLines(p3, oFile.go)
#   })
# })
# 
# close(oFile.go)

######## produce figures from selected pathways
df = read.csv(file='results/GO_ni_KOvsWT.csv', stringsAsFactors = F, header=T, row.names = 1)
cvScan = scan(what = character(), sep='\n')

table(df$Term %in% cvScan)
df = df[df$Term %in% cvScan, ]
dim(df)

# select a few pathways ~ 8
df = df[order(df$Pvalue, decreasing = F),]
df = df[1:8,]
## https://stackoverflow.com/questions/20241065/r-barplot-wrapping-long-text-labels
# Core wrapping function
wrap.it <- function(x, len)
{ 
  sapply(x, function(y) paste(strwrap(y, len), 
                              collapse = "\n"), 
         USE.NAMES = FALSE)
}


# Call this function with a list or vector
wrap.labels <- function(x, len)
{
  if (is.list(x))
  {
    lapply(x, wrap.it, len)
  } else {
    wrap.it(x, len)
  }
}

m = -1* log(df$Pvalue)
names(m) = df$Term
iCol = 'grey'
par(mar=c(5,6,4,2)+0.1)
a = names(m)
range(m)
summary(m)
m[m>19] = 19
barplot(m, names.arg=wrap.labels(a, 15), horiz=T, las=1, cex.names=0.7,
        width=1, main='GO Over-represented Pathways: Not Induced KO VS WT', xlab='-log Pvalue', col=iCol,
        xlim=c(0, 19), xaxt='n', sub='Cholesterol')
axis(side = 1, at = c(0:18, 19), labels = c(0:18, '>19'))

# legend('topright', legend = c('Upregulated', 'Downregulated'), 
#        fill = c('pink', 'skyblue'), bty = 'n')
