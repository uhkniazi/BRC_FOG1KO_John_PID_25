# File: 16_keywordGO2GeneMapping.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: using specific keywords to search database ids, reverse map the ids to genes
# Date: 27/8/2020

source('header.R')

###########################################################
############ load the count matrix and normalise
###########################################################
library(RMySQL)

db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)

# select the right table using data and project id
q = paste0('select MetaFile.* from MetaFile
           where (MetaFile.idData = 46) AND (MetaFile.comment like "%count%")')
dfSample = dbGetQuery(db, q)
dfSample
n = paste0(dfSample$location, dfSample$name)
load(n)

## load the metadata i.e. covariates
q = paste0('select Sample.* from Sample where Sample.idData = 46')
dfSample = dbGetQuery(db, q)
dim(dfSample)
dfSample = na.omit(dfSample)
dim(dfSample)
head(dfSample)
# close connection after getting data
dbDisconnect(db)

## make count matrix
names(lCounts)
mCounts = do.call(cbind, lCounts)
colnames(mCounts) = names(lCounts)

# sanity check
identical(dfSample$id, as.integer(colnames(mCounts)))

mData = mCounts
dim(mData)

## each lane has a technical replicate, identify those by
## parsing the title field
d = dfSample$title
d = strsplit(d, '_')
fReplicates = factor(sapply(d, function(x) x[5]))
levels(fReplicates)
## sanity check
f = factor(dfSample$group1):factor(dfSample$group2):factor(dfSample$group3)
levels(f)
table(f, fReplicates)
dfSample$fReplicates = factor(fReplicates)
# combine the technical replicates
i = seq_along(1:ncol(mData))
m = tapply(i, dfSample$fReplicates, function(x) {
  return(x)
})

mData = sapply(m, function(x){
  return(rowSums(mCounts[,x]))
})
# get a shorter version of dfSample after adding technical replicates
dfSample.2 = dfSample[sapply(m, function(x) return(x[1])), ]
identical(colnames(mData), as.character(dfSample.2$fReplicates))
dim(dfSample.2)
dfSample.2 = droplevels.data.frame(dfSample.2)

## merge the second set of technical replicates 
## at the library prep level
str(dfSample.2)
f = as.character(dfSample.2$fReplicates)
f = gsub('(.+)-\\d+$', replacement = '\\1', f)
f = factor(f); levels(f)
i = seq_along(1:ncol(mData))
m = tapply(i, f, function(x) {
  return(x)
})

mData = sapply(m, function(x){
  return(rowSums(mData[,x]))
})

# get a shorter version of dfSample after adding technical replicates
dfSample.2 = dfSample.2[sapply(m, function(x) return(x[1])), ]
identical(colnames(mData), levels(f))
dim(dfSample.2)
dfSample.2 = droplevels.data.frame(dfSample.2)

## normalise the data
# drop the rows where average across rows is less than 3
i = rowMeans(mData)
table( i < 3)
mData = mData[!(i< 3),]
dim(mData)

library(DESeq2)
sf = estimateSizeFactorsForMatrix(mData)
mData.norm = sweep(mData, 2, sf, '/')
###########################################################

## load the common binary matrix of DE genes created in earlier results
dfCommonGenes = read.csv('results/commonDEGenes.xls', header=T, row.names=1)
head(dfCommonGenes)

#### To translate the data into a more meaningful biological context and to 
# characterize more thoroughly sets of functionally related genes
# organize the differentially expressed datasets into gene ontology groupings (figures)?
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

## perform analysis on selected combinations
colnames(dfCommonGenes)
bInduced = apply(dfCommonGenes[,c(1, 3, 5)], 1, function(x) all(x == T))
bIndVsNi = apply(dfCommonGenes[,c(2, 4, 6, 7)], 1, function(x) all(x == T))
bNinduced = apply(dfCommonGenes[,c(8, 9, 10)], 1, function(x) all(x == T))

lFilters = list(induced=bInduced, indVSni=bIndVsNi, notInd=bNinduced)
lGO.results = lapply(lFilters, function(flag){
  return(goTest(rownames(dfCommonGenes)[flag]))
})

names(lGO.results)

oFile.go = file('results/GO_groups.csv', 'wt')
temp = sapply(names(lGO.results), function(group){
  p1 = paste('Contrast Comparison group ', group)
  df = summary(lGO.results[[group]])
  p2 = paste(colnames(df), collapse = ',')
  writeLines(p1, oFile.go)
  writeLines(p2, oFile.go)
  sapply(1:10, function(x){
    p3 = gsub(',', replacement = '-', df[x,])
    p3 = paste(p3, collapse=',')
    writeLines(p3, oFile.go)
  })
})

close(oFile.go)

#######################################################
####### find particular types of genes from pathway keyword search
dfGO = AnnotationDbi::select(org.Mm.eg.db, keys = 'Zfpm1', columns = c('GO'), keytype = 'SYMBOL')
dfGO = AnnotationDbi::select(org.Mm.eg.db, keys = rownames(dfCommonGenes), columns = c('GO'), keytype = 'ENTREZID')
dfGO = dfGO[dfGO$ONTOLOGY == 'BP', ]
dfGO = na.omit(dfGO)
dim(dfGO)

library(GO.db)
columns(GO.db)
dfGO = AnnotationDbi::select(GO.db, keys=as.character(unique(dfGO$GO)), columns=columns(GO.db), keytype='GOID')
dim(dfGO)
## keyword search
i = grep('chromatin', dfGO$DEFINITION, ignore.case = T)
length(i)
temp = dfGO[i,]

i = grep('blood cell', dfGO$DEFINITION, ignore.case = T)
length(i)
temp = rbind(temp, dfGO[i,])
dfGO.keyword = temp

### work back to original gene list
dfGO = AnnotationDbi::select(org.Mm.eg.db, keys = dfGO.keyword$GOID, keytype = c('GO'), columns = 'ENTREZID')
dim(dfGO)
head(dfGO)
table(dfGO$ENTREZID %in% rownames(dfCommonGenes))
cvGenes.keyword = dfGO$ENTREZID[(dfGO$ENTREZID %in% rownames(dfCommonGenes))]
cvGenes.keyword = unique(cvGenes.keyword)
length(cvGenes.keyword)
# ## go up to stan section to load the d.bk dataframe
# d.bk = d[as.character(d$split.2) %in% cvGenes.keyword,]
# d.bk = droplevels.data.frame(d.bk)
# library(org.Mm.eg.db)
# df = AnnotationDbi::select(org.Mm.eg.db, keys = as.character(d.bk$split.2), columns = 'SYMBOL', keytype = 'ENTREZID')
# i = match(as.character(d.bk$split.2), df$ENTREZID)
# df = df[i,]
# d.bk$SYMBOL = df$SYMBOL
# identical(as.character(d.bk$split.2), df$ENTREZID)
# head(d.bk)
# d.bk$coef = colMeans(mCoef[,d.bk$cols])
# temp = gsub('(ND:|D:).+', '\\1', as.character(d.bk$fBatch))
# temp = gsub(':', '', temp)
# d.bk$differentiated = factor(temp, levels=c('ND', 'D'))
# 
# temp = gsub('ND:|D:', '', as.character(d.bk$fBatch))
# d.bk$groups = factor(temp)
# 
# library(lattice)
# xyplot(coef ~ differentiated | SYMBOL, groups=groups, data=d.bk, type=c('l', 'p'), scales=list(relation='free', x=list(cex=0.7), y=list(cex=0.7)), 
#        ylab='Model Estimated log Deflections from Intercept', main=list(label='Tooth Development Genes DE expressed at 3 time points in WT', cex=0.8),
#        auto.key=list(columns=3))
# 
# ## cluster the data on trends of expression
# m = split(d.bk$coef, f = d.bk$SYMBOL)
# m = do.call(rbind, m)
# hc = hclust(dist(t(scale(t(m)))))
# plot(hc, main='clustered')
# c = cutree(hc, k = 3)
# table(c)
# 
# ## expand the cluster variable after matching it with symbol
# i = match(as.character(d.bk$SYMBOL), names(c))
# d.bk$SYMBOL.cluster = factor(c[i]):factor(d.bk$SYMBOL)
# 
# xyplot(coef ~ differentiated | SYMBOL.cluster, groups=groups, data=d.bk, type=c('l', 'p'), scales=list(relation='free', x=list(cex=0.7), y=list(cex=0.7)), 
#        ylab='Model Estimated log Deflections from Intercept', main=list(label='Zfpm1 (Fog-1) and Erythrocyte development genes expressed in data', cex=0.8),
#        auto.key=list(columns=3))


# dim(m)
# 
# ## set names for columns
# cn = split(d.bk$fBatch, f = d.bk$SYMBOL)
# colnames(m) = as.character(cn[[1]])
# head(m)
# df = stack(data.frame(m))
# head(df)
# 
# df = data.frame(df, clustering=factor(c), symbols=rownames(m))
# head(df)
# #df$ind = factor(gsub('X', '', df$ind))
# l = factor(df$clustering:df$symbols)
# l = levels(l)
# l = gsub('\\d:', '', l)
# xyplot(values ~ ind | clustering:symbols, data=df, type=c('l', 'p'), scales=list(relation='free', x=list(cex=0.7), y=list(cex=0.7)), 
#        ylab='Model Estimated log Deflections from Intercept', main=list(label='Tooth Development Genes DE expressed at 3 time points in WT', cex=0.8),
#        strip=strip.custom(factor.levels=l))

#######################################################


################################################################################
