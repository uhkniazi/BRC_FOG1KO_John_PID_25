# File: 23_GSEA2GeneMapping.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: using specific keywords to search database ids, reverse map the ids to genes
# Date: 12/4/2021

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

## normalise the data
# drop the rows where average across rows is less than 3
i = rowMeans(mData)
table( i < 3)
# FALSE  TRUE 
# 12364 12230 
mData = mData[!(i< 3),]
dim(mData)
# [1] 12364     48

library(DESeq2)
sf = estimateSizeFactorsForMatrix(mData)
mData.norm = sweep(mData, 2, sf, '/')

identical(colnames(mData.norm), as.character(dfSample.2$fReplicates))
###########################################################
# paste pathway names to search the mapping genes
cvPathways = scan(what=character())

library(gage)
oMsigGS.c5 = readList('dataExternal/gmt/mouse.c5.all.v7.2.entrez.gmt')

f = names(oMsigGS.c5) %in% cvPathways
oMsigGS.c5 = oMsigGS.c5[f]

# create a data frame and add symbols for genes
i = names(oMsigGS.c5)
dfOut = data.frame(pathway=NA, entrezid=NA)
for(x in seq_along(i)) {
  df = data.frame(i[x], (oMsigGS.c5[i[x]]))
  colnames(df) = colnames(dfOut)
  dfOut = rbind(dfOut, df)
}
dim(dfOut)
dfOut = na.omit(dfOut)
dim(dfOut)

library(org.Mm.eg.db)
df = AnnotationDbi::select(org.Mm.eg.db, keys = dfOut$entrezid,
                           columns = c('SYMBOL'), keytype = 'ENTREZID')
identical(df$ENTREZID, dfOut$entrezid)
dfOut$Symbol = df$SYMBOL

#######################################################
####### find particular types of genes from pathway keyword search
####### or do a reverse lookup after selecting a gene of interest
library(GO.db)
columns(GO.db)

## search the full GO DB OR
dfGO = AnnotationDbi::select(GO.db, keys=keys(GO.db), columns=columns(GO.db), keytype='GOID')
dim(dfGO)

## OR do a reverse lookup using a seed gene list
dfGO = AnnotationDbi::select(org.Mm.eg.db, keys = '27409', columns = c('GO'), keytype = 'ENTREZID')
#dfGO = dfGO[dfGO$ONTOLOGY == 'BP', ]
dfGO = na.omit(dfGO)
dim(dfGO)
dfGO = AnnotationDbi::select(GO.db, keys=as.character(unique(dfGO$GO)), columns=columns(GO.db), keytype='GOID')
dim(dfGO)
dfGO.keyword = dfGO

## do a keyword search to reduce number of pathways to selected type
i = grep('Cholesterol|lipid', dfGO$TERM, ignore.case = T)
length(i)
dfGO.keyword = dfGO[i,]
# temp = dfGO[i,]
# i = grep('metaboli', temp$TERM, ignore.case = T)
# length(i)
# temp$TERM[i]
# dfGO.keyword = temp[i,]

### work back to original gene list
dfGO = AnnotationDbi::select(org.Mm.eg.db, keys = dfGO.keyword$GOID, keytype = c('GO'), columns = c('ENTREZID', 'SYMBOL', 'GENENAME'))
dim(dfGO)
head(dfGO)
table(dfGO$ENTREZID %in% rownames(mData.norm))
dfGO = dfGO[dfGO$ENTREZID %in% rownames(mData.norm), ]

table(rownames(mData.norm) %in% dfGO$ENTREZID)

write.csv(dfGO, file='results/cholensterol_lipid.csv')
write.csv(dfOut, file='results/GSEA2Genes.csv')
### to add here
# 
# cvGenes.keyword = dfGO$ENTREZID[(dfGO$ENTREZID %in% rownames(dfCommonGenes))]
# cvGenes.keyword = unique(cvGenes.keyword)
# length(cvGenes.keyword)
# # ## go up to stan section to load the d.bk dataframe
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
