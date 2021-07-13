# File: 26_bamFilterAndMerge.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: filter a set of bam files for a region of interest and merge files
# Date: 12/7/2021

source('header.R')

library(RMySQL)
require(Rsamtools)

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(GenomicRanges)

# get the exons into GRangesList object
oGRLgenes = exonsBy(TxDb.Mmusculus.UCSC.mm10.knownGene, by = 'gene')
oGRquery = range(oGRLgenes[['13018']])
# expand the query region size
oGRquery = resize(oGRquery, width = width(oGRquery)+10000, fix='center')
# sanity check if region overlaps
table(overlapsAny(oGRLgenes[['13018']], oGRquery))

##### connect to mysql database to get samples
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)
dbListFields(db, 'File')
dbListFields(db, 'Sample')
# get the query
g_did
q = paste0('select Sample.id as sid, Sample.group1, Sample.group2, Sample.title, File.* from Sample, File
           where (Sample.idData = 46) AND (File.idSample = Sample.id AND File.type like "quality 10 sorted bam")')
dfSample = dbGetQuery(db, q)
nrow(dfSample)
dfSample
# close connection after getting data
dbDisconnect(db)
# remove any whitespace from the names
dfSample$name = gsub(" ", "", dfSample$name, fixed = T)
dfSample$title = gsub(" ", "", dfSample$title, fixed = T)
i = grep('q10_sort.bam', dfSample$name)
dfSample = dfSample[i,]
#### set working directory to appropriate location with bam files
setwd('dataExternal/remote/Aligned/')
csFiles = list.files('.', pattern = '*.bam$', recursive = T)
# check if these files match the file names in database
table(dfSample$name %in% csFiles)

dfSample = dfSample[dfSample$name %in% csFiles, ]
dim(dfSample)
#csFiles = dfSample$name
## create factors for grouping files
f1 = factor(dfSample$group1)
f2 = rep('WT', times=nrow(dfSample))
f2[!dfSample$group2 == 'MELWT'] = 'KO'
f2 = factor(f2, levels = c('WT', 'KO'))
fGroups = factor(f2:f1)
levels(fGroups)


dir.create('output')
dir.create('merged')
# ### testing reading one file
# param = ScanBamParam(which = oGRquery)
# d = filterBam(dfSample$name[2], destination = 'output/testing_2.bam', param=param)
# dir.create('merged')
# d2 = mergeBam(c('output/testing.bam', 'output/testing_2.bam'),
#               'merged/merged.bam', indexDestination=T, overwrite=T)
## filter the bam files 
param = ScanBamParam(which = oGRquery)
d = sapply(dfSample$name, function(x){
  filterBam(x, destination = paste0('output/', x, '.bam'), param=param)
})

## split the file list by groups of interest
lFiles = split(paste0(dfSample$name, '.bam'), fGroups)
lRet = sapply(names(lFiles), function(x){
  mergeBam(paste0('output/', lFiles[[x]]), 
           paste0('merged/', x, '.bam'), indexDestination=T)
})
setwd(gcswd)
