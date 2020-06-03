# File: 01_createDbEntry.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: list the samples and create appropriate db entries
# Date: 2/6/2020


## set variables and source libraries
source('header.R')

## connect to mysql database 
library('RMySQL')

##### connect to mysql database to get samples
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)

# sample and file table
dbGetQuery(db, paste('describe Sample;'))
cSampleCol = dbGetQuery(db, paste('describe Sample;'))$Field[-1]

dbGetQuery(db, paste('describe File;'))
cFileCol = dbGetQuery(db, paste('describe File;'))$Field[-1]

setwd(gcRemoteDir)
setwd('raw/')

# list the files
cvFiles = list.files(pattern = 'fq.gz', recursive = T)

# each sample has 2 files 
fSplit = gsub('_[1|2].fq.gz', '', cvFiles)

lFiles = split(cvFiles, fSplit)

setwd(gcswd)
## load the metadata file
dfMeta = read.csv('dataExternal/metaData.csv', header=T, stringsAsFactors = F)
## perform some acrobatics to format file names
## due to the nature of the metadata file recieved
str(dfMeta)
fn.1 = paste(dfMeta[,1], dfMeta[,2], dfMeta[,3], dfMeta[,4], dfMeta[,5], sep='_')
fn.2 = paste(fn.1, dfMeta[,6], sep='')
fn.2 = gsub('NA', '', fn.2)
fn = paste(fn.2, c('1.fq.gz', '2.fq.gz'), sep='_')

dfMeta = dfMeta[,-c(1:7)]
dfMeta$Name = fn
cn = colnames(dfMeta)
# remove any white space
for (i in seq_along(1:ncol(dfMeta))) dfMeta[,cn[i]] = gsub(' ', '', dfMeta[,cn[i]])

# sanity check
table(as.character(dfMeta$Name) %in% cvFiles)
table(cvFiles %in% as.character(dfMeta$Name))

## order the table in the same sequence as file names
i = match(cvFiles, dfMeta$Name)
dfMeta = dfMeta[i,]
identical(as.character(dfMeta$Name), cvFiles)
dfMeta = dfMeta[,-5]
dfMeta$fSplit = fSplit

## extract reduced sample table by removing one pair of the fastq file
cvSample = unique(dfMeta$fSplit)
i = match(cvSample, dfMeta$fSplit)
dfMeta.sam = dfMeta[i,]
str(dfMeta.sam)
xtabs( ~ Sample_ID + Lesional + Patient_ID, data=dfMeta.sam)
xtabs( ~ Sequencing_Run + Lesional + Patient_ID, data=dfMeta.sam)
xtabs( ~ Sequencing_Run + Sequencing_Lane + Patient_ID, data=dfMeta.sam)

## create the entry for samples
cSampleCol

dfSamples = data.frame(idProject=g_pid, idData=g_did, title=dfMeta.sam$Sample_ID, 
                       description= paste('sequencing lane', as.character(dfMeta.sam$Sequencing_Lane),
                                          'group1 is Treatment',
                                          'group2 is patient batch',
                                          'group3 is sequencing run',
                                          'use replicate pairs as technical replicates identifier', sep=';'),
                       group1 = dfMeta.sam$Lesional, group2= dfMeta.sam$Patient_ID, group3=dfMeta.sam$Sequencing_Run)
# write this data to the database
rownames(dfSamples) = NULL

### NOTE: Do not execute this anymore as entry created
# write this table to database
#dbWriteTable(db, name='Sample', value=dfSamples, append=T, row.names=F)
# get this table again from database with ids added
g_did
dfSamples = dbGetQuery(db, paste('select * from Sample where Sample.idData = 43;'))

# create entries for these files in the database
dbListTables(db)
cn = dbListFields(db, 'File')[-1]
cn
identical(names(lFiles), dfMeta.sam$fSplit)
names(lFiles) = dfSamples$id

# get the names of the samples
temp = lapply(as.character(dfSamples$id), function(x){
  # get the file names
  df = data.frame(name=lFiles[[x]], type='fastq', idSample=dfSamples[as.character(dfSamples$id) == x, 'id'])
  return(df)
})

dfFiles = do.call(rbind, temp)
rownames(dfFiles) = NULL

# write this table to database
## note: do not execute as it is already done
#dbWriteTable(db, name='File', value=dfFiles, append=T, row.names=F)

dbDisconnect(db)
