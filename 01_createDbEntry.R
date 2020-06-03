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

setwd('dataExternal/remote/raw/')

# list the files
cvFiles = list.files(pattern = 'fq.gz', recursive = T)

# each sample has 2 files 
fSplit = gsub('_[1|2].fq.gz', '', cvFiles)

lFiles = split(cvFiles, fSplit)

setwd(gcswd)
## load the metadata file
dfMeta = read.csv('dataExternal/metaData.csv', header=T, stringsAsFactors = F)
str(dfMeta)
cn = colnames(dfMeta)
# remove any white space
for (i in seq_along(1:ncol(dfMeta))) dfMeta[,cn[i]] = gsub(' ', '', dfMeta[,cn[i]])
dim(dfMeta)
# sanity check
table(as.character(dfMeta$Name) %in% cvFiles)
table(cvFiles %in% as.character(dfMeta$Name))

## order the table in the same sequence as file names
i = match(cvFiles, dfMeta$Name)
dfMeta = dfMeta[i,]
identical(as.character(dfMeta$Name), cvFiles)
dfMeta$fSplit = fSplit

## extract reduced sample table by removing one pair of the fastq file
## this is because the metadata file provided matches each row of the file
## to each fastq file, and there are 2 fastq files per sample
cvSample = unique(dfMeta$fSplit)
i = match(cvSample, dfMeta$fSplit)
dfMeta.sam = dfMeta[i,]
str(dfMeta.sam)
xtabs( ~ fSplit + Lane, data=dfMeta.sam)
xtabs( ~ Technical_Replicate + Lane, data=dfMeta.sam)
xtabs( ~ Technical_Replicate + Biological_Replicate, data=dfMeta.sam)
xtabs( ~ Cell_Line + Induction_state, data=dfMeta.sam)
xtabs( ~ Cell_Line + Biological_Replicate + Technical_Replicate + Induction_state, data=dfMeta.sam)

## create the entry for samples
cSampleCol

dfSamples = data.frame(idProject=g_pid, idData=g_did, title=dfMeta.sam$fSplit,
                       location='rosalind scratch and a copy with john',
                       description= paste('sequencing lane', as.character(dfMeta.sam$Lane),
                                          'group1 is Treatment',
                                          'group2 is cell line',
                                          'group3 is combination of Biological and Technical Replicate',
                                          'samples labelled with NA for treatment are not part of this experiment', sep=';'),
                       group1 = dfMeta.sam$Induction_state,
                       group2= dfMeta.sam$Cell_Line,
                       group3=paste(dfMeta.sam$Biological_Replicate, dfMeta.sam$Technical_Replicate, sep='_'))
# write this data to the database
rownames(dfSamples) = NULL

### NOTE: Do not execute this anymore as entry created
# write this table to database
#dbWriteTable(db, name='Sample', value=dfSamples, append=T, row.names=F)
# get this table again from database with ids added
g_did
dfSamples = dbGetQuery(db, paste('select * from Sample where Sample.idData = 46;'))

# create entries for these files in the database
dbListTables(db)
cn = dbListFields(db, 'File')[-1]
cn
table(names(lFiles) %in% dfSamples$title)
i = match(names(lFiles), dfSamples$title)
dfSamples = dfSamples[i,]
identical(names(lFiles), dfSamples$title)
names(lFiles) = dfSamples$id

# get the names of the samples
temp = lapply(as.character(dfSamples$id), function(x){
  # get the file names
  df = data.frame(name=lFiles[[x]], type='fastq', idSample=dfSamples[as.character(dfSamples$id) == x, 'id'])
  return(df)
})

dfFiles = do.call(rbind, temp)
rownames(dfFiles) = NULL

# sanity check
fn = gsub('_[1|2].fq.gz', '', dfFiles$name)
sapply(dfSamples$id, function(x) {
  table(fn[dfFiles$idSample == x] %in% dfSamples$title[dfSamples$id == x])
})
# write this table to database
## note: do not execute as it is already done
#dbWriteTable(db, name='File', value=dfFiles, append=T, row.names=F)

dbDisconnect(db)
