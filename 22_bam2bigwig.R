# File: 22_bam2bigwig.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: use the bam file to calculate coverage and write a bigwig file output.
#       some ideas taken from following places:
#       https://research.stowers.org/cws/CompGenomics/Tutorial/bam2bigwig.html
#       https://github.com/Przemol/seqplots/issues/5
#       https://github.com/uhkniazi/CBamQuality/blob/437fad3f760a4297322b7d0654fee94e210a7b68/CBamQuality.R#L57
# Date: 7/4/2021

source('header.R')

library(RMySQL)
require(Rsamtools)
require(GenomicAlignments)
require(rtracklayer)

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
dir.create('output')

bam2bw = function(bf, out){
  # load the bam file as GAlignment object
  ga = readGAlignments(bf)
  cat('read galignments from ', bf, '\n')
  cov = coverage(ga)
  out = paste(out, ".bw", sep="")
  cat(paste("exporting to bigwig", out, "\n", sep=" "))
  export.bw(cov, con=out)
}

dfSample$out = paste0('output/', paste(dfSample$group1, dfSample$group2, 
                     dfSample$idSample, dfSample$title, sep = '_'))

sapply(seq_along(1:nrow(dfSample)), function(x) bam2bw(dfSample$name[x], dfSample$out[x]))


