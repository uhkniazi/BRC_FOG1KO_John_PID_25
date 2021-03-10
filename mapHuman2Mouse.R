# File: mapHuman2Mouse.R
# Auth: umar.niazi@kcl.ac.uk
# Desc: map human ids and mouse ids to convert human gmt files from msigdb to mouse
# Date: 10/3/2021

library(biomaRt)

## see here for some usage examples
## https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/biomaRt.html#selecting-a-biomart-database-and-dataset
l = listMarts()
l

cMart = l$biomart[1]
oMart = useMart(cMart)
l = listDatasets(oMart)
## find human and mouse datasets
i = grep('human', l$description, ignore.case = T)
i = c(i, grep('mouse', l$description, ignore.case = T))
l[i,]

oHuman = useMart(cMart, dataset="hsapiens_gene_ensembl")
oMouse = useMart(cMart, dataset="mmusculus_gene_ensembl") 

# search for the attributes in the mart
temp = searchAttributes(oHuman, 'entrez')
temp

temp = searchAttributes(oMouse, 'entrez')
temp

dfMap = getLDS(attributes = c('entrezgene_id', 'hgnc_symbol'), filters = 'entrezgene_id', 
               values = '3569', mart=oHuman,
               attributesL = c('entrezgene_id', 'mgi_symbol'),
               martL = oMouse)

colnames(dfMap) = c('human', 'HGNC.symbol', 'mouse', 'MGI.symbol')

dfMap

### read the gmt file using gage 
library(gage)

oMsigGS = readList('dataExternal/gmt/c8.all.v7.2.entrez.gmt')

cvHuman = unique(unlist(oMsigGS))
length(cvHuman)
dfMap = getLDS(attributes = c('entrezgene_id', 'hgnc_symbol'), filters = 'entrezgene_id', 
               values = cvHuman, mart=oHuman,
               attributesL = c('entrezgene_id', 'mgi_symbol'),
               martL = oMouse)

colnames(dfMap) = c('human', 'HGNC.symbol', 'mouse', 'MGI.symbol')
dfMap = na.omit(dfMap)
dfMap = dfMap[!(duplicated(dfMap$human)),]
rownames(dfMap) = dfMap$human

oMsigMapped = lapply(oMsigGS, function(x) {
  return(as.character(na.omit(dfMap[x, 'mouse'])))
})

# write results to file
lGmt.sub = oMsigMapped[!sapply(oMsigMapped, is.null)]

oFile = file('dataExternal/gmt/mouse.c8.all.v7.2.entrez.gmt', 'wt')

x = names(lGmt.sub)
for(y in seq_along(x)){
  p = paste(x[y], 'human_to_mouse_mapping_biomaRt_v_2.46.3', paste(lGmt.sub[[x[y]]], collapse='\t'), sep='\t')
  writeLines(p, oFile)
}

close(oFile)


