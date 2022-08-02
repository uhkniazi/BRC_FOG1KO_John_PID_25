# small random tasks

# merge the requested files
df.ni = read.csv('results/NotInduced_KOvsWT_merged.xls', header=T, 
                 stringsAsFactors = F, row.names=1)
df.ind = read.csv('results/induced_KOvsWT_merged.xls', header=T, 
                  stringsAsFactors = F, row.names=1)

cvRn = Reduce(intersect, list(rownames(df.ni), rownames(df.ind)))
head(cvRn)

dfCombined = df.ni[cvRn, -5]
colnames(dfCombined)
colnames(dfCombined)[c(3, 4)] = paste0('ni_KOvsWT_', colnames(dfCombined)[c(3, 4)])
colnames(df.ind)
identical(rownames(dfCombined), rownames(df.ind))

dfCombined = cbind(dfCombined, df.ind[cvRn, c(3, 4)])
colnames(dfCombined)[c(7, 8)] = paste0('Ind_KOvsWT_', colnames(dfCombined)[c(7, 8)])
write.csv(dfCombined, 
          file='results/Filtered and merged KO - logFC Pvalue in one sheet.xls')

cvSel = scan(what=character())
f = dfCombined$SYMBOL %in% cvSel
table(f)
dfCombined.sub = dfCombined[f,]
dim(dfCombined.sub)
write.csv(dfCombined.sub, 
          file='results/Filtered and merged KO - logFC Pvalue in one sheet - Qaigen ROS.xls')
