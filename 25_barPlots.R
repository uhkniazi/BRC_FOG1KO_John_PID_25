# File: 25_barPlots.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: create barplots for gsea results
# Date: 1/6/2021

source('header.R')

df = read.csv(file.choose(), header=T, stringsAsFactors = F)
head(df)
up = df[,2] + runif(nrow(df), 1e-6, 1e-5)
down = df[,3] + runif(nrow(df), 1e-6, 1e-5)
which(up <= 0.05)
which(down <= 0.05)
dfPlot = data.frame(pv=c(df$induced_KOvsWT.down[which(down <= 0.05)],
                         df$induced_KOvsWT.up[which(up <= 0.05)]))
rownames(dfPlot) = df$Pathway.Label[c(which(down <= 0.05), which(up <= 0.05))]
dfPlot$iCol = c(rep('skyblue', times=6), rep('pink', time=3))
dfPlot$pv = dfPlot$pv + runif(1, 1e-4, 1e-3)
dfPlot$pv = -1 * log(dfPlot$pv)

# iCol = rep('grey', times=length(m))
# iCol[z > 0] = 'pink'
# iCol[z < 0] = 'skyblue'
dfPlot = na.omit(dfPlot)
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
m = dfPlot$pv
names(m) = rownames(dfPlot)
iCol = as.character(dfPlot$iCol)
par(mar=c(5,6,4,2)+0.1)
a = names(m)
barplot(m, names.arg=wrap.labels(a, 15), horiz=T, las=1, cex.names=0.7,
        width=1, main='GSEA pathways: Induced KO VS WT', xlab='-log Pvalue', col=iCol,
        xlim=c(0, 8), xaxt='n')
axis(side = 1, at = c(0:6, 6.9), labels = c(0:6, '>6.9'))

legend('topright', legend = c('Upregulated', 'Downregulated'), 
       fill = c('pink', 'skyblue'), bty = 'n')
