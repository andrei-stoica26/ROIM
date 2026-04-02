metrics <- c('nCount_RNA', 'nCount_ATAC', 'percent.mt', 
             'blacklistFrac', 'TSS.enrichment', 'nucleosome_signal',
             'nFeature_RNA', 'nFeature_ATAC')
cutoffs <- c(1000, 1200, 15, 0.01, 3, 2, 500, 1000)
ids <- c('0', '12', '24', 'ctr')

for (i in seq_along(metrics)){
    metric <- metrics[i]
    cutoff <- cutoffs[i]
    for (j in seq_along(ids)){
        fileName <- paste0(metric, '-', ids[j], '.pdf')
        p <- cutoffVlnPlot(seurats[[j]], metric, cutoff)
        pdf(fileName, width = 8, height = 8)
        print(p)
        dev.off()
    }
}

