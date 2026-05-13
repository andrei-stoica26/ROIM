library(Seurat)
library(qs2)
library(ggplot2)
library(hammers)

seuratObj <- qs_read('annotatedSeurat.qs2')
p <- DimPlot(seuratObj, group.by='orig.ident') + ggtitle('Experimental conditions')
devPlot(p)
DimPlot(seuratObj, group.by='celltype', label=TRUE, repel=TRUE, label.size=3) + NoLegend()


miniSeurat <- qs_read('miniSeurat.qs2')
seurats <- SplitObject(miniSeurat, 'orig.ident')

selGenes <- c('NFIA', 'NFIX', 'SOX2','LHX2', 'HMGA1', 'ASCL1', 'SMARCA5', 'YAP1',
              'SOX9', 'STAT3', 'E2F3', 'FOXN4', 'MYB', 'FOXO3', 'SOX5', 'NFIB')