library(Seurat)
library(qs2)
library(ggplot2)
library(hammers)

source('nnet_subsets.R')

seuratObj <- qs_read('annotatedSeurat.qs2')
p <- DimPlot(seuratObj, group.by='orig.ident') + ggtitle('Experimental conditions')
devPlot(p)
DimPlot(seuratObj, group.by='celltype', label=TRUE, repel=TRUE, label.size=3) + NoLegend()


miniSeurat <- qs_read('miniSeurat.qs2')
seurats <- SplitObject(miniSeurat, 'orig.ident')

names(seurats)

selGenes <- c('NFIA', 'NFIX', 'SOX2','LHX2', 'HMGA1', 'ASCL1', 'SMARCA5', 'YAP1',
              'SOX9', 'STAT3', 'E2F3', 'FOXN4', 'MYB', 'FOXO3', 'SOX5', 'NFIB')

netSeurats <- lapply(seurats, function(obj) prepareNet(obj, selGenes))

qs_save(netSeurats, 'netSeurats.qs2')
netPlots <- mapply(function(x, y) createNetplots(x, y),
                   netSeurats, 
                   list(seq(20),
                        c(1, 2, 5, 7, 8, 10, 14, 16, 19),
                        seq(20),
                        c(3, 4, 6, 14, 20)),
                   SIMPLIFY=FALSE)

