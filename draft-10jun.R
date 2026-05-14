library(Seurat)
library(qs2)
library(ggplot2)
library(hammers)
library(tradeSeq)

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
################################################################################

miniSeurat <- qs_read('miniSeurat.qs2')
seurats <- SplitObject(miniSeurat, 'orig.ident')

gam <- qs_read('mgcGam.qs2')
res <- associationTest(gam)
res <- res[order(res$waldStat, decreasing=TRUE),]
topGenes <- rownames(subset(res, waldStat > 1000))

obj <- seurats[['0h']]
net0h <- prepareNet(obj, topGenes)

obj <- seurats[['12h']]
net12h <- prepareNet(obj, topGenes)

obj <- seurats[['24h']]
net24h <- prepareNet(obj, topGenes)

obj <- seurats[['Control']]
netCtr <- prepareNet(obj, topGenes)

netSeuratsTradeseq <- list(net0h, net12h, net24h, netCtr)
names(netSeuratsTradeseq) <- names(seurats)

netPlots <- mapply(function(x, y) createNetplots(x, y),
                   netSeuratsTradeseq, 
                   list(seq(20),
                        seq(20),
                        seq(20),
                        seq(20)),
                   SIMPLIFY=FALSE)

netPlots[['0h']][[20]]
netPlots[['12h']][[20]]
netPlots[['24h']][[20]]
netPlots[['Control']][[20]]
dev.off()
