library(Seurat)
library(Signac)
library(qs2)
library(ggplot2)
library(hammers)
library(tradeSeq)
library(CSOA)
library(henna)
library(scLang)
source('nnet_subsets.R')

seuratObj <- qs_read('annotatedSeurat.qs2')
p <- DimPlot(seuratObj, group.by='orig.ident') + ggtitle('Experimental conditions')
devPlot(p)
DimPlot(seuratObj, group.by='celltype', label=TRUE, repel=TRUE, label.size=3) + NoLegend()


miniSeurat <- qs_read('miniSeurat.qs2')
allMarkers <- FindAllMarkers(miniSeurat, group.by='orig.ident', only.pos=T, 
                             logfc.threshold=1.5, min.pct=0.2)

miniSeurat$orig.ident <- factor(miniSeurat$orig.ident, levels=c('Control', '0h', '12h', '24h'))
p <- DoHeatmap(miniSeurat, rownames(allMarkers), group.by='orig.ident', angle=0, vjust=0.2)
devPlot(p)

featureWes(miniSeurat, 'FTH1', idClass='orig.ident')

gam <- qs_read('mgcGam.qs2')
res <- associationTest(gam)
res <- res[order(res$waldStat, decreasing=TRUE),]
topGenes <- rownames(subset(res, waldStat > 500))
p <- DoHeatmap(miniSeurat, topGenes, group.by='orig.ident', angle=0, vjust=0.2)
devPlot(p)

################################################################################

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

set.seed(123)

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

rm(net0h)
rm(net12h)
rm(net24h)
rm(netCtr)

netPlots <- mapply(function(x, y) createNetplots(x, y),
                   netSeuratsTradeseq, 
                   list(seq(20),
                        seq(20),
                        seq(20),
                        c(seq(16), seq(18, 20))),
                   SIMPLIFY=FALSE)

cond <- 'Control'

conds <- c('Control', '0h', '12h', '24h')
cells <- list(seq(19), seq(20), seq(20), seq(20))

invisible(mapply(function(x, y){
    for (i in y){
        p <- netPlots[[x]][[i]]
        p <- centerTitle(p, paste0(x, ' - Metacell ', i))
        devPlot(p)
    }
}, conds, cells, SIMPLIFY=FALSE))

featureWes(miniSeurat, 'TTR', idClass='orig.ident')


################################################################################

source('csoa_networks.R')

miniSeurat <- qs_read('miniSeurat.qs2')
a <- FindMarkers(miniSeurat, group.by="orig.ident", ident.1="Control")
p <- volcanoPlot(a, labLogFCThr = 1,  labPvalThr=1e-50)
p
a <- FindMarkers(miniSeurat, group.by="orig.ident", ident.1="0h")
p <- volcanoPlot(a)
a <- FindMarkers(miniSeurat, group.by="orig.ident", ident.1="12h")
p <- volcanoPlot(a)
a <- FindMarkers(miniSeurat, group.by="orig.ident", ident.1="24h")
p <- volcanoPlot(a)
devPlot(p)

csoaNets <- buildCSOANetworks(miniSeurat)
qs_save(csoaNets,'csoaNets.qs2')


filteredNets <- lapply(csoaNets, function(x) breakWeakTies(x, cutoff=0.05, doConnComp=F))


geneRadialPlot(filteredNets[[3]], groupLegendTitle='Component')

?geneRadialPlot
