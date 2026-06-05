library(Seurat)
library(Signac)
library(qs2)
library(ggplot2)
library(hammers)
library(tradeSeq)
library(CSOA)
library(henna)
library(scLang)
library(enResUtils)

source('csoa_networks.R')
source('visualization.R')
source('nnet_subsets.R')

seuratObj <- qs_read('annotatedSeurat.qs2')

seuratObj$celltype <- factor(seuratObj$celltype, 
                             levels=c('Muller glial cells', 
                                      'Retinal progenitor cells',
                                      'Gliogenic progenitors',
                                      'Neurogenic progenitors',
                                      'Proliferative RPC',
                                      'Photoreceptor precursors',
                                      'Rods',
                                      'Bipolar cells',
                                      'Amacrine cells',
                                      'Cones',
                                      'ISL1+ retinal ganglion cells',
                                      'HTR2C+ retinal ganglion cells',
                                      'NEUROD6+ retinal ganglion cells',
                                      'ATF5+ horizontal cells',
                                      'ONECUT1+ horizontal cells',
                                      'RPE cells'))

celltypeCols <- c('Muller glial cells' = rgb(247/255,147/255,30/255), 
                  'Retinal progenitor cells' = rgb(210/255,0/255,0/255), 
                  'Gliogenic progenitors' = rgb(248/255,115/255,106/255), 
                  'Neurogenic progenitors' = rgb(180/255,0/255,150/255),
                  'Proliferative RPC' = rgb(200/255,50/255,15/255), 
                  'Photoreceptor precursors' = rgb(163/255,100/255,140/255), 
                  'Rods' = rgb(169/255,169/255,169/255), 
                  'Bipolar cells' = rgb(163/255,165/255,0/255), 
                  'Amacrine cells' = rgb(157/255,115/255,194/255), 
                  'Cones' = rgb(230/255,134/255,201/255), 
                  'ISL1+ retinal ganglion cells' = rgb(140/255,198/255,63/255),
                  'HTR2C+ retinal ganglion cells' = rgb(192/255,193/255,48/255), 
                  'NEUROD6+ retinal ganglion cells' = rgb(74/255,176/255,91/255), 
                  'ATF5+ horizontal cells' = rgb(97/255,156/255,255/255), 
                  'ONECUT1+ horizontal cells' = rgb(3/255,161/255,198/255), 
                  'RPE cells'= rgb(129/255,70/255,58/255))

p <- DimPlot(seuratObj, group.by='celltype', cols=celltypeCols) +
    theme(
        legend.key.size=unit(0.4, 'cm'),
        legend.text=element_text(size=8)
    )
p <- centerTitle(p, 'Cell type annotation')
devPlot(p)

p <- DimPlot(seuratObj, group.by='orig.ident') +
    theme(
        legend.key.size=unit(0.4, 'cm'),
        legend.text=element_text(size=8)
    )
p <- centerTitle(p, 'Experimental conditions')
devPlot(p)

################################################################################

miniSeurat <- qs_read('miniSeurat.qs2')

allMarkers <- FindAllMarkers(miniSeurat, group.by='orig.ident', only.pos=T, 
                             logfc.threshold=1.5, min.pct=0.2)

miniSeurat$orig.ident <- factor(miniSeurat$orig.ident, levels=c('Control', '0h', '12h', '24h'))
p <- DoHeatmap(miniSeurat, rownames(allMarkers), group.by='orig.ident', angle=0, vjust=0.2)
p <- centerTitle(p, "Top differentially expressed genes")
devPlot(p)

gam <- qs_read('mgcGam.qs2')
res <- associationTest(gam)
res <- res[order(res$waldStat, decreasing=TRUE),]
topGenes <- rownames(subset(res, waldStat > 500))
p <- DoHeatmap(miniSeurat, topGenes, group.by='orig.ident', angle=0, vjust=0.2)
p <- centerTitle(p, "Top genes varying along pseudotime")
devPlot(p)

################################################################################

markersCtr <- FindMarkers(miniSeurat, group.by="orig.ident", ident.1="Control")
p <- volcanoPlot(markersCtr, "Top differentially expressed genes - Control",
                 labLogFCThr = 1,  labPvalThr=1e-50)
devPlot(p)

markers0h <- FindMarkers(miniSeurat, group.by="orig.ident", ident.1="0h")
p <- volcanoPlot(markers0h, "Top differentially expressed genes - 0 hours",
                 labLogFCThr = 1,  labPvalThr=1e-40)
devPlot(p)

markers12h <- FindMarkers(miniSeurat, group.by="orig.ident", ident.1="12h")
p <- volcanoPlot(markers12h, "Top differentially expressed genes - 12 hours",
                 logFCThr=0.8, labLogFCThr = 0.8,  labPvalThr=1e-20)
devPlot(p)

markers24h <- FindMarkers(miniSeurat, group.by="orig.ident", ident.1="24h")
p <- volcanoPlot(markers24h, "Top differentially expressed genes - 24 hours",
                 logFCThr=0.8, labLogFCThr = 0.8,  labPvalThr=1e-15)
devPlot(p)

################################################################################
csoaNets <- qs_read('csoaNets.qs2')
filteredNets <- lapply(csoaNets, function(x) breakWeakTies(x, cutoff=0.05, doConnComp=F))

p <- networkPlot(filteredNets[[4]], "Top overlap genes - Control")
devPlot(p)
p <- geneRadialPlot(filteredNets[[4]], "Top overlap genes - Control", groupLegendTitle='Component')
devPlot(p)

p <- networkPlot(filteredNets[[1]], "Top overlap genes - 0 hours")
devPlot(p)
p <- geneRadialPlot(filteredNets[[1]], "Top overlap genes - 0 hours", groupLegendTitle='Component')
devPlot(p)

p <- networkPlot(filteredNets[[2]], "Top overlap genes - 12 hours")
devPlot(p)
p <- geneRadialPlot(filteredNets[[2]], "Top overlap genes - 12 hours", groupLegendTitle='Component')
devPlot(p)

p <- networkPlot(filteredNets[[3]], "Top overlap genes - 24 hours")
devPlot(p)
p <- geneRadialPlot(filteredNets[[3]], "Top overlap genes - 24 hours", groupLegendTitle='Component')
devPlot(p)

m <- genesER(overlapGenes(filteredNets[[4]]), 'human')
p <- newCnetplot(m, "Top enriched GO terms - Control")
devPlot(p)

m <- genesER(overlapGenes(filteredNets[[1]]), 'human')
p <- newCnetplot(m, "Top enriched GO terms - 0 hours")
devPlot(p)

m <- genesER(overlapGenes(filteredNets[[2]]), 'human')
p <- newCnetplot(m, "Top enriched GO terms - 12 hours")
devPlot(p)

m <- genesER(overlapGenes(filteredNets[[3]]), 'human')
p <- newCnetplot(m, "Top enriched GO terms - 24 hours")
devPlot(p)

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


conds <- c('Control', '0h', '12h', '24h')
cells <- list(seq(19), seq(20), seq(20), seq(20))

invisible(mapply(function(x, y){
    for (i in y){
        p <- netPlots[[x]][[i]]
        p <- centerTitle(p, paste0(x, ' - Metacell ', i))
        devPlot(p)
    }
}, conds, cells, SIMPLIFY=FALSE))
