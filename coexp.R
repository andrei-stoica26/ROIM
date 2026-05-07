library(Seurat)
library(Signac)
library(qs2)
library(scLang)
library(hammers)
library(CSOA)
library(henna)
library(withr)

source('csoa_networks.R')

miniSeurat <- qs_read('miniSeurat.qs2')

res <- buildCSOANetworks(miniSeurat)
res2 <- lapply(res, function(x) 
    breakWeakTies(x, cutoff=0.1, doConnComp=T))
allGenes <- overlapGenes(res2[[1]], unique(res2[[1]]$component))


miniSeurat <- runCSOA(miniSeurat, list(CSOA = allGenes[[4]]))
FeaturePlot(miniSeurat, 'CSOA')

w <- miniSeurat$CSOA

identical(v, w)

View(miniSeurat[[]])
FeaturePlot(miniSeurat, 'CSOA')
plots <- with_seed(1, lapply(res2, networkPlot))
plots[[1]]
View(res2[[1]])
devPlot(plots)

res2 <- lapply(res, function(x) 
    breakWeakTies(x, cutoff=0.45, doConnComp=T))
plots <- with_seed(1, lapply(res2, networkPlot))
