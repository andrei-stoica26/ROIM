library(Seurat)
library(Signac)
library(qs2)
library(scLang)
library(hammers)
library(CSOA)
library(henna)
library(withr)


miniSeurat <- qs_read('miniSeurat.qs2')
res <- buildCSOANetworks(miniSeurat)
res2 <- lapply(res, function(x) 
    breakWeakTies(x, cutoff=0.1, doConnComp=T))
plots <- with_seed(1, lapply(res2, networkPlot))
devPlot(plots)

res2 <- lapply(res, function(x) 
    breakWeakTies(x, cutoff=0.45, doConnComp=T))
plots <- with_seed(1, lapply(res2, networkPlot))
plots[[4]]
