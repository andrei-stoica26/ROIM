library(Seurat)
library(Signac)
library(qs2)
library(scLang)
library(hammers)
library(CSOA)


miniSeurat <- qs_read('miniSeurat.qs2')
res <- buildCSOANetworks()
