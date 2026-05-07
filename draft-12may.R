library(Seurat)
library(Signac)
library(qs2)
library(slingshot)
library(ggplot2)

source('trajectory_analysis_tools.R')
miniSeurat <- qs_read('miniSeurat.qs2')
sce <- as.SingleCellExperiment(miniSeurat)
sce <- slingshot(sce,
                 clusterLabels = 'orig.ident',
                 reducedDim = 'UMAP',
                 start.clus = 'Control',
                 dist.method = 'simple')

miniSeurat <- addLineages(miniSeurat, sce)
miniSeurat[['Pseudotime']] <- miniSeurat[['Lineage1']]
qs_save(miniSeurat, 'miniSeuratWP.qs2')
