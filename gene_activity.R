library(Seurat)
library(Signac)
library(qs2)

source('quality_control.R')
source('tools.R')

miniSeurat <- qs_read('MGCSeurat005.qs2')
DimPlot(miniSeurat, group.by='orig.ident')
activityMat <- GeneActivity(miniSeurat, assay='ATAC', biotypes=NULL)
qs_save(activityMat, 'MGCSeurat005GANullBiotypes.qs2')
dim(activityMat)

activityMat <- GeneActivity(miniSeurat, assay='ATAC')
qs_save(activityMat, 'MGCSeurat005GA.qs2')
activityMat <- qs_read('MGCSeurat005GA.qs2')

miniSeurat[['GeneActivity']] <- CreateAssayObject(activityMat)
DefaultAssay(miniSeurat) <- 'GeneActivity'
miniSeurat <- processSeurat(miniSeurat, assay='GeneActivity', minFeatCells=1, 
                            cutoff=0.01)
p <- DimPlot(miniSeurat, group.by='orig.ident')
devPlot(p)


activityMat <- qs_read('MGCSeurat005GANullBiotypes.qs2')
miniSeurat[['GeneActivity']] <- CreateAssayObject(activityMat)
DefaultAssay(miniSeurat) <- 'GeneActivity'
miniSeurat <- processSeurat(miniSeurat, assay='GeneActivity')
DimPlot(miniSeurat, group.by='orig.ident')
