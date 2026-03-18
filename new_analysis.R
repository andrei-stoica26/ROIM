library(Seurat)
library(Signac)
library(qs2)

source('quality_control.R')
source('tools.R')

seuratObj <- qs_read('MGCSeurat005.qs2')
DimPlot(seuratObj, group.by='orig.ident')
activityMat <- GeneActivity(seuratObj, assay='ATAC', biotypes=NULL)
qs_save(activityMat, 'MGCSeurat005GANullBiotypes.qs2')
dim(activityMat)

activityMat <- GeneActivity(seuratObj, assay='ATAC')
qs_save(activityMat, 'MGCSeurat005GA.qs2')
seuratObj[['GeneActivity']] <- CreateAssayObject(activityMat)
DefaultAssay(seuratObj) <- 'GeneActivity'
seuratObj <- processSeurat(seuratObj, assay='GeneActivity')
DimPlot(seuratObj, group.by='orig.ident')

activityMat <- qs_read('MGCSeurat005GANullBiotypes.qs2')
seuratObj[['GeneActivity']] <- CreateAssayObject(activityMat)
DefaultAssay(seuratObj) <- 'GeneActivity'
seuratObj <- processSeurat(seuratObj, assay='GeneActivity')
DimPlot(seuratObj, group.by='orig.ident')
