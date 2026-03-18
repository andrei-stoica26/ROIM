source('load_packages.R')
source('doublets.R')
source('quality_control.R')
source('preprocessing_tools.R')
source('processing_tools.R')
source('trajectory_analysis_tools.R')
source('tools.R')

seuratObj <- qs_read('annotatedSeurat.qs2')
DimPlot(seuratObj, group.by='celltype', label=T, repel=T, label.size=3)

miniSeurat <- subset(seuratObj, celltype == 'Muller glial cells')
miniSeurat <- removeRareFeatures(miniSeurat, 1)
miniSeurat <- removeRareFeatures(miniSeurat, 1, 'ATAC')
miniSeurat <- basicDimRed(miniSeurat)
miniSeurat <- jointUMAP(miniSeurat, FALSE, 
                        cutoff=0.005, 
                        repulsion.strength=0.1, 
                        spread=0.5)
DimPlot(miniSeurat, group.by='orig.ident')

miniSeurat <- subset(seuratObj, celltype == 'Retinal progenitor cells')
miniSeurat <- removeRareFeatures(miniSeurat, 1)
miniSeurat <- removeRareFeatures(miniSeurat, 1, 'ATAC')
miniSeurat <- basicDimRed(miniSeurat)
miniSeurat <- jointUMAP(miniSeurat, FALSE, 0.01)
DimPlot(miniSeurat, group.by='orig.ident')

