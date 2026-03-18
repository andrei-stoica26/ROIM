source('load_packages.R')
source('doublets.R')
source('quality_control.R')
source('preprocessing_tools.R')
source('processing_tools.R')
source('trajectory_analysis_tools.R')
source('tools.R')

#############################Work in progress###################################

seuratObj <- qs_read('annotatedSeurat.qs2')
DimPlot(seuratObj, group.by='celltype', label=T, repel=T, label.size=3)

miniSeurat <- subset(seuratObj, celltype == 'Muller glial cells')
miniSeurat <- removeRareFeatures(miniSeurat, 1)
miniSeurat <- removeRareFeatures(miniSeurat, 1, 'ATAC')
miniSeurat <- basicDimRed(miniSeurat)
miniSeurat <- jointUMAP(miniSeurat, FALSE, 
                        cutoff=0.001, 
                        repulsion.strength=0.001,
                        spread=0.3)



