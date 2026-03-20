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
miniSeurat <- subset(seuratObj, celltype %in% 
                         c('Muller glial cells', 'Retinal progenitor cells'))
miniSeurat <- removeRareFeatures(miniSeurat, 1)
miniSeurat <- removeRareFeatures(miniSeurat, 1, 'ATAC')
miniSeurat <- basicDimRed(miniSeurat)
miniSeurat <- jointUMAP(miniSeurat, FALSE, 
                        cutoff=0.005)

DimPlot(miniSeurat, group.by='orig.ident', label=T)
DimPlot(miniSeurat, group.by='celltype')

sce <- as.SingleCellExperiment(miniSeurat)
sce <- slingshot(sce,
                 clusterLabels = 'orig.ident',
                 reducedDim = 'UMAP',
                 start.clus = 'Control',
                 dist.method = 'simple')

miniSeurat <- addLineages(miniSeurat, sce)
mat <- extractLineageMat(sce, 'Lineage1')
df <- points2Seg(mat)
p <- singleLineagePlot(miniSeurat, sce, 'Lineage1', 'orig.ident')

p1 <- featureWes(miniSeurat, 'Lineage1', idClass='orig.ident') + 
    labs(color='Pseudotime') +
    geom_segment(data=df, aes(x=x, y=y, xend=xEnd, yend=yEnd),
                 arrow = arrow(length = unit(0.1, "cm")))

p1 <- featurePlot(miniSeurat, 'Lineage1', 
                  palette=paletteer_c("grDevices::Plasma", 30), 
                  pointSize=0.8) + 
    labs(color='Pseudotime') +
    geom_segment(data=df, aes(x=x, y=y, xend=xEnd, yend=yEnd),
                 arrow = arrow(length = unit(0.1, "cm")))
