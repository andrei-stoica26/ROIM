#############################Work in progress###################################

seuratObj <- qs_read('annotatedSeurat.qs2')
p <- DimPlot(seuratObj, group.by='celltype', label=T, repel=T, label.size=3) + 
    NoLegend() + ggtitle('Cell type annotation')
devPlot(p)

miniSeurat <- subset(seuratObj, celltype == 'Muller glial cells')
miniSeurat <- removeRareFeatures(miniSeurat, 1)
miniSeurat <- removeRareFeatures(miniSeurat, 1, 'ATAC')
miniSeurat <- basicDimRed(miniSeurat)
miniSeurat <- jointUMAP(miniSeurat, FALSE, 
                        cutoff=0.001, repulsion.strength=0.3, spread=0.3)

p <- DimPlot(miniSeurat, group.by='orig.ident', label=T, repel=T, label.size=3)
devPlot(p)

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
devPlot(p)

p <- featurePlot(miniSeurat, 'Lineage1', 
                 palette=paletteer_c("grDevices::Plasma", 30), 
                 pointSize=0.8) + 
    labs(color='Pseudotime') +
    geom_segment(data=df, aes(x=x, y=y, xend=xEnd, yend=yEnd),
                 arrow = arrow(length = unit(0.1, "cm")))

devPlot(p)

################################################################################
sce <- as.SingleCellExperiment(miniSeurat)
sce <- slingshot(sce,
                 clusterLabels = 'orig.ident',
                 reducedDim = 'UMAP',
                 start.clus = 'Control',
                 end.clus = '0h',
                 dist.method = 'simple')
miniSeurat <- addLineages(miniSeurat, sce)
mat <- extractLineageMat(sce, 'Lineage1')
df <- points2Seg(mat)
p <- singleLineagePlot(miniSeurat, sce, 'Lineage1', 'orig.ident')
devPlot(p)

p <- featurePlot(miniSeurat, 'Lineage1', 
                 palette=paletteer_c("grDevices::Plasma", 30), 
                 pointSize=0.8) + 
    labs(color='Pseudotime') +
    geom_segment(data=df, aes(x=x, y=y, xend=xEnd, yend=yEnd),
                 arrow = arrow(length = unit(0.1, "cm")))

devPlot(p)

################################################################################

seuratObj[['ATAC']] <- NULL
miniSeurat <- subset(seuratObj, celltype=='Muller glial cells')
miniSeurat <- processSeurat(miniSeurat, minFeatCells=1, cutoff=0.005)
p <- DimPlot(miniSeurat, group.by='orig.ident', label = TRUE, label.size=3) + NoLegend()
miniSeurat <- FindNeighbors(miniSeurat, reduction = "umap", dims = 1:2)
miniSeurat <- FindClusters(miniSeurat, resolution=0.1)

miniSeurat <- subset(miniSeurat, seurat_clusters != 3)
miniSeurat <- processSeurat(miniSeurat, minFeatCells=1, cutoff=0.005)
p <- DimPlot(miniSeurat, group.by='orig.ident')
devPlot(p)

a <- FindMarkers(miniSeurat, ident.1=3, only.pos=TRUE, 
                 min.diff.pct=0.3, 
                 logfc.threshold=2)

Idents(seuratObj) <- 'celltype'
genes <- rownames(a)[seq(4)]

invisible(lapply(genes, function(gene){
    p <- FeaturePlot(miniSeurat, gene)
    devPlot(p)
    p <- FeaturePlot(seuratObj, gene, label=TRUE, repel=TRUE, label.size=3)
    devPlot(p)
}))

m <- genesER(rownames(a), 'human')
p <- newCnetplot(m)
devPlot(p)

a <- FindMarkers(miniSeurat, ident.1=c(0, 1, 2), only.pos=TRUE, 
                 min.diff.pct=0.2, 
                 logfc.threshold=0.1)

p <- VlnPlot(miniSeurat, 'nCount_RNA')
devPlot(p)
p <- VlnPlot(miniSeurat, 'nFeature_RNA')
devPlot(p)

p <- VlnPlot(seuratObj, 'nCount_RNA', group.by='celltype') + NoLegend()
devPlot(p)
p <- VlnPlot(seuratObj, 'nFeature_RNA', , group.by='celltype') + NoLegend()
devPlot(p)

################################################################################

DimPlot(miniSeurat, label=TRUE)
cells <- colnames(subset(miniSeurat, seurat_clusters != 3))

seuratObj <- qs_read('annotatedSeurat.qs2')

miniSeurat <- subset(seuratObj, cells=cells)
miniSeurat <- removeRareFeatures(miniSeurat, 1)
miniSeurat <- removeRareFeatures(miniSeurat, 1, 'ATAC')
miniSeurat <- basicDimRed(miniSeurat)
miniSeurat <- jointUMAP(miniSeurat, FALSE, cutoff=0.005)

p <- DimPlot(miniSeurat, label=TRUE, group.by='orig.ident')
devPlot(p)

###############################################################################
miniSeurat <- qs_read('MGCSeurat005.qs2')
p <- DimPlot(miniSeurat, group.by='orig.ident')
devPlot(p)
p <- FeaturePlot(miniSeurat, 'EYS')
devPlot(p)

################################################################################

miniSeurat <- subset(seuratObj, celltype %in% 
                         c('Muller glial cells', 'Retinal progenitor cells'))
miniSeurat <- removeRareFeatures(miniSeurat, 1)
miniSeurat <- removeRareFeatures(miniSeurat, 1, 'ATAC')
miniSeurat <- basicDimRed(miniSeurat)
miniSeurat <- jointUMAP(miniSeurat, FALSE, 
                        cutoff=0.005)

#qs_save(miniSeurat, 'mgrpSeurat.qs2')
miniSeurat <- qs_read('mgrpSeurat.qs2')
p <- DimPlot(miniSeurat, group.by='celltype') + ggtitle('Cell type annotation - RPC and MGC')
devPlot(p)

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
devPlot(p)


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
devPlot(p1)

################################################################################
miniSeurat <- qs_read('mgrpSeurat.qs2')
miniSeurat <- jointUMAP(miniSeurat, cutoff=0.02)
FeaturePlot(miniSeurat, 'RLBP1')
FeaturePlot(miniSeurat, 'TF')
FeaturePlot(miniSeurat, 'WIF1')
View(miniSeurat)
miniSeurat <- FindClusters(miniSeurat, resolution=0.7, graph.name = "wknn")
p <- DimPlot(miniSeurat, label=T)

devPlot(p)

p <- DotPlot(miniSeurat, c('CADPS', 'ZNF385B', 'FSTL5', 'EYS'))
devPlot(p)

v <- clusterMean(miniSeurat, c('RLBP1', 'TF', 'WIF1'), seq(0, 14))
View(v)

a <- FindMarkers(miniSeurat, ident.1 = c(5, 14, 8, 4), only.pos=T, min.diff.pct=0.2)
View(a)

microSeurat <- subset(miniSeurat, seurat_clusters %in% c(5, 14, 8, 4))
#microSeurat <- processSeurat(microSeurat, minFeatCells=1, cutoff=0.005)

microSeurat <- removeRareFeatures(microSeurat, 1)
microSeurat <- removeRareFeatures(microSeurat, 1, 'ATAC')
microSeurat <- basicDimRed(microSeurat)
microSeurat <- jointUMAP(microSeurat, FALSE, cutoff=0.02)
p <- DimPlot(microSeurat, group.by='orig.ident', label=TRUE) + labs(title='Newly identified MGC')
devPlot(p)
