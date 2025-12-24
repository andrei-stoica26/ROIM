#################First trajectory: Control - 0h - 12h - 24h#####################
miniSeurat <- qs_read('MGCSeurat005.qs2')

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

p1 <- featureWes(miniSeurat, 'Lineage1', idClass='orig.ident',
                 labelSize=4) + labs(color='Pseudotime') +
    geom_segment(data=df, aes(x=x, y=y, xend=xEnd, yend=yEnd),
                 arrow = arrow(length = unit(0.1, "cm")))

seurats <- SplitObject(miniSeurat, 'orig.ident')
names(seurats)
#gam <- computeGam(miniSeurat, sce, 'mgcGam')
gam <- qs_read('mgcGam.qs2')
res <- associationTest(gam)

w <- createResultsTable(seurats, 
                        res, 
                        'Genes associated with pseudotime - Ctrl-0h-12h-24h')

featureWes(miniSeurat, 'CRYAB', idClass = 'orig.ident')
DimPlot(miniSeurat, group.by='orig.ident')

#################Second trajectory: Control - 24h - 12h - 0h#####################
miniSeurat <- qs_read('MGCSeurat005.qs2')

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

p1 <- featureWes(miniSeurat, 'Lineage1', idClass='orig.ident',
                 labelSize=4) + labs(color='Pseudotime') +
    geom_segment(data=df, aes(x=x, y=y, xend=xEnd, yend=yEnd),
                 arrow = arrow(length = unit(0.1, "cm")))

#gam <- computeGam(miniSeurat, sce, 'mgcGam_ctrl-24-12-0')
seurats <- SplitObject(miniSeurat, 'orig.ident')
names(seurats)
gam <- qs_read('mgcGam_ctrl-24-12-0.qs2')
res <- associationTest(gam)

w <- createResultsTable(seurats, 
                        res, 
                        'Genes associated with pseudotime - Ctrl-24h-12h-0h')
