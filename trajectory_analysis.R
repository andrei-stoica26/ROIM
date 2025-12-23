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

pseudotime <- slingPseudotime(sce)
cellWeights <- slingCurveWeights(sce)

counts <- LayerData(miniSeurat, assay='RNA', layer='counts')

x <- Sys.time()
BPPARAM <- BiocParallel::bpparam()
BPPARAM$workers <- 60

gam <- fitGAM(counts,
              pseudotime=pseudotime,
              cellWeights=cellWeights,
              parallel=TRUE,
              BPPARAM=BPPARAM)
qs_save(gam, 'mgcGam.qs2')

y <- Sys.time()
print(y - x)

gam <- qs_read('mgcGam.qs2')
res <- associationTest(gam)
View(res)
a <- subset(res, pvalue < 0.05 & meanLogFC > 0.5 & waldStat > 20)
View(a)

pseudotime <- slingPseudotime(sce)[, 1]
ord <- order(pseudotime, na.last = NA)
sigGenes <- rownames(res)[seq(20)]
expr <- log1p(counts(sce)[sigGenes, ord])
exprScaled <- t(scale(t(expr)))

pheatmap(
    exprScaled,
    cluster_cols = FALSE,
    cluster_rows = TRUE,
    show_colnames = FALSE,
    main = "Genes varying along Slingshot pseudotime"
)
