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
a <- subset(res, pvalue < 1e-04 & meanLogFC > 0.8 & waldStat > 80)
sigGenes <- rownames(a)

miniSeurat$pseudotimeBin <- cut(miniSeurat$Lineage1, breaks=4)
DoHeatmap(
    miniSeurat,
    features = intersect(sigGenes, VariableFeatures(miniSeurat)),
    group.by = "pseudotimeBin",
    label = FALSE
)

pt <- slingPseudotime(sce)[, 1] 
ord <- order(pt)
expr <- log1p(counts(sce)[sigGenes, ord])
exprScaled <- t(scale(t(expr)))

exprOrdered <- exprScaled[, ord] 


genePeak <- apply(exprOrdered, 1, function(gene) {
    sm <- zoo::rollmean(gene, 3, fill = NA)
    which.max(sm)
})
ordGenes <- order(genePeak)
exprOrdered <- exprOrdered[ordGenes, ]


ha <- HeatmapAnnotation(
    Pseudotime = pt[ord],
    col = list(Pseudotime = colorRamp2(c(min(pt), max(pt)), c("blue", "red")))
)


Heatmap(
    exprOrdered,
    top_annotation = ha,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    show_column_names = FALSE,
    show_row_names = FALSE,
    column_title = "Genes varying along Slingshot pseudotime"
)

