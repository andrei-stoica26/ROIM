library(ggalluvial)
library(scLang)
library(hammers)

seuratObj <- qs_read('annotatedSeurat.qs2')

p <- DimPlot(seuratObj, group.by='celltype', label=TRUE, repel=TRUE, 
             label.size=3) + NoLegend()
p <- centerTitle(p, 'Cell type annotation')
devPlot(p)

df <- repAnalysis(seuratObj, 'orig.ident', 'celltype')
p <- pvalRiverPlot(df, title='Overrepresented cell types')
devPlot(p)

df <- scColPairPercs(seuratObj, 'orig.ident', 'celltype')
p <- ggplot(df, aes(x=orig.ident, 
                    y=perc, 
                    fill=celltype, 
                    stratum=celltype, 
                    alluvium=celltype)) +                 
    geom_col(width = 0.6,color="white")+
    geom_flow(width = 0.6, alpha=0.2, knot.pos = 0.1, color="white") +  
    theme_void()+ 
    theme(axis.text.x=element_text(size=10,vjust = 5))
p <- centerTitle(p, 'Changes in cell type representation')
devPlot(p)


df <- scColPairRatio(seuratObj, 'orig.ident', 'celltype')
df <- df[, c(1, 2, 6)]
p <- tilePlot(df, title='Observed over expected ratio', 
              doMelt=FALSE, showNumbers=FALSE, 
              palette='Plasma', 
              reverseColors=FALSE, 
              legendTitle='Observed over\nexpected ratio', 
              tileBoundaryWidth=0)
devPlot(p)

################################################################################

miniSeurat <- qs_read('miniSeurat.qs2')

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
