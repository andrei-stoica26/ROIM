library(ggalluvial)
library(scLang)
library(hammers)

source('add_to_hammers.R')

seuratObj <- qs_read('annotatedSeurat.qs2')

p <- DimPlot(seuratObj, group.by='celltype', label=TRUE, repel=TRUE, 
             label.size=3) + NoLegend()
p <- centerTitle(p, 'Cell type annotation')
devPlot(p)

df <- repAnalysis(seuratObj, 'orig.ident', 'celltype')
p <- pvalRiverPlot(df, title='Overrepresented cell types')
devPlot(p)

################################################################################

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

################################################################################


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

df <- scColPairRatio(seuratObj, 'orig.ident', 'celltype')
df <- df[, c(1, 2, 6)]

df <- reshape2::melt(df)
df$significant <- case_when(
    df$value < 1 ~ "±",
    df$value > 1 & df$value < 1.5 ~ "+",
    df$value > 1.5 & df$value < 2 ~ "++",
    df$value > 2 ~ "+++",
)

p <- ggplot(df, aes(x = orig.ident, y = celltype, label = significant)) +
    geom_tile(aes(fill = value), color="white", size=0.5) +
    geom_text(color = "black") + 
    theme_classic() +
    scale_fill_viridis_c(option = "plasma", direction = 1, alpha = 1) +
    theme(axis.title = element_blank()) +
    labs(fill='Observed over\nexpected ratio')
p <- centerTitle(p, 'Observed over expected ratio')
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

gam <- computeGam(miniSeurat, sce, 'mgcGam')

gam <- qs_read('mgcGam.qs2')
res <- associationTest(gam)
topRes <- subset(res, pvalue < 0.05 & meanLogFC >= 2 & waldStat >= 100)

p <- pseudotimeHeatmapPlot(sce, rownames(topRes))


p <- featureWes(miniSeurat, 'HSPA6', idClass='orig.ident')
devPlot(p)





seuratObj$celltype <- factor(seuratObj$celltype, 
                                levels=c('Muller glial cells', 
                                         'Retinal progenitor cells',
                                         'Gliogenic progenitors',
                                         'Neurogenic progenitors',
                                         'Proliferative RPC',
                                         'Photoreceptor precursors',
                                         'Rods',
                                         'Bipolar cells',
                                         'Amacrine cells',
                                         'Cones',
                                         'ISL1+ retinal ganglion cells',
                                         'HTR2C+ retinal ganglion cells',
                                         'NEUROD6+ retinal ganglion cells',
                                         'ATF5+ horizontal cells',
                                         'ONECUT1+ horizontal cells',
                                         'RPE cells'))

celltypeCols <- c('Muller glial cells' = rgb(247/255,147/255,30/255), 
                  'Retinal progenitor cells' = rgb(210/255,0/255,0/255), 
                  'Gliogenic progenitors' = rgb(248/255,115/255,106/255), 
                  'Neurogenic progenitors' = rgb(180/255,0/255,150/255),
                  'Proliferative RPC' = rgb(200/255,50/255,15/255), 
                  'Photoreceptor precursors' = rgb(163/255,100/255,140/255), 
                  'Rods' = rgb(169/255,169/255,169/255), 
                  'Bipolar cells' = rgb(163/255,165/255,0/255), 
                  'Amacrine cells' = rgb(157/255,115/255,194/255), 
                  'Cones' = rgb(230/255,134/255,201/255), 
                  'ISL1+ retinal ganglion cells' = rgb(140/255,198/255,63/255),
                  'HTR2C+ retinal ganglion cells' = rgb(192/255,193/255,48/255), 
                  'NEUROD6+ retinal ganglion cells' = rgb(74/255,176/255,91/255), 
                  'ATF5+ horizontal cells' = rgb(97/255,156/255,255/255), 
                  'ONECUT1+ horizontal cells' = rgb(3/255,161/255,198/255), 
                  'RPE cells'= rgb(129/255,70/255,58/255))




rgb(183/255,76/255,171/255)
rgb(0/255,114/255,189/255)
p <- DimPlot(seuratObj, group.by='celltype', cols=celltypeCols) +
    theme(
        legend.key.size=unit(0.4, 'cm'),
        legend.text=element_text(size=8)
    )
p <- centerTitle(p, 'Cell type annotation')
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

genes <- c('NFIA', 'NFIX', 'SOX2', 'LHX2', 'HMGA1', 'ASCL1', 'SMARCA5', 'YAP1',
           'SOX9', 'STAT3', 'E2F3', 'FOXN4', 'MYB', 'FOXO3', 'SOX5', 'NFIB')

p <- pseudotimeHeatmapPlot(sce, genes)


p7 <- centerTitle(p7, "Expression of selected genes along pseudotime")
devPlot(p7)

FeaturePlot(miniSeurat, 'STAT3')
VlnPlot(miniSeurat, 'FTL', group.by='orig.ident')
FeaturePlot(miniSeurat, 'Lineage1')

################################################################################

psDF <- miniSeurat[['Lineage1']]
psDF <- psDF[order(psDF$Lineage1), drop=FALSE, ]
View(psDF)

gam <- qs_read('mgcGam.qs2')
res <- associationTest(gam)
topRes <- subset(res, pvalue < 0.05 & meanLogFC >= 2 & waldStat >= 100)
genes <- rownames(topRes)

mat <- scExpMat(miniSeurat, genes=genes)
mat <- mat[, rownames(psDF)]
df <- reshape2::melt(mat)
colnames(df) <- c('Gene', 'Cell', 'Expression')

breaks <- colnames(mat)[seq(1, dim(mat)[2], length.out=10)]
labels <- round(seq(0, max(psDF), length.out=10), 2)

p <- ggplot() + geom_tile(data=df, mapping=aes(x=Cell, y=Gene, fill=Expression)) +
    scale_fill_viridis() +
    scale_x_discrete(breaks=breaks, labels=labels) +
    labs(x='Pseudotime', y=NULL, fill='Expression level')

featureWes(miniSeurat, 'FTL', idClass='orig.ident')

featureWes(miniSeurat, 'Lineage1', idClass='orig.ident')

################################################################################

genes <- c('NFIA', 'NFIX', 'SOX2', 'LHX2', 'HMGA1', 'ASCL1', 'SMARCA5', 'YAP1',
           'SOX9', 'STAT3', 'E2F3', 'FOXN4', 'MYB', 'FOXO3', 'SOX5', 'NFIB')
mat <- scExpMat(miniSeurat, genes=genes)
#mat <- cluster_matrix(mat)
df <- reshape2::melt(mat)
View(df)
levels(df$Var1) <- rev(genes)
breaks <- colnames(mat)[seq(1, dim(mat)[2], length.out=10)]
labels <- round(seq(0, max(psDF), length.out=10), 2)

groupDF <- data.frame(Var1=genes, Group=c(rep('Rest', 4),
                                          rep('Reactivity', 6),
                                          rep('Proliferation', 3),
                                          rep('Restore rest', 3)))
levels(groupDF$Var1) <- rev(genes)

groupPlot <- ggplot() +
    geom_tile(data=groupDF,aes(x=1, y=Var1, fill=Group), width=1) +
    scale_fill_manual(name = "Group", values=wes_palette("GrandBudapest1")) +
    easy_remove_axes() + NoLegend() + theme(plot.margin=margin(0, 0, 0, 20))

heatPlot <- ggplot() + geom_tile(data=df, mapping=aes(x=Var2, y=Var1, fill=value)) +
    scale_fill_viridis() +
    scale_x_discrete(breaks=breaks, labels=labels) +
    labs(x='Pseudotime', y=NULL, fill='Expression level')

p7 <- (groupPlot + heatPlot + plot_layout(widths = c(0.02, 1)))

VlnPlot(miniSeurat, 'NFIB', group.by='orig.ident')

featureWes(miniSeurat, 'FTL', idClass='orig.ident')
