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
topRes <- subset(res, pvalue < 0.05 & meanLogFC >= 2 & waldStat >= 150)

seurats <- SplitObject(miniSeurat, 'orig.ident')



w <- w[order(w$meanLogFC, decreasing=T), ]

p <- featureWes(miniSeurat, 'HSPA6', idClass='orig.ident')
devPlot(p)

p <- featureWes(miniSeurat, 'HSPH1', idClass='orig.ident')
devPlot(p)

p <- featureWes(miniSeurat, 'FTL', idClass='orig.ident')
devPlot(p)

w <- createResultsTable(seurats, 
                        res, 
                        'Genes associated with pseudotime.csv')
w <- w[order(w$waldStat, decreasing=T), ]
genes <- rownames(w)[seq(30)]
p <- featureWes(miniSeurat, genes[4], idClass='orig.ident')
devPlot(p)
invisible(lapply(genes, message))
m <- genesER(genes, 'human')
p <- newCnetplot(m)
devPlot(p)

w <- w[order(w$meanLogFC, decreasing=T), ]
genes <- rownames(w)[seq(30)]
p <- featureWes(miniSeurat, genes[2], idClass='orig.ident')
devPlot(p)
p <- featureWes(miniSeurat, 'PNLDC1', idClass='orig.ident')
devPlot(p)

################################################################################


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

gam <- qs_read('mgcGam.qs2')
res <- associationTest(gam)
topRes <- subset(res, pvalue < 0.05 & meanLogFC >= 2 & waldStat >= 150)
View(topRes)
genes <- rownames(topRes)
p <- contExpHeatmap(miniSeurat, genes)
p <- centerTitle(p, 'Genes strongly varying along pseudotime')
devPlot(p)

genes <- c('NFIA', 'NFIX', 'SOX2', 'LHX2', 'HMGA1', 'ASCL1', 'SMARCA5', 'YAP1',
           'SOX9', 'STAT3', 'E2F3', 'FOXN4', 'MYB', 'FOXO3', 'SOX5', 'NFIB')
p <- contExpHeatmap(miniSeurat, genes)
p <- centerTitle(p, 'Selected genes along pseudotime')
devPlot(p)

x <- 13
p <- featureWes(miniSeurat, 'EGFLAM', idClass='orig.ident')
devPlot(p)
p <- featureWes(miniSeurat, 'PEX5L', idClass='orig.ident')
devPlot(p)
invisible(lapply(genes, message))


m <- genesER(genes, 'human')
View(m@result)

topRes <- subset(res, pvalue < 0.05 & meanLogFC >= 2 & waldStat >= 100)
genes <- rownames(topRes)
m <- genesER(genes, 'human')
invisible(lapply(genes, message))



p <- newCnetplot(m)
devPlot(p)

activityMat <- GeneActivity(miniSeurat, assay='ATAC')
qs_save(activityMat, 'miniSeuratGeneActivity.qs2')
