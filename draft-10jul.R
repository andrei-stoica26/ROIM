library(Seurat)
library(Signac)
library(qs2)
library(hammers)
library(henna)
library(scLang)
library(ggplot2)
library(ggalluvial)
library(CSOA)

source('nnet_subsets.R')

seuratObj <- qs_read('annotatedSeuratNew.qs2')

seuratObj$celltype <- factor(seuratObj$celltype, 
                             levels=c('Muller glial cells', 
                                      'Retinal progenitor cells',
                                      'Gliogenic progenitors',
                                      'Neurogenic progenitors',
                                      'Proliferative RPC',
                                      'Photoreceptor precursors',
                                      'Rods',
                                      'Bipolar cells',
                                      'NRP1+ amacrine cells',
                                      'HTR2C+ amacrine cells',
                                      'ATF5+ amacrine cells',
                                      'Cones',
                                      'ISL1+ retinal ganglion cells',
                                      'TBR1+ retinal ganglion cells',
                                      'Horizontal cells',
                                      'RPE cells'))

celltypeCols <- c('Muller glial cells' = rgb(247/255,147/255,30/255), 
                  'Retinal progenitor cells' = rgb(210/255,0/255,0/255), 
                  'Gliogenic progenitors' = rgb(248/255,115/255,106/255), 
                  'Neurogenic progenitors' = rgb(180/255,0/255,150/255),
                  'Proliferative RPC' = rgb(200/255,50/255,15/255), 
                  'Photoreceptor precursors' = rgb(163/255,100/255,140/255), 
                  'Rods' = rgb(169/255,169/255,169/255), 
                  'Bipolar cells' = rgb(163/255,165/255,0/255), 
                  'NRP1+ amacrine cells' = rgb(137/255,95/255,174/255), 
                  'HTR2C+ amacrine cells' = rgb(170/255,120/255,196/255),
                  'ATF5+ amacrine cells' = rgb(200/255,125/255,198/255),
                  'Cones' = rgb(230/255,134/255,201/255), 
                  'ISL1+ retinal ganglion cells' = rgb(140/255,198/255,63/255),
                  'TBR1+ retinal ganglion cells' = rgb(74/255,176/255,91/255), 
                  'Horizontal cells' = rgb(3/255,161/255,198/255), 
                  'RPE cells'= rgb(129/255,70/255,58/255))

p <- DimPlot(seuratObj, group.by='celltype', cols=celltypeCols) +
    theme(
        legend.key.size=unit(0.4, 'cm'),
        legend.text=element_text(size=8)
    )
p <- centerTitle(p, 'Cell type annotation')
devPlot(p)

p <- featureWes(seuratObj, 'TBR1', idClass='celltype', repel=TRUE)
devPlot(p)

p <- featureWes(seuratObj, 'DLGAP1', idClass='celltype', repel=TRUE)
devPlot(p)

p <- featureWes(seuratObj, 'ATF5', idClass='celltype', repel=TRUE)
devPlot(p)

p <- featureWes(seuratObj, 'HTR2C', idClass='celltype', repel=TRUE)
devPlot(p)
################################################################################


df <- scColPairPercs(seuratObj, 'orig.ident', 'celltype')

df$orig.ident <- factor(df$orig.ident, 
                        levels=c('Control', '0h', '12h', '24h'))
p <- ggplot(df, aes(x=orig.ident, 
                    y=perc, 
                    fill=celltype, 
                    stratum=celltype, 
                    alluvium=celltype)) +                 
    geom_col(width = 0.6,color="white")+
    geom_flow(width = 0.6, alpha=0.2, knot.pos = 0.1, color="white") +  
    theme_void()+ 
    theme(axis.text.x=element_text(size=10,vjust = 5)) + scale_fill_discrete(palette=celltypeCols)
p <- centerTitle(p, 'Changes in cell type representation')

devPlot(p)

################################################################################

set.seed(123)

miniSeurat <- qs_read('miniSeurat.qs2')
seurats <- SplitObject(miniSeurat, 'orig.ident')

selGenes <- c('NFIA', 'NFIX', 'SOX2','LHX2', 'HMGA1', 'ASCL1', 'SMARCA5', 'YAP1',
              'SOX9', 'STAT3', 'E2F3', 'FOXN4', 'MYB', 'FOXO3', 'SOX5', 'NFIB')

obj <- seurats[['0h']]
net0h <- prepareNet(obj, selGenes)

obj <- seurats[['12h']]
net12h <- prepareNet(obj, selGenes)

obj <- seurats[['24h']]
net24h <- prepareNet(obj, selGenes)

obj <- seurats[['Control']]
netCtr <- prepareNet(obj, selGenes)

netSeurats <- list(net0h, net12h, net24h, netCtr)
names(netSeurats) <- names(seurats)

rm(net0h)
rm(net12h)
rm(net24h)
rm(netCtr)

v[[6]]
v <- createNetplots(netSeurats[[4]], c(seq(2, 5), 7, 9))
netPlots <- mapply(function(x, y) createNetplots(x, y),
                   netSeurats, 
                   list(seq(20),
                        seq(19),
                        seq(20),
                        c(seq(2, 5), 7, 9)),
                   SIMPLIFY=FALSE)
