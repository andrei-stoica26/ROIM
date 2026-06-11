library(Seurat)
library(Signac)
library(qs2)

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
                                      'Amacrine cells',
                                      'Cones',
                                      'Retinal ganglion cells',
                                      'HTR2C+ amacrine-like cells',
                                      'TBR1+ retinal ganglion-like cells',
                                      'ATF5+ horizontal-like cells',
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
                  'Amacrine cells' = rgb(157/255,115/255,194/255), 
                  'Cones' = rgb(230/255,134/255,201/255), 
                  'Retinal ganglion cells' = rgb(140/255,198/255,63/255),
                  'HTR2C+ amacrine-like cells' = rgb(192/255,193/255,48/255), 
                  'TBR1+ retinal ganglion-like cells' = rgb(74/255,176/255,91/255), 
                  'ATF5+ horizontal-like cells' = rgb(97/255,156/255,255/255), 
                  'Horizontal cells' = rgb(3/255,161/255,198/255), 
                  'RPE cells'= rgb(129/255,70/255,58/255))

p <- DimPlot(seuratObj, group.by='celltype', cols=celltypeCols) +
    theme(
        legend.key.size=unit(0.4, 'cm'),
        legend.text=element_text(size=8)
    )
p <- centerTitle(p, 'Cell type annotation')
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
    theme(axis.text.x=element_text(size=10,vjust = 5)) + scale_fill_discrete(palette=celltypeCols)
p <- centerTitle(p, 'Changes in cell type representation')

devPlot(p)
