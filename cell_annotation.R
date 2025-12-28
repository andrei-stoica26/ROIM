seuratObj <- qs_read('preclusteringSeurat.qs2')
seuratObj <- FindClusters(seuratObj, resolution=3.6, graph.name='wknn')
DimPlot(seuratObj, label=TRUE, label.size=3) + NoLegend()

seuratObj <- addMetadataCategory(seuratObj, 'seurat_clusters', 'celltype',
                                 list(c(0, 1, 8, 28, 32, 16, 17, 31, 35, 10),
                                      c(2, 39, 5, 20, 26, 37, 38, 
                                        33, 14, 22, 9, 11, 24, 30, 41),
                                      3, 
                                      4, 
                                      c(6, 29, 18),
                                      7,
                                      c(12, 15, 25),
                                      c(13, 44),
                                      19,
                                      21,
                                      23, 
                                      27,
                                      34, 
                                      36,
                                      40, 
                                      42, 
                                      43),
                                 c('Cones',
                                   'Retinal progenitor cells',
                                   'Bipolar cells',
                                   'RPE cells',
                                   'Rods',
                                   'Muller glial cells',
                                   'Amacrine cells',
                                   'Photoreceptor precursors',
                                   'Gliogenic progenitors',
                                   'TMEM179+ retinal ganglion cells',
                                   'Horizontal cells',
                                   'Neurogenic progenitors',
                                   'HTR2C+ retinal ganglion cells',
                                   'ISL1+ retinal ganglion cells',
                                   'NEUROD6+ retinal ganglion cells',
                                   'Mesenchymal-like cells',
                                   'Myocyte-like cells'),
                                 'topMarkers',
                                 c('PDE6H/MPP4/RIMS2',
                                   'SFRP2/DAPL1/HES1',
                                   'CA10/CSMD3/VSX1',
                                   'TRPM1/CD96/TYR',
                                   'NRL/NR2E3/ROM1',
                                   'RLBP1/TF/WIF1',
                                   'PRDM13/TFAP2A/TFAP2B', 
                                   'PRDM1/VXN/LHX4',
                                   'TOP2A/CENPF/MKI67',
                                   'TMEM179/KCNIP1/KCNC2',
                                   'ONECUT1/ONECUT2/ONECUT3',
                                   'HES6/ATOH7/INSM1',
                                   'HTR2C/FYB2/EMX2',
                                   'ISL1/POU4F1/POU4F2',
                                   'NEUROD6/NEUROD2/POU3F3',
                                   'COL1A1/COL3A1/COL1A2',
                                   'TNNT2/DES/MYL1'))

qs_save(seuratObj, 'annotatedSeurat.qs2')

p <- DimPlot(seuratObj, group.by='celltype', label=TRUE, label.size=2.5, 
             repel=TRUE) + NoLegend() + ggtitle('Cell type annotation')
devPlot(p)

markers <- unique(unlist(str_split(unique(seuratObj$topMarkers), '/')))
markers <- c('TFAP2A', 'VSX1', 'PDE6H', 'TOP2A', 'ONECUT1', 'HTR2C', 'ISL1', 
             'COL1A2', 'TF', 'DES', 'NEUROD6', 'HES6','PRDM1', 'SFRP2', 
             'NRL','TRPM1', 'TMEM179')
p <- DotPlot(seuratObj, markers, group.by='celltype') + RotatedAxis() +
    theme(axis.title=element_text(size=10),
          axis.text=element_text(size=10),
          legend.title=element_text(size=10),
          legend.text=element_text(size=10))
devPlot(p)

a <- FindMarkers(seuratObj, group.by='celltype', ident.1='Mesenchymal-like cells',
                 only.pos=TRUE, logfc.threshold=5, min.pct=0.6)
a$pct.ratio <- a$pct.1 / a$pct.2

p <- DotPlot(seuratObj, rownames(a), group.by='celltype') + RotatedAxis() +
    theme(axis.title=element_text(size=10),
          axis.text=element_text(size=10),
          legend.title=element_text(size=10),
          legend.text=element_text(size=10))
devPlot(p)

a <- FindMarkers(seuratObj, group.by='celltype', ident.1='Myocyte-like cells',
                 only.pos=TRUE, logfc.threshold=5, min.pct=0.15)
a$pct.ratio <- a$pct.1 / a$pct.2
View(a)

p <- DotPlot(seuratObj, rownames(a), group.by='celltype') + RotatedAxis() +
    theme(axis.title=element_text(size=10),
          axis.text=element_text(size=10),
          legend.title=element_text(size=10),
          legend.text=element_text(size=10))
devPlot(p)

markers <- c('TOP2A', 'MKI67', 'CENPF', 'HES5', 'SOX9', 'HES6', 'ATOH7')
p <- DotPlot(seuratObj, markers, group.by='celltype') + RotatedAxis() +
    theme(axis.title=element_text(size=10),
          axis.text=element_text(size=10),
          legend.title=element_text(size=10),
          legend.text=element_text(size=10))
devPlot(p)



clusterMean(seuratObj, c('SLC1A3', 'TF', 'WIF1', 'CRYAB', 'RLBP1'), 
            c(7, 19, 24), F)
