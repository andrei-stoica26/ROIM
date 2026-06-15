seuratObj <- qs_read('preclusteringSeurat.qs2')
seuratObj <- FindClusters(seuratObj, resolution=4, graph.name='wknn')
DimPlot(seuratObj, label=T) + NoLegend()

seuratObj <- addMetadataCategory(seuratObj, 'seurat_clusters', 'celltype',
                                 list(0,
                                      1,
                                      c(2, 43, 17, 21, 14, 40, 24),
                                      c(3, 4, 8, 26, 37, 28, 12, 36, 32, 18, 22, 41, 35),
                                      c(5, 33, 9, 15, 10, 46, 25, 27, 30),
                                      c(6, 16, 29),
                                      c(7, 11),
                                      c(13, 20),
                                      c(19, 45),
                                      23,
                                      31,
                                      34,
                                      38,
                                      39,
                                      42,
                                      44),
                                 c('Rods',
                                   'RPE cells',
                                   'Retinal progenitor cells',
                                   'Cones',
                                   'Muller glial cells',
                                   'Amacrine cells',
                                   'Photoreceptor precursors',
                                   'Bipolar cells',
                                   'Gliogenic progenitors',
                                   'Proliferative RPC',
                                   'Horizontal cells',
                                   'Neurogenic progenitors',
                                   'HTR2C+ amacrine-like cells',
                                   'Retinal ganglion cells',
                                   'ATF5+ amacrine-like cells',
                                   'TBR1+ retinal ganglion-like cells'),
                                 'topMarkers',
                                 c('NRL/NR2E3/ROM1',
                                   'TRPM1/CD96/TYR',
                                   'SFRP2/DAPL1/HES1',
                                   'PDE6H/MPP4/RIMS2',
                                   'RLBP1/TF/CD44',
                                   'PRDM13/TFAP2A/TFAP2B',
                                   'PRDM1/VXN/LHX4',
                                   'CA10/CSMD3/VSX1',
                                   'MOXD1/FLT1/MYO3A',
                                   'TOP2A/CENPF/MKI67',
                                   'ONECUT1/ONECUT2/ONECUT3',
                                   'HES6/ATOH7/INSM1',
                                   'HTR2C/GRM8/FYB2',
                                   'ISL1/POU4F1/POU4F2',
                                   'ATF5/UNC5B/DLGAP1',
                                   'TBR1/NTRK3/LRRC7'))

qs_save(seuratObj, 'annotatedSeuratNew.qs2')

DimPlot(seuratObj, group.by='celltype', label=T, repel=T, label.size=3) + NoLegend()

a <- FindMarkers(seuratObj, group.by='celltype', ident.1='Photoreceptor precursors',
                 ident.2 = 'ATF5+ horizontal-like cells', only.pos=TRUE, logfc.threshold=1)

View(a)
b <- FindMarkers(seuratObj, group.by='celltype', ident.2='Photoreceptor precursors',
                 ident.1 = 'ATF5+ horizontal-like cells', only.pos=TRUE, logfc.threshold=1)

View(b)


v <- clusterMean(seuratObj, 
                 c('SFRP2', 'DAPL1', 'HES1'),
                 c(43, 2, 17, 21, 24, 14, 42, 45, 19, 30, 5, 10, 15, 33, 9, 25, 23, 46, 27, 34))
View(v)

w <- clusterMean(seuratObj, 
                 c('SOX2', 'SLC1A3', 'CD44', 'SLITRK2', 'COL4A2'),
                 c(43, 2, 17, 21, 24, 14, 42, 45, 19, 30, 5, 10, 15, 33, 9, 25, 23, 46, 27, 34))
View(w)

w <- clusterMean(seuratObj, 
                 c('RLBP1', 'TF', 'CD44', 'SLITRK2'),
                 c(43, 2, 17, 21, 24, 14, 42, 45, 19, 30, 5, 10, 15, 33, 9, 25, 23, 46, 27, 34))
View(w)

v

a <- FindMarkers(seuratObj, ident.1=42, logfc.threshold=1, only.pos=TRUE)
View(a)

FeaturePlot(seuratObj, 'GBE1')

DimPlot(seuratObj, label=T) + NoLegend()