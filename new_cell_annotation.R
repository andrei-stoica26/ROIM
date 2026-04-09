seuratObj <- qs_read('preclusteringSeuratFull.qs2')
seuratObj <- FindClusters(seuratObj, resolution=3.6, graph.name='wknn')
DimPlot(seuratObj, label=T) + NoLegend()
seuratObj <- addMetadataCategory(seuratObj, 'seurat_clusters', 'celltype',
                                 list(19, 41, 42
                                 ),
                                 c('Mixed cells',
                                   'Myocyte-like cells',
                                   'Mesenchymal-like cells'),
                                 'topMarkers',
                                 c('KCNH7/TMEM179/SGCZ',
                                   'TNNT2/DES/MYL1',
                                   'COL1A1/COL3A1/COL1A2'))

seuratObj <- subset(seuratObj, !seurat_clusters %in% c(19, 41, 42))
seuratObj <- removeRareFeatures(seuratObj, 10, 'RNA')
seuratObj <- removeRareFeatures(seuratObj, 10, 'ATAC')
seuratObj <- basicDimRed(seuratObj)
seuratObj <- jointUMAP(seuratObj)
qs_save(seuratObj,'preclusteringSeurat.qs2')

################################################################################

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
                                   'ONECUT1+ horizontal cells',
                                   'Neurogenic progenitors',
                                   'HTR2C+ retinal ganglion cells',
                                   'ISL1+ retinal ganglion cells',
                                   'ATF5+ horizontal cells',
                                   'NEUROD6+ retinal ganglion cells'),
                                 'topMarkers',
                                 c('NRL/NR2E3/ROM1',
                                   'TRPM1/CD96/TYR',
                                   'SFRP2/DAPL1/HES1',
                                   'PDE6H/MPP4/RIMS2',
                                   'RLBP1/TF/WIF1',
                                   'PRDM13/TFAP2A/TFAP2B',
                                   'PRDM1/VXN/LHX4',
                                   'CA10/CSMD3/VSX1',
                                   'MOXD1/FLT1/MYO3A',
                                   'TOP2A/CENPF/MKI67',
                                   'ONECUT1/ONECUT2/ONECUT3',
                                   'HES6/ATOH7/INSM1',
                                   'HTR2C/FYB2/EMX2',
                                   'ISL1/POU4F1/POU4F2',
                                   'ATF5/UNC5B/GBE1',
                                   'NEUROD6/NEUROD2/POU3F3'))

qs_save(seuratObj, 'annotatedSeurat.qs2')

DimPlot(seuratObj, label=T, repel=T, label.size=3, group.by='celltype') + NoLegend()

miniSeurat <- subset(seuratObj, celltype=='Muller glial cells')
miniSeurat <- removeRareFeatures(miniSeurat, 1)
miniSeurat <- removeRareFeatures(miniSeurat, 1, 'ATAC')
miniSeurat <- basicDimRed(miniSeurat)
miniSeurat <- jointUMAP(miniSeurat, FALSE, cutoff=0.05)
qs_save(miniSeurat, 'miniSeurat.qs2')

DimPlot(miniSeurat, group.by='orig.ident', label=TRUE)







