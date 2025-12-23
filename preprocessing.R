annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))

folders <- paste0('sc-', c('0', '12', '24', 'ctr'))
idents <- c('0h', '12h', '24h', 'Control')

filePath <- paste0('/home/Andrei/ROInjuryData/Results/', folders, '/outs')
seurats <- mapply(function(folder, ident) 
    jointReadSeurats(folder, '../ROInjuryData/Results/', ident, annotation), 
    folders, 
    idents)

seurats <- lapply(folders, function(folder) 
    qs_read(paste0(folder, 'RawSeurat.qs2')))

seurats <- lapply(seurats, addBasicQCInfo)
invisible(mapply(function(x, y) 
    qs_save(x, paste0(y, 'QCSeurat.qs2')), seurats, folders))

seurats <- lapply(folders, function(folder) 
    qs_read(paste0(folder, 'QCSeurat.qs2')))
seurats <- lapply(seurats, basicQC)

#doublets <- mapply(predictDoublets, seurats, folders, SIMPLIFY=FALSE)
doublets <- lapply(folders, function(folder) 
    qs_read(paste0(folder, 'DoubletsRNA.qs2')))
seurats <- mapply(addDoublets, seurats, doublets, SIMPLIFY=FALSE)
seurats <- lapply(seurats, function(seuratObj) removeDoublets(seuratObj, 
                                                              unitClass))
seuratObj <- mergeSeurats(seurats, idents)
qs_save(seuratObj, 'mergedSeurat.qs2')

seuratObj <- qs_read('mergedSeurat.qs2')
seuratObj <- removeRareFeatures(seuratObj, 10, 'RNA')
seuratObj <- removeRareFeatures(seuratObj, 10, 'ATAC')
seuratObj <- basicDimRed(seuratObj)
seuratObj <- jointUMAP(seuratObj)
qs_save(seuratObj,'preclusteringSeurat.qs2')

seuratObj <- qs_read('preclusteringSeurat.qs2')
seuratObj <- FindClusters(seuratObj, resolution=3, graph.name='wknn')
DimPlot(seuratObj, label=TRUE) + NoLegend()

v <- clusterMean(seuratObj, 
                 c('SLC1A3', 'RLBP1','SLITRK2'),
                 c(7, 13, 23, 9))

miniSeurat <- subset(seuratObj, seurat_clusters==7)
miniSeurat <- basicDimRed(miniSeurat)
miniSeurat <- jointUMAP(miniSeurat, FALSE, 0.1)
DimPlot(miniSeurat)







