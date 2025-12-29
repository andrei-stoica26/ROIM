library(CelliD)

miniSeurat <- qs_read('MGCSeurat005.qs2')


runMCA <- function (X, nmcs = 50, features = NULL, reduction.name = "mca", 
          slot = "data", assay = DefaultAssay(X), ...) 
{
    InitAssay <- DefaultAssay(X)
    DefaultAssay(X) <- assay
    data_matrix <- as.matrix(GetAssayData(X, layer = slot))
    MCA <- RunMCA(X = data_matrix, nmcs = nmcs, features = features)
    geneEmb <- MCA$featuresCoordinates
    cellEmb <- MCA$cellsCoordinates
    stdev <- MCA$stdev
    X <- CelliD:::setDimMCSlot.Seurat(X = X, cellEmb = cellEmb, geneEmb = geneEmb, 
                      stdev = stdev, reduction.name = reduction.name)
    DefaultAssay(X) <- InitAssay
    return(X)
}

miniSeurat <- runMCA(miniSeurat)
qs_save(miniSeurat, 'MGCSeurat005MCA.qs2')
DimPlotMC(miniSeurat, 
          reduction = "mca", 
          group.by = "orig.ident", 
          features = c('AURKA', 'AURKB'), as.text = TRUE,
          size.feature = 2,
          size.feature.text = 3)

?DimPlotMC
