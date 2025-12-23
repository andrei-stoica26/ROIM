basicDimRed <- function(seuratObj){
    DefaultAssay(seuratObj) <- 'RNA'
    seuratObj <- NormalizeData(seuratObj)
    seuratObj <- FindVariableFeatures(seuratObj)
    seuratObj <- ScaleData(seuratObj)
    seuratObj <- RunPCA(seuratObj)
    
    DefaultAssay(seuratObj) <- 'ATAC'
    seuratObj <- FindTopFeatures(seuratObj, min.cutoff=5)
    seuratObj <- RunTFIDF(seuratObj)
    seuratObj <- RunSVD(seuratObj)
    
    DefaultAssay(seuratObj) <- 'RNA'
    return(seuratObj)
}

jointUMAP <- function(seuratObj, 
                      useHarmony = TRUE, 
                      cutoff = NULL,
                      dimsList = list(1:50, 2:40),
                      seed = 1){
    if (useHarmony){
        seuratObj <- with_seed(seed, 
                               RunHarmony(seuratObj,
                                          group.by.vars = 'orig.ident',
                                          reduction.use = 'pca',
                                          reduction.save = 'harmony_rna', 
                                          assay.use = 'RNA',
                                          project.dim = FALSE))
        seuratObj <- with_seed(seed,
                               RunHarmony(seuratObj,
                                          group.by.vars = 'orig.ident',
                                          reduction.use = 'lsi',
                                          reduction.save = 'harmony_atac',
                                          assay.use = 'ATAC',
                                          project.dim = FALSE))
        if (!is.null(cutoff)){
            umapDimsRNA <- chooseUMAPDims(seuratObj, reduction="harmony_rna",
                                          cutoff)
            umapDimsATAC <- chooseUMAPDims(seuratObj, reduction="harmony_atac",
                                           cutoff)
            message('Using ', umapDimsRNA, ' dimensions for RNA and ',
                    umapDimsATAC, ' dimensions for ATAC.')
            dimsList <- list(seq(umapDimsRNA),
                             seq(2, umapDimsATAC))
        }
        
        seuratObj <- FindMultiModalNeighbors(seuratObj,
                                             reduction.list=list('harmony_rna',
                                                                 'harmony_atac'),
                                             dims.list=dimsList,
                                             verbose=TRUE)
    } else{
        if (!is.null(cutoff)){
            umapDimsRNA <- chooseUMAPDims(seuratObj, reduction="pca", cutoff)
            umapDimsATAC <- chooseUMAPDims(seuratObj, reduction="lsi", cutoff)
            message('Using ', umapDimsRNA, ' dimensions for RNA and ',
                    umapDimsATAC, ' dimensions for ATAC.')
            dimsList <- list(seq(umapDimsRNA),
                             seq(2, umapDimsATAC))
            }
        seuratObj <- FindMultiModalNeighbors(seuratObj,
                                             reduction.list=list('pca', 'lsi'),
                                             dims.list=dimsList,
                                             verbose=TRUE)
        
    }
    
    seuratObj <- RunUMAP(seuratObj, nn.name = "weighted.nn", verbose=TRUE)
    DefaultAssay(seuratObj) <- 'RNA'
    return(seuratObj)
}