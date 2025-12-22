jointReadSeurats <- function(folder, folderPrefix, ident, annotation){
    folderPath <- file.path(folderPrefix, folder, 'outs')
    counts <- Read10X(file.path(folderPath, 'filtered_feature_bc_matrix'))
    fragpath <- file.path(folderPath, 'atac_fragments.tsv.gz')
    seuratObj <- CreateSeuratObject(
        counts = counts$`Gene Expression`,
        project = ident,
        assay = "RNA"
    )
    seuratObj[["ATAC"]] <- CreateChromatinAssay(
        counts = counts$Peaks,
        sep = c(":", "-"),
        fragments = fragpath,
        annotation = annotation
    )
    qs_save(seuratObj, paste0(folder, 'RawSeurat.qs2'))
    return(seuratObj)
}

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

jointIntegration <- function(seuratObj){
    seuratObj <- RunHarmony(seuratObj, group.by.vars = 'orig.ident', 
                            reduction.use = 'pca', 
                            reduction.save = 'harmony_rna', assay.use = 'RNA',
                            project.dim = FALSE)
    seuratObj <- RunHarmony(seuratObj, group.by.vars = 'orig.ident', 
                            reduction.use = 'lsi', 
                            reduction.save = 'harmony_atac', assay.use = 'ATAC',
                            project.dim = FALSE)
    seuratObj <- FindMultiModalNeighbors(seuratObj, 
                                         reduction.list = list('harmony_rna', 
                                                               'harmony_atac'), 
                                         dims.list = list(1:50, 2:40),
                                         verbose = TRUE)
    seuratObj <- RunUMAP(seuratObj, nn.name = "weighted.nn", verbose=TRUE)
    DefaultAssay(seuratObj) <- 'RNA'
    return(seuratObj)
}
