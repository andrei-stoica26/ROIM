processSeurat <- function(seuratObj,
                          minFeatCells = 10,
                          assay = 'RNA',
                          varsToRegress = NULL,
                          cutoff = 0.1,
                          useHarmony=FALSE){
    seuratObj <- removeRareFeatures (seuratObj, minFeatCells, assay)
    seuratObj <- NormalizeData(seuratObj)
    seuratObj <- FindVariableFeatures(seuratObj)
    seuratObj <- ScaleData(seuratObj, vars.to.regress=varsToRegress)
    seuratObj <- RunPCA(seuratObj)
    if (useHarmony){
        seuratObj <- RunHarmony(seuratObj, group.by.vars = 'orig.ident',
                                assay.use = 'RNA')
        nUMAPDims <- chooseUMAPDims(seuratObj, reduction="harmony", cutoff)
        message(nUMAPDims, ' harmony dimensions will be used for UMAP.')
        seuratObj <- RunUMAP(seuratObj, dims=seq(nUMAPDims), 
                             reduction="harmony")
        
    } else {
        nUMAPDims <- chooseUMAPDims(seuratObj, 'pca', cutoff)
        message(nUMAPDims, ' PCA dimensions will be used for UMAP.')
        seuratObj <- RunUMAP(seuratObj, dims=seq(nUMAPDims))
    }
    
    return(seuratObj)
}

clusterMean <- function(seuratObj, genes, clusters, doNormalize = T){
    message('Filtering expression matrix...')
    expression <- as.matrix(LayerData(seuratObj, layer='data')[genes, ])
    df <- data.table::transpose(data.frame(lapply(clusters, function(x){
        message('Assesing mean expression of selected genes in cluster ', x, '...')
        clusterExp <- expression[, which(seuratObj$seurat_clusters == x), drop=FALSE]
        return(rowMeans(clusterExp))
    })))
    
    if(doNormalize){
        colList <- as.list(df)
        maxima <- apply(df, 2, max)
        df <- data.frame(mapply(function(x, y) x / y, colList, maxima))
    }
    
    rownames(df) <- clusters
    colnames(df) <- genes
    
    df$Mean <- rowMeans(df)
    df <- df[order(df$Mean, decreasing = T), ]
    df <- round(df, 2)
    return(df)
}


