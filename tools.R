processSeurat <- function(seuratObj,
                          assay = 'RNA',
                          minFeatCells = 10,
                          varsToRegress = NULL){
    seuratObj <- removeRareFeatures (seuratObj, minFeatCells, assay)
    seuratObj <- NormalizeData(seuratObj)
    seuratObj <- FindVariableFeatures(seuratObj)
    seuratObj <- ScaleData(seuratObj, vars.to.regress=varsToRegress)
    seuratObj <- RunPCA(seuratObj)
    nUMAPDims <- chooseUMAPDims(seuratObj)
    message(nUMAPDims, ' PCA dimensions will be used for UMAP.')
    seuratObj <- RunUMAP(seuratObj, dims=seq(nUMAPDims))
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