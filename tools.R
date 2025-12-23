clusterMean <- function(seuratObj, genes, clusters, doNormalize = T){
    message('Filtering expression matrix...')
    expression <- as.matrix(LayerData(seuratObj, layer='data')[genes, ])
    df <- data.table::transpose(data.frame(lapply(clusters, function(x){
        message('Assesing mean expression of selected genes in cluster ', x, '...')
        clusterExp <- expression[, which(seuratObj$seurat_clusters == x), drop = F]
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