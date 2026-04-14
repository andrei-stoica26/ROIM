contExpHeatmap <- function(seuratObj, genes,
                           colStr='Lineage1', 
                           xLab='Pseudotime',
                           fillLab='Expression level'){
    contDF <- metadataDF(seuratObj)[, colStr, drop=FALSE]
    contDF <- contDF[order(contDF[[colStr]]), drop=FALSE, ]
    
    mat <- scExpMat(seuratObj, genes=genes)
    mat <- mat[, rownames(contDF)]
    df <- reshape2::melt(mat)
    colnames(df) <- c('Gene', 'Cell', 'Expression')
    
    breaks <- colnames(mat)[seq(1, dim(mat)[2], length.out=10)]
    labels <- round(seq(min(contDF), max(contDF), length.out=10), 2)
    
    p <- ggplot() + geom_tile(data=df, mapping=aes(x=Cell, y=Gene, fill=Expression)) +
        scale_fill_viridis() +
        scale_x_discrete(breaks=breaks, labels=labels) +
        labs(x=xLab, y=NULL, fill=fillLab)
    
    return(p)
}

newCnetplot <- function(enrichmentResult, title = NULL, nCategories = 10){
    p <- enrichplot::cnetplot(enrichmentResult, showCategory = nCategories,
                              color_category = 'purple',
                              color_item = 'red',
                              color_edge = 'thistle') + NoLegend()
    p <- centerTitle(p, title)
    return(p)
}