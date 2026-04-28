buildCSOANetworks <- function(seuratObj, col='orig.ident'){
    seurats <- SplitObject(miniSeurat, col)
    seurats <- lapply(seurats, function(x){
        x <- removeRareGenes(x, 10)
        x <- FindVariableFeatures(x)
        return(x)
    })
    mats <- lapply(seurats, function(seuratObj) scExpMat(seuratObj))
    res <- mapply(function(x, y){
        overlapDF <- generateOverlaps(x, pairs=getPairs(VariableFeatures(y)))
        overlapDF <- CSOA:::prefilterOverlaps(overlapDF)
        return(overlapDF)
    },
    mats,
    seurats,
    SIMPLIFY=FALSE)
    return(res)
}