chooseUMAPDims <- function(seuratObj, reduction){
    pct <- seuratObj[[reduction]]@stdev / sum(seuratObj[[reduction]]@stdev) * 100
    nUMAPDims <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), 
                      decreasing = TRUE)[1] + 1
    return(nUMAPDims)
}