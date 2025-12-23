extractLineageMat <- function(slingshotObj, lineage)
    return(slingCurves(slingshotObj)[[lineage]]$s)

points2Seg <- function(mat){
    df <- data.frame(x = mat[-nrow(mat), 1],
                     y = mat[-nrow(mat), 2],
                     xEnd = mat[-1, 1],
                     yEnd = mat[-1, 2])
    return(df)
}

singleLineagePlot <- function(seuratObj, slingshotObj, lineage,
                              groupBy='celltype', ...){
    mat <- extractLineageMat(slingshotObj, lineage)
    df <- points2Seg(mat)
    p <- DimPlot(seuratObj, group.by=groupBy, label=TRUE, ...) + NoLegend() +
        geom_segment(data=df, aes(x=x, y=y, xend=xEnd, yend=yEnd),
                     arrow = arrow(length = unit(0.1, "cm")))
    p <- henna::centerTitle(p, lineage)
    return(p)
}

addSlingshotResults <- function(seuratObj, slingshotObj, colStr, fun){
    nLineages <- ncol(fun(slingshotObj))
    for (i in 1:nLineages)
        seuratObj[[]][[paste0(colStr, i)]] <- tidyr::replace_na(fun(slingshotObj)[, i], NaN)
    return (seuratObj)
}

addLineages <- function(seuratObj, slingshotObj)
    return(addSlingshotResults(seuratObj, slingshotObj, "Lineage", slingPseudotime))

addCurveweights <- function(seuratObj, slingshotObj)
    return(addSlingshotResults(seuratObj, slingshotObj, "Curveweight", slingCurveWeights))