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
                              groupBy='celltype', label=TRUE, 
                              legend=FALSE, ...){
    mat <- extractLineageMat(slingshotObj, lineage)
    df <- points2Seg(mat)
    p <- DimPlot(seuratObj, group.by=groupBy, label=label, ...)  +
        geom_segment(data=df, aes(x=x, y=y, xend=xEnd, yend=yEnd),
                     arrow = arrow(length = unit(0.1, "cm")))
    if (!legend)
        p <- p + NoLegend()
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
    return(addSlingshotResults(seuratObj, 
                               slingshotObj, 
                               "Lineage", 
                               slingPseudotime))

addCurveweights <- function(seuratObj, slingshotObj)
    return(addSlingshotResults(seuratObj, 
                               slingshotObj, 
                               "Curveweight", 
                               slingCurveWeights))

computeGam <- function(seuratObj, sce, fileName){
    pseudotime <- slingPseudotime(sce)
    cellWeights <- slingCurveWeights(sce)
    counts <- LayerData(seuratObj, assay='RNA', layer='counts')
    x <- Sys.time()
    BPPARAM <- BiocParallel::bpparam()
    BPPARAM$workers <- 60
    
    gam <- fitGAM(counts,
                  pseudotime=pseudotime,
                  cellWeights=cellWeights,
                  parallel=TRUE,
                  BPPARAM=BPPARAM)
    fileName <- paste0(fileName, '.qs2')
    message('Saving file: ', fileName, '...')
    qs_save(gam, fileName)
    y <- Sys.time()
    print(y - x)
    return(gam)
}

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

createResultsTable <- function(seurats, res, fileName){
    sigRes <- subset(res, pvalue < 1e-2)
    sigRes$df <- c()
    sigRes <- addRanks(sigRes, signs = c(-1, 1, -1))
    genes <- rownames(sigRes)
    fractionCols <- paste0('Fraction', names(seurats))
    for (i in seq_along(names(seurats))){
        fractionCol <- fractionCols[i]
        condition <- names(seurats)[i]
        nCells <- ncol(seurats[[i]])
        fraction <- round(genePresence(seurats[[i]], genes)[genes, ]$nCells / 
                              nCells, 2)
        sigRes[fractionCol] <- fraction
    }
    sigRes$Fraction_Mean <- round(rowMeans(sigRes[, fractionCols]), 2)
    sigRes$Fraction_SD <- round(apply(sigRes[, fractionCols], 1, sd), 2)
    sigRes$Fraction_Maxdif <- apply(sigRes[, fractionCols], 1, 
                                    function(x) max(x) - min(x))
    sigRes <- sigRes[order(sigRes$Fraction_Maxdif, decreasing=TRUE), ]
    filename <- paste0(fileName, '.csv')
    message('Saving file: ', fileName, '.')
    write.csv(sigRes, fileName)
    return(sigRes)
}
