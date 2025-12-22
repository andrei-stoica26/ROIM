predictDoublets <- function(seuratObj, dataset, nRuns = 100, dbr = 0.02, BPPARAM = bpparam(), 
                            assay = 'RNA', ...)
{
    DefaultAssay(seuratObj) <- assay
    message("Converting Seurat object to SingleCellExperiment...")
    sce_seurat <- as.SingleCellExperiment(seuratObj)
    message("Predicting doublets...")
    doublets <- bplapply(seq_len(nRuns), function(i) {
        message(paste0("Predicting doublets: run ", i, "..."))
        res <- scDblFinder(sce_seurat, dbr = dbr, ...)$scDblFinder.class
        message(paste0("Run ", i, " completed."))
        return(res)
    }, BPPARAM = BPPARAM)
    doublets <- data.frame(doublets)
    colnames(doublets) <- paste0("unit.class", seq_len(nRuns))
    fileName <- paste0(dataset, 'Doublets', assay, '.qs2')
    message(paste0("Saving file: ", fileName, "..."))
    qs_save(doublets, fileName)
    return(doublets)
}

doubletListings <- function(df)
    return(apply(df, 1, function(x) sum (x == 'doublet')))

doubletCounts <- function(listings, nRuns = 100)
    return(sapply(0:nRuns, function(x) sum(listings == x)))

addDoublets <- function(seuratObj, doubletsDF, colStr = 'unitClass'){
    listings <- doubletListings(doubletsDF)
    nRuns <- ncol(doubletsDF)
    dblMean <- sum(listings) / nRuns
    message(paste0("The average number of doublets predicted per run was ",  dblMean))
    dblCounts <- doubletCounts(listings, nRuns)
    revCumsum <- revcumsum(dblCounts)
    
    dblCutoff <- which.min(abs(revCumsum - dblMean))
    message(paste0("Cell units predicted to be doublets in at least ",  dblCutoff, " runs will be regarded as doublets."))
    message(paste0("Taking this cutoff implies ", revCumsum[dblCutoff], " predicted doublets."))
    
    unitClass <- rep('singlet', ncol(seuratObj))
    unitClass[which(listings >= dblCutoff)] <- 'doublet'
    seuratObj[[colStr]] <- factor(x=unitClass, levels=c("singlet", "doublet"))
    return(seuratObj)
}