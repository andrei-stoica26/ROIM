addBasicQCInfo <- function(seuratObj){
    DefaultAssay(seuratObj) <- "RNA"
    message("Adding percentage of mitochondrial genes...")
    seuratObj[["percent.mt"]] <- PercentageFeatureSet(seuratObj, pattern = "^MT-")
    
    DefaultAssay(seuratObj) <- "ATAC"
    seuratObj <- NucleosomeSignal(seuratObj)
    seuratObj <- TSSEnrichment(seuratObj)
    seuratObj$blacklistFrac <- FractionCountsInRegion(seuratObj, assay='ATAC', 
                                                      regions=blacklist_hg38_unified)
    
    return(seuratObj)
}

cutoffVlnPlot <- function(seuratObj, col, cutoff, color='blue', linewidth=0.5)
    return(VlnPlot(seuratObj, col) + 
               geom_hline(yintercept=cutoff, color=color, linewidth=linewidth))

basicQC <- function(seuratObj, 
                    cutoffs = c(1000, 1200, 15, 0.01, 3, 2, 500, 1000)){
    nCells <- ncol(seuratObj)
    message('Number of cells before filtering: ', nCells, '.')
    seuratObj <- subset(
        seuratObj, 
        nCount_RNA >= cutoffs[1] &
            nCount_ATAC >= cutoffs[2] &
            percent.mt < cutoffs[3] &
            blacklistFrac < cutoffs[4] &
            TSS.enrichment > cutoffs[5] &
            nucleosome_signal < cutoffs[6] &        
            nFeature_RNA >= cutoffs[7] &
            nFeature_ATAC >= cutoffs[8]
    )
    message('Number of cells after filtering: ', ncol(seuratObj), '.')
    perc <- round((1 - ncol(seuratObj) / nCells) * 100, 2)
    message('Percentage of filtered cells: ', perc, '%.')
    return(seuratObj)
}

findRareFeatures <- function(seuratObj, nCells = 1){
    assayData <- LayerData(seuratObj, layer = "counts")
    df <- data.frame(Count = rowSums(assayData != 0))
    df <- subset(df, Count < nCells)
    return(df)
}

removeRareFeatures <- function(seuratObj, nCells = 10, assay = 'RNA'){
    DefaultAssay(seuratObj) <- assay
    rareFeatures <- findRareFeatures(seuratObj, nCells)
    freqFeatures <- setdiff(rownames(seuratObj), rownames(rareFeatures))
    message(paste0('Removing features found in less than ', nCells, ' cells...'))
    seuratObj[[assay]] <- subset(seuratObj[[assay]], features = freqFeatures)
    message(paste0(length(rownames(rareFeatures))), ' rare features removed.')
    return(seuratObj)
}