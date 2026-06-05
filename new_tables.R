library(hammers)
library(Seurat)
library(Signac)
library(qs2)

createFractionTable <- function(seuratObj, celltype){
    seurats <- SplitObject(seuratObj, "orig.ident")
    message('Generating results...')
    df <- Reduce(cbind, lapply(seurats, function(seurat){
        v <- genePresence(seurat)
        v$Fraction <- round(v$nCells / dim(seurat)[2], 2)
        v <- v[sort(rownames(v)), ]
        return(v[, 3, drop=FALSE])
    }))
    
    colnames(df) <- paste0(celltype, '_fraction_', names(seurats))
    df[[paste0(celltype, '_Mean')]] <- round(rowMeans(df), 2)
    df[[paste0(celltype, '_SD')]] <- round(apply(df[, seq(4)], 1, sd), 2)
    df[[paste0(celltype, '_maxdif')]] <- apply(df[, seq(4)], 1, function(x) max(x) - min(x))
    df <- df[order(df[[paste0(celltype, '_maxdif')]], decreasing=TRUE), ]
    return(df)
}

write_summaries <- function(seuratObj){
    types <- unique(seuratObj$celltype)
    shortTypes <- c('BC', 'RPC', 'MGC', 'AC', 'Cones', 'Rods', 'ONECUT1+_HC',
                    'NEUROD6+_RGC', 'Gliogenic_prog', 'Photoreceptor_prec',
                    'Prol_RPC', 'RPE', 'ATF5+_HC', 'Neurogenic_prog', 'ISL1+_RGC',
                    'HTR2C+_RGC')
    tables <- mapply(function(type, shortType){
        subSeurat <- subset(seuratObj, celltype==type)
        df <- createFractionTable(subSeurat, shortType)
        write.csv(df, paste0('exp_', type, '.csv'))
        return(df)
    }, types, shortTypes, SIMPLIFY=FALSE)
    
    res <- Reduce(cbind, tables)
    write.csv(res, 'exp_summary.csv')
    return(res)
}

seuratObj <- qs_read('annotatedSeurat.qs2')
res <- write_summaries(seuratObj)

################################################################################
miniSeurat <- qs_read('miniSeurat.qs2')
for (id in c('Control', '0h', '12h', '24h')){
    m <- FindMarkers(miniSeurat, group.by='orig.ident', ident.1=id, only.pos=TRUE)
    write.csv(m, paste0('markers_', id, '.csv'))
}


