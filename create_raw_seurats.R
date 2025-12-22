annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))

folders <- paste0('sc-', c('0', '12', '24', 'ctr'))
idents <- c('0h', '12h', '24h', 'Control')
filePath <- paste0('/home/Andrei/ROInjuryData/Results/', folders, '/outs')
seurats <- mapply(function(folder, ident) 
    jointReadSeurats(folder, '../ROInjuryData/Results/', ident, annotation), 
    folders, 
    idents)

doublets <- mapply(predictDoublets, seurats, folders, SIMPLIFY=FALSE)
