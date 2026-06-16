library(motifmatchr)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Seurat)
library(Signac)
library(qs2)

miniSeurat <- qs_read('miniSeurat.qs2')
pwm <- getMatrixSet(
    x = JASPAR2020,
    opts = list(species = 9606, all_versions = FALSE)
)
DefaultAssay(miniSeurat) <- 'ATAC'
miniSeurat <- AddMotifs(miniSeurat, genome = BSgenome.Hsapiens.UCSC.hg38, 
                        pfm = pwm)
qs_save(miniSeurat, 'miniSeuratWithMotifs.qs2')

################################################################################

findEnrichedMotifs <- function(seuratObj, id1){
    daPeaks <- FindMarkers(
        seuratObj,
        group.by = 'orig.ident',
        ident.1 = id1,
        only.pos = TRUE,
    )
    topDaPeaks <- rownames(daPeaks[daPeaks$p_val < 0.005 & 
                                       daPeaks$pct.1 > 0.1, ])
    enrichedMotifs <- FindMotifs(miniSeurat, topDaPeaks)
    return(enrichedMotifs)
}

################################################################################

miniSeurat <- qs_read('miniSeuratWithMotifs.qs2')
enrichedMotifs <- findEnrichedMotifs(miniSeurat, '24h')
MotifPlot(miniSeurat, head(rownames(enrichedMotifs)))
