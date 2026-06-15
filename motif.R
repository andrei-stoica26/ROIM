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

# add motif information
DefaultAssay(miniSeurat) <- 'ATAC'
miniSeurat <- AddMotifs(miniSeurat, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pwm)
qs_save(miniSeurat, 'miniSeuratWithMotifs.qs2')

miniSeurat <- qs_read('miniSeuratWithMotifs.qs2')
da_peaks <- FindMarkers(
    object = miniSeurat,
    group.by = 'orig.ident',
    ident.1 = '0h',
    only.pos = TRUE,
    min.pct = 0.05,
)

top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005 & da_peaks$pct.1 > 0.2, ])
enriched.motifs <- FindMotifs(
    object = miniSeurat,
    
    features = top.da.peak
)

MotifPlot(
    object = miniSeurat,
    motifs = head(rownames(enriched.motifs))
)
