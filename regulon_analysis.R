library(dorothea)
library(decoupleR)

miniSeurat <- qs_read('MGCSeurat005.qs2')

data(dorothea_hs, package = "dorothea")
colnames(dorothea_hs)
net <- dorothea_hs %>%
    dplyr::filter(confidence %in% c("A", "B", "C")) %>%
    dplyr::select(source = tf, target = target, mor = mor)

runMLM <- function(seuratObj, regulons){
    mat <- as.matrix(Seurat::GetAssayData(seuratObj, assay = 'RNA', 
                                          layer = "data"))
    tf_activities <- run_mlm(mat, network = net, .source="source", .target="target")
    dorothea_assay <- Seurat::CreateAssayObject(data = tf_activities)
    Seurat::Key(dorothea_assay) <- "dorothea_"
    seuratObj[["dorothea"]] <- dorothea_assay
    return(seuratObj)
}

miniSeurat <- runMLM(miniSeurat, regulons)
qs_save(miniSeurat, 'MGCSeurat005MLM.qs2')
