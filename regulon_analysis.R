library(dorothea)
library(decoupleR)
library(tidyr)
library(tibble)
library(dplyr)
library(scLang)
library(patchwork)
library(Seurat)
library(qs2)
library(Signac)

miniSeurat <- qs_read('MGCSeurat005.qs2')

data(dorothea_hs, package = "dorothea")
colnames(dorothea_hs)
net <- dorothea_hs %>%
    dplyr::filter(confidence %in% c("A", "B", "C")) %>%
    dplyr::select(source = tf, target = target, mor = mor)

runMLM <- function(seuratObj, regulons){
    mat <- as.matrix(GetAssayData(seuratObj, assay = 'RNA',  layer = "data"))
    tf_activities <- run_mlm(mat, network = net, .source="source", .target="target")
    tfMat <- tf_activities %>%
        filter(statistic == "mlm") %>% 
        select(source, condition, score) %>%
        mutate(score = as.numeric(score)) 
    tfMat <- tfMat %>%
        pivot_wider(
            names_from = condition,
            values_from = score
        ) %>%
        column_to_rownames(var = "source")
    tfMat <- as(as.matrix(tfMat), 'dgCMatrix')
    dorothea_assay <- Seurat::CreateAssayObject(data = tfMat)
    Seurat::Key(dorothea_assay) <- "dorothea_"
    seuratObj[["dorothea"]] <- dorothea_assay
    return(seuratObj)
}

seurats <- SplitObject(miniSeurat, 'orig.ident')
regSeurats <- lapply(seq_along(seurats), function(i){
    message('Running MLM: ', names(seurats)[i], '...')
    seuratObj <- runMLM(seurats[[i]], net)
    DefaultAssay(seuratObj) <- 'dorothea'
    return(seuratObj)
})

names(regSeurats)
for (i in seq_along(regSeurats)){
    expMat <- scExpMat(regSeurats[[i]])
    df <- data.frame(Regulon = rownames(regSeurats[[i]]),
                     MeanActivity = rowMeans(expMat))
    write.csv(df, paste0('regulons', names(seurats)[i], '.csv'), 
                         row.names=FALSE)
}

rownames(seuratObj)

featureWes(seuratObj, feature = "HSF1", idClass='orig.ident')

violinPlot(miniSeurat, feature = "EBF1")

miniSeurat <- runMLM(miniSeurat, regulons)
qs_save(miniSeurat, 'MGCSeurat005MLM.qs2')

DefaultAssay(miniSeurat) <- 'dorothea'
a <- FindMarkers(seuratObj, group.by='orig.ident', ident.1='Control', only.pos=TRUE)
View(a)

reg <- 'HSF1'
regPlot <- function(seuratObj, reg1, reg2){
    p1 <- featureWes(seuratObj, feature=reg1, idClass='orig.ident', repel=TRUE)
    p2 <- violinPlot(seuratObj, feature=reg1)
    p3 <- featureWes(seuratObj, feature=reg2, idClass='orig.ident', repel=TRUE)
    p4 <- violinPlot(seuratObj, feature=reg2)
    p <- (p1 | p2) / (p3 | p4)
    return(p)
}

p <- regPlot(miniSeurat, 'HSF1', 'MAZ')
devPlot(p)
p <- regPlot(miniSeurat, 'MITF', 'NFE2L2')
devPlot(p)
p <- regPlot(miniSeurat, 'ETS2', 'SOX9')
devPlot(p)
p <- regPlot(miniSeurat, 'NR2F2', 'PRDM14')
devPlot(p)
 

