library(CelliD)

miniSeurat <- qs_read('MGCSeurat005.qs2')


runMCA <- function (X, nmcs = 50, features = NULL, reduction.name = "mca", 
          slot = "data", assay = DefaultAssay(X), ...) 
{
    InitAssay <- DefaultAssay(X)
    DefaultAssay(X) <- assay
    data_matrix <- as.matrix(GetAssayData(X, layer = slot))
    MCA <- RunMCA(X = data_matrix, nmcs = nmcs, features = features)
    geneEmb <- MCA$featuresCoordinates
    cellEmb <- MCA$cellsCoordinates
    stdev <- MCA$stdev
    X <- CelliD:::setDimMCSlot.Seurat(X = X, cellEmb = cellEmb, geneEmb = geneEmb, 
                      stdev = stdev, reduction.name = reduction.name)
    DefaultAssay(X) <- InitAssay
    return(X)
}

miniSeurat <- runMCA(miniSeurat)
qs_save(miniSeurat, 'MGCSeurat005MCA.qs2')

cgs <- GetCellGeneSet(miniSeurat, n.features=20)
df <- as.data.frame(table(unlist(cgs)))
df <- df[order(df$Freq, decreasing=TRUE), ]
View(df)

genes <- df[seq(90), 1]
extra <- grep('ENSG0', genes, value=TRUE)
genes <- setdiff(genes, extra)
length(genes)
p <- genesDimPlot(miniSeurat, genes, groupBy='orig.ident') + 
    ggtitle('Centers of mass - Top highly cell-specific genes')
devPlot(p)

featureWes(miniSeurat, 'TMEM88')

m <- scExpMat(miniSeurat, genes=genes)
m <- data.frame(cmdscale(dist(m)))
p <- densityPlot(m, 'MDS plot - Top highly cell-specific genes', drawNN=FALSE)
devPlot(p)
