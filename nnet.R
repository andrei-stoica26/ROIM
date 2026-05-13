library(NeighbourNet)
library(ggplot2)
library(Seurat)
library(Signac)
library(qs2)


obj <- qs_read('miniSeurat.qs2')
genes <- select.gene(obj, min.cells=10)
obj <- obj |>
    prepare.seurat(genes = genes$genes) |>   # scale + PCA
    prepare.graph() |>                       # 30-NN graph
    select.cell() |>                         # subsample 
    prepare.reg(predictors = genes$tfs,      # local variance scaffolding
                responses  = genes$targets)

selGenes <- c('NFIA', 'NFIX', 'SOX2','LHX2', 'HMGA1', 'ASCL1', 'SMARCA5', 'YAP1',
           'SOX9', 'STAT3', 'E2F3', 'FOXN4', 'MYB', 'FOXO3', 'SOX5', 'NFIB')

obj <- run.nn.reg(obj, 
                  responses = selGenes, 
                  return.p.val = TRUE) |> build.meta.network() 
qs_save(obj, 'miniSeuratNnreg.qs2')

ctr.genes <- select.central.genes(obj)  
obj <- prepare.visualise(obj, central.genes = ctr.genes)

c(2, 3, 5, 6, 12, 13, 14, 16, 19)

cut <- mean(apply(obj@misc$NNet.mod$meta.network$p.val[,,1], 1, max))
visualise.network(obj, 1, meta.network = TRUE, cutoff = cut,
                  radius = c(.4,.7,.85,1), pie.radius = .04,
                  text.size = 3)
