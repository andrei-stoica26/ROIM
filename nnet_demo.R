library(NeighbourNet)
library(ggplot2)
library(Seurat)

load("luad.rda")

genes  <- select.gene(obj, min.cells = 10) # QC → TF / target lists

obj <- obj |>
    prepare.seurat(genes = genes$genes) |>   # scale + PCA
    prepare.graph() |>                       # 30-NN graph
    select.cell() |>                         # subsample 
    prepare.reg(predictors = genes$tfs,      # local variance scaffolding
                responses  = genes$targets)

top10 <- head(genes$targets, 10)           # demo: first 10 targets
obj   <- run.nn.reg(obj, responses = top10, return.p.val = TRUE) |>
    build.meta.network() 

ctr.genes <- select.central.genes(obj)  
obj <- prepare.visualise(obj, central.genes = ctr.genes)

visualise.network(obj, 2, 
                  radius = c(.4,.7,.85,1), pie.radius = .04,
                  text.size = 5)

cut <- mean(apply(obj@misc$NNet.mod$meta.network$p.val[,,1], 1, max))
visualise.network(obj, 3, meta.network = TRUE, cutoff = cut,
                  radius = c(.4,.7,.85,1), pie.radius = .04,
                  text.size = 5)