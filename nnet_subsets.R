library(NeighbourNet)
library(ggplot2)
library(Seurat)
library(Signac)
library(qs2)

prepareNet <- function(obj, selGenes){
    genes <- select.gene(obj, min.cells=10)
    obj <- obj |>
        prepare.seurat(genes = genes$genes) |>   # scale + PCA
        prepare.graph() |>                       # 30-NN graph
        prepare.reg(predictors = genes$tfs,      # local variance scaffolding
                    responses  = genes$targets)
    
    obj <- run.nn.reg(obj, 
                      responses = selGenes, 
                      return.p.val = TRUE) |> build.meta.network()
    
    ctr.genes <- select.central.genes(obj)  
    obj <- prepare.visualise(obj, central.genes = ctr.genes)
    return(obj)
}

createNetplots <- function(obj, metacells=seq(20)){
    cut <- mean(apply(obj@misc$NNet.mod$meta.network$p.val[,,1], 1, max))
    message('Creating network plots for Seurat object...')
    plots <- lapply(metacells, function(i){
        message('Plotting metacell ', i, '...')
        visualise.network(obj, i, meta.network = TRUE, cutoff = cut,
                          radius = c(.4,.7,.85,1), pie.radius = .04,
                          text.size = 2.5)
    })
    return(plots)
}

