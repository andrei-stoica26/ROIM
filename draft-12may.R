library(Seurat)
library(Signac)
library(qs2)
library(slingshot)
library(ggplot2)
library(scLang)
library(viridis)
library(tradeSeq)
library(NeighbourNet)
library(henna)
library(hammers)
library(wesanderson)
library(ggeasy)
library(patchwork)

source('trajectory_analysis_tools.R')
source('visualization.R')

miniSeurat <- qs_read('miniSeurat.qs2')


sce <- as.SingleCellExperiment(miniSeurat)
sce <- slingshot(sce,
                 clusterLabels = 'orig.ident',
                 reducedDim = 'UMAP',
                 start.clus = 'Control',
                 dist.method = 'simple')

miniSeurat <- addLineages(miniSeurat, sce)
miniSeurat[['Pseudotime']] <- miniSeurat[['Lineage1']]

gam <- qs_read('mgcGam.qs2')
res <- associationTest(gam)
res <- res[order(res$waldStat, decreasing=TRUE),]
topGenes <- rownames(subset(res, waldStat > 1000))
View(res)

obj <- miniSeurat
genes <- select.gene(obj, min.cells=10)
obj <- obj |>
    prepare.seurat(genes = genes$genes) |>   # scale + PCA
    prepare.graph() |>                       # 30-NN graph
    select.cell() |>                         # subsample 
    prepare.reg(predictors = genes$tfs,      # local variance scaffolding
                responses  = genes$targets)

obj <- run.nn.reg(obj, 
                  responses = topGenes, 
                  return.p.val = TRUE) |> build.meta.network() 

ctr.genes <- select.central.genes(obj)  
obj <- prepare.visualise(obj, central.genes = ctr.genes)

cut <- mean(apply(obj@misc$NNet.mod$meta.network$p.val[,,1], 1, max))

x <- 20
p <- visualise.network(obj, x, meta.network = TRUE, cutoff = cut,
                  radius = c(.4,.7,.85,1), pie.radius = .04,
                  text.size = 3) 
p <- centerTitle(p, paste0('Regulatory network - Metacell ', x))
devPlot(p)

################################################################################
selGenes <- c('NFIA', 'NFIX', 'SOX2','LHX2', 'HMGA1', 'ASCL1', 'SMARCA5', 
              'YAP1','SOX9', 'STAT3', 'E2F3', 'FOXN4', 'MYB', 'FOXO3', 
              'SOX5', 'NFIB')

groupDF <- data.frame(Var1=selGenes, Group=c(rep('Rest', 4),
                                          rep('Reactivity', 6),
                                          rep('Proliferation', 3),
                                          rep('Restore rest', 3)))
levels(groupDF$Var1) <- rev(selGenes)

groupPlot <- ggplot() +
    geom_tile(data=groupDF,aes(x=1, y=Var1, fill=Group), width=1) +
    scale_fill_manual(name = "Group", values=wes_palette("GrandBudapest1")) +
    easy_remove_axes() + NoLegend() + theme(plot.margin=margin(0, 0, 0, 20))

heatPlot <- contExpHeatmap (miniSeurat, selGenes)

p <- (groupPlot + heatPlot + plot_layout(widths = c(0.02, 1)))
p <- centerTitle(p, 'Expression of selected genes along pseudotime')
devPlot(p)

obj <- miniSeurat
genes <- select.gene(obj, min.cells=10)
obj <- obj |>
    prepare.seurat(genes = genes$genes) |>   # scale + PCA
    prepare.graph() |>                       # 30-NN graph
    select.cell() |>                         # subsample 
    prepare.reg(predictors = genes$tfs,      # local variance scaffolding
                responses  = genes$targets)

obj <- run.nn.reg(obj, 
                  responses = selGenes, 
                  return.p.val = TRUE) |> build.meta.network() 

ctr.genes <- select.central.genes(obj)  
obj <- prepare.visualise(obj, central.genes = ctr.genes)

cut <- mean(apply(obj@misc$NNet.mod$meta.network$p.val[,,1], 1, max))

x <- 2
p <- visualise.network(obj, x, meta.network = TRUE, cutoff = cut,
                       radius = c(.4,.7,.85,1), pie.radius = .04,
                       text.size = 3) 
p <- centerTitle(p, paste0('Regulatory network - Metacell ', x))
devPlot(p)

