seuratObj <- qs_read('annotatedSeurat.qs2')

df <- scColPairPercs(seuratObj, 'orig.ident', 'celltype')
df$orig.ident <- factor(df$orig.ident, 
                        levels=c('Control', '0h', '12h', '24h'))
df$celltype <- factor(df$celltype,
                      levels=c('Muller glial cells', 
                               'Retinal progenitor cells',
                               'Gliogenic progenitors',
                               'Neurogenic progenitors',
                               'Proliferative RPC',
                               'Photoreceptor precursors',
                               'Rods',
                               'Bipolar cells',
                               'Amacrine cells',
                               'Cones',
                               'ISL1+ retinal ganglion cells',
                               'HTR2C+ retinal ganglion cells',
                               'NEUROD6+ retinal ganglion cells',
                               'ATF5+ horizontal cells',
                               'ONECUT1+ horizontal cells',
                               'RPE cells'))


celltypeCols <- c('Muller glial cells' = rgb(247/255,147/255,30/255), 
                  'Retinal progenitor cells' = rgb(210/255,0/255,0/255), 
                  'Gliogenic progenitors' = rgb(248/255,115/255,106/255), 
                  'Neurogenic progenitors' = rgb(180/255,0/255,150/255),
                  'Proliferative RPC' = rgb(200/255,50/255,15/255), 
                  'Photoreceptor precursors' = rgb(163/255,100/255,140/255), 
                  'Rods' = rgb(169/255,169/255,169/255), 
                  'Bipolar cells' = rgb(163/255,165/255,0/255), 
                  'Amacrine cells' = rgb(157/255,115/255,194/255), 
                  'Cones' = rgb(230/255,134/255,201/255), 
                  'ISL1+ retinal ganglion cells' = rgb(140/255,198/255,63/255),
                  'HTR2C+ retinal ganglion cells' = rgb(192/255,193/255,48/255), 
                  'NEUROD6+ retinal ganglion cells' = rgb(74/255,176/255,91/255), 
                  'ATF5+ horizontal cells' = rgb(97/255,156/255,255/255), 
                  'ONECUT1+ horizontal cells' = rgb(3/255,161/255,198/255), 
                  'RPE cells'= rgb(129/255,70/255,58/255))

p <- ggplot(df, aes(x=orig.ident, 
                    y=perc, 
                    fill=celltype, 
                    stratum=celltype, 
                    alluvium=celltype)) +                 
    geom_col(width = 0.6,color="white")+
    geom_flow(width = 0.6, alpha=0.2, knot.pos = 0.1, color="white") +  
    theme_void()+ 
    theme(axis.text.x=element_text(size=10,vjust = 5)) + scale_fill_discrete(palette=celltypeCols)
p <- centerTitle(p, 'Changes in cell type representation')

p

devPlot(p)

################################################################################

df <- scColPairRatio(seuratObj, 'orig.ident', 'celltype')
df <- df[, c(1, 2, 6)]

df <- reshape2::melt(df)
df$significant <- case_when(
    df$value < 1 ~ "±",
    df$value > 1 & df$value < 1.5 ~ "+",
    df$value > 1.5 & df$value < 2 ~ "++",
    df$value > 2 ~ "+++",
)

df$orig.ident <- factor(df$orig.ident, 
                        levels=c('Control', '0h', '12h', '24h'))
df$celltype <- factor(df$celltype,
                      levels=rev(c('Muller glial cells', 
                               'Retinal progenitor cells',
                               'Gliogenic progenitors',
                               'Neurogenic progenitors',
                               'Proliferative RPC',
                               'Photoreceptor precursors',
                               'Rods',
                               'Bipolar cells',
                               'Amacrine cells',
                               'Cones',
                               'ISL1+ retinal ganglion cells',
                               'HTR2C+ retinal ganglion cells',
                               'NEUROD6+ retinal ganglion cells',
                               'ATF5+ horizontal cells',
                               'ONECUT1+ horizontal cells',
                               'RPE cells')))

p <- ggplot(df, aes(x = orig.ident, y = celltype, label = significant)) +
    geom_tile(aes(fill = value), color="white", linewidth=0.5) +
    geom_text(color = "black") + 
    theme_classic() +
    scale_fill_viridis_c(option = "plasma", direction = 1, alpha = 1) +
    theme(axis.title = element_blank()) +
    labs(fill='Observed over\nexpected ratio')
p <- centerTitle(p, 'Observed over expected ratio')

devPlot(p)
