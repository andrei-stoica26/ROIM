library(ggalluvial)
library(scLang)
library(hammers)

df <- scColPairPercs(seuratObj, 'orig.ident', 'celltype')

p <- ggplot(df, aes(x=orig.ident, 
                    y=perc, 
                    fill=celltype, 
                    stratum=celltype, 
                    alluvium=celltype)) +                 
    geom_col(width = 0.6,color="white")+
    geom_flow(width = 0.6, alpha=0.2, knot.pos = 0.1, color="white") +  
    theme_void()+ 
    theme(axis.text.x=element_text(size=10,vjust = 5))

df <- repAnalysis(seuratObj, 'orig.ident', 'celltype')
pvalRiverPlot(df)



