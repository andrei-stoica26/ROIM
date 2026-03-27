newCnetplot <- function(enrichmentResult, title = NULL, nCategories = 10){
    p <- enrichplot::cnetplot(enrichmentResult, showCategory = nCategories,
                              color_category = 'purple',
                              color_item = 'red',
                              color_edge = 'thistle') + NoLegend()
    p <- centerTitle(p, title)
    return(p)
}