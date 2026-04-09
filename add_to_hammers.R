scColPairRatio <- function (scObj, col1, col2, sigDigits = 2){
    df <- scColPairPercs(scObj, col1, col2)
    totals <- vapply(unique(df[, 2]), function(x) 
        sum(df[df[, 2] == x, ]$n), integer(1))
    total <- sum(totals)
    df$totalPerc <- round(totals[df[, 2]] / total * 100, sigDigits)
    df$ratio <- df$perc / df$totalPerc
    return(df)
}