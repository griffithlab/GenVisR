#' Assign NA samples a gene
#'
#' Replace NA values in a gene column with the top gene name
#' @name waterfall_NA2gene
#' @param x a data frame in anno format
#' @return a data frame with NA values in a gene column coerced to the top gene
#' name

waterfall_NA2gene <- function(x)
{
    # Get The gene with the most Mutations and add the NA samples to that gene
    # (ensures that the NAs are added in as gene with most mutations will always
    # be plotted)
    top_gene <- na.omit(rev(x$gene))[1]
    x$gene <- replace(x$gene, is.na(x$gene), top_gene)
    return(x)
}
