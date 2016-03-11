#' mutation sample cutoff gene based
#'
#' Subset a internal mutSpec file keeping only samples within the specified gene
#'  list
#' @name waterfall_geneAlt
#' @param x a data frame in long format with columns 'gene', 'trv_type'
#' @param genes character vector listing genes to plot
#' @return a subset data frame

waterfall_geneAlt <- function(x, genes)
{
    message("Removing genes not in: ", toString(genes))
    # Perform quality checks
    if(typeof(genes) != 'character' & class(genes) != 'character')
    {
        memo <- paste0("argument supplied to plotGenes is not a character ",
                       "vector, attempting to coerce")
        warning(memo)
        genes <- as.character(genes)
    }

    if(!all(toupper(genes) %in% toupper(x$gene)))
    {
        memo <- paste0("genes supplied in plotGenes contains an element not ",
                       "found in x or it's subsequent subsets")
        warning(memo)
    }
    genes <- c(genes, NA)
    x <- x[(toupper(x$gene) %in% toupper(genes)), ]

    return(x)
}
