#' Assign NA samples a gene
#'
#' Replace NA values in a gene column with the top gene name
#' @name waterfall_NA2gene
#' @param x a data frame in anno format
#' @return a data frame with NA values in a gene column coerced to the top gene
#' name
#' @importFrom stats na.omit

waterfall_NA2gene <- function(x)
{
    # Get The gene with the most Mutations and add the NA samples to that gene
    # (ensures that the NAs are added in as gene with most mutations will always
    # be plotted) i.e. makes sure that samples are plotted, happens with rmvSilent param
    
    # find top gene
    top_gene <- stats::na.omit(rev(x$gene))[1]
    
    # set the trv_type to NA if gene is NA (makes sure that wile a sample is plotted the cell is empty)
    x$trv_type <- replace(x$trv_type, is.na(x$gene), NA)
    x$gene <- replace(x$gene, is.na(x$gene), top_gene)
    
    return(x)
}
