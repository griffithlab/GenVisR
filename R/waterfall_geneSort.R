#' sort waterfall file by gene
#'
#' order a waterfall file ranking genes with more mutations higher if a gene
#' order is unspecified.
#' @name waterfall_geneSort
#' @param x Data frame with columns names "gene", "trv_type".
#' @param geneOrder Character vector specifying the order in which to plot
#' genes.
#' @return Character vector of ordered genes

waterfall_geneSort <- function(x, geneOrder=NULL)
{
    if(!is.null(geneOrder))
    {
        geneOrder <- as.character(unique(geneOrder))
        # if there are any genes in geneOrder not in x, remove those
        gene_to_rmv <- geneOrder[!geneOrder %in% unique(x$gene)]
        if(length(gene_to_rmv) == 0 & length(geneOrder) >= 1)
        {
            return(rev(geneOrder))
        } else if(length(gene_to_rmv) == length(geneOrder)) {
            memo <- paste0("Did not find any genes supplied to the parameter",
                           " geneOrder in the input supplied to x, perhaps ",
                           " there is a case issue?")
            warning(memo)
        } else if(length(gene_to_rmv) != length(geneOrder)) {
            memo <- paste0("The following genes were not found in the input", 
                           " supplied to parameter x: ", toString(gene_to_rmv),
                           ", removing these from geneOrder!")
            warning(memo)
            gene_order <- geneOrder[geneOrder %in% unique(x$gene)]
            return(rev(gene_order))
        }
    }
    
    # order based on mutation frequency
    gene_mutation_table <- table(x[,c('gene', 'trv_type')])
    gene_order <- names(sort(rowSums(gene_mutation_table)))
    return(gene_order)
}
