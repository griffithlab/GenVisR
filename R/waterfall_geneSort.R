#' sort MAF File by gene
#' 
#' order a MAF file ranking genes with more mutations higher
#' @name waterfall_geneSort
#' @param x a data frame in annotation format
#' @return a data frame sorted such that genes with more mutations are ordered 
#' first

waterfall_geneSort <- function(x)
{
    gene_mutation_table <- table(x[,c('gene', 'trv_type')])
    gene_order <- names(sort(rowSums(gene_mutation_table)))
    return(gene_order)
}