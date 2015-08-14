#' sort MAF File by gene
#' 
#' order a MAF file ranking genes with more mutations higher
#' @name mutSpec.gene_sort
#' @param data_frame a data frame in MAF format
#' @return a data frame sorted such that genes with more mutations are ordered 
#' first

mutSpec.gene_sort <- function(data_frame)
{
    gene_mutation_table <- table(data_frame[,c('gene', 'trv_type')])
    gene_order <- names(sort(rowSums(gene_mutation_table)))
    return(gene_order)
}