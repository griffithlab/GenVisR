#' sort MAF File by gene
#' 
#' order a MAF file ranking genes with more mutations higher
#' @name gene_sort
#' @param data_frame a data frame in MAF format
#' @return a data frame sorted such that genes with more mutations are ordered first

gene_sort <- function(data_frame)
{
  
  #########################################################################################################################################
  ####### Create a table of genes and their mutation status, sum the mutation totals for each gene and sort giving a vector of gene #######
  ####### names sorted by mutation frequency, note that the function expects a data frame in long format with columns 'gene', 'trv_type' ##
  #########################################################################################################################################
  
  gene_mutation_table <- table(data_frame[,c('gene', 'trv_type')])
  gene_order <- names(sort(rowSums(gene_mutation_table)))
  return(gene_order)
}