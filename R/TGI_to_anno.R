#' Convert TGI File
#' 
#' Convert columns of a mutation annotation file "TGI" into a format recognizable by internal functions
#' @name TGI_to_anno
#' @param x a data frame in TGI internal format
#' @return a data frame coerced from TGI to internal annotation format

TGI_to_anno <- function(x)
{
  ##################################################################################################################
  ############## Function to take a MAF file and coerce it into a format recognizable by other functions ###########
  ##################################################################################################################
  
  x <- x[,c('sample', 'gene_name', 'trv_type')]
  colnames(x) <- c('sample', 'gene', 'trv_type')
  
  return(x)
}