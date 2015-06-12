#' Convert MGI File
#' 
#' Convert columns of a mutation annotation file "MGI" into a format recognizable by internal functions
#' @name MGI_to_anno
#' @param x a data frame in MGI internal format
#' @return a data frame coerced from MGI to internal annotation format

MGI_to_anno <- function(x)
{
  ##################################################################################################################
  ############## Function to take a MAF file and coerce it into a format recognizable by other functions ###########
  ##################################################################################################################
  if(!all(c('sample', 'gene_name', 'trv_type') %in% colnames(x)))
  {
    stop("Did not detect correct column names, check file_type flag")
  }
  
  x <- x[,c('sample', 'gene_name', 'trv_type')]
  colnames(x) <- c('sample', 'gene', 'trv_type')
  
  return(x)
}