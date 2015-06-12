#' Convert MAF File
#' 
#' Convert columns of a mutation annotation file "MAF" into a format recognizable by internal functions
#' @name MAF_to_anno
#' @param x a data frame in MAF format
#' @return a data frame coerced from MAF to TGI format

MAF_to_anno <- function(x)
{
  ##################################################################################################################
  ############## Function to take a MAF file and coerce it into a format recognizable by other functions ###########
  ##################################################################################################################
  if(!all(c('Tumor_Sample_Barcode', 'Hugo_Symbol', 'Variant_Classification') %in% colnames(x)))
  {
    stop("Did not detect correct column names, check file_type flag")
  }

  x <- x[,c('Tumor_Sample_Barcode', 'Hugo_Symbol', 'Variant_Classification')]
  colnames(x) <- c('sample', 'gene', 'trv_type')
  
  return(x)
}