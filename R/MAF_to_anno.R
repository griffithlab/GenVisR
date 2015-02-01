#' Convert MAF File
#' 
#' Convert columns of a mutation annotation file "MAF" into a format recognizable by internal functions
#' @name MAF_to_anno
#' @param data_frame a data frame in MAF format
#' @return a data frame coerced from MAF to TGI format

MAF_to_anno <- function(data_frame)
{
  ##################################################################################################################
  ############## Function to take a MAF file and coerce it into a format recognizable by other functions ###########
  ##################################################################################################################
  
  data_frame <- data_frame[,c('Tumor_Sample_Barcode', 'Hugo_Symbol', 'Variant_Classification')]
  colnames(data_frame) <- c('sample', 'gene', 'trv_type')
  
  return(data_frame)
}