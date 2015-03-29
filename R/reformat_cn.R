#' reformat CN data frame
#' 
#' given a CN data frame reformat it into a data frame build_CN_plot can understand
#' @name reformat_cn
#' @param data_frame a data frame with columns Chr, Coord, Tumor, Normal, Diff, p_value
#' @return object of class data frame

reformat_cn <- function(data_frame)
{
  colnames(data_frame) <- c('Chr', 'Coord', 'Tumor', 'Normal', 'Diff', 'p_value')
  
  return(data_frame)
}