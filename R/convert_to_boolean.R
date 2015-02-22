#' Convert counts to a boolean
#' 
#' Convert a data frame such that mutation type are converted to boolean instead of counts 
#' @name convert_to_boolean
#' @param data_frame a data frame in MAF format
#' @return a data frame with values converted to boolean

convert_to_boolean <- function(x)
{
  
  ###########################################################################################################
  ################ Function to convert counts to a boolean based on if there or not #########################
  ###########################################################################################################
  
  # if greater than 1 return 1
  i = which(x > 1)
  x[i] = 1
  return(x)
}