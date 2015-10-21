#' Remove multinucleotide codes
#' 
#' Given a data frame with columns reference and variants remove all
#' multinucleotides from data
#' @name TvTi_rmMnuc
#' @param x Object of class data frame containing columns 'reference', 'variant'
#' @return Object of class data frame with multi nucleotide codes removed

TvTi_rmMnuc <- function(x)
{
  original_size <- nrow(x)
  
  # Find and multi nucleotide codes
  x <- x[grep('[ACGT]{2,}', x$reference, perl=TRUE, invert=TRUE),]
  x <- x[grep('[ACGT]{2,}', x$variant, perl=TRUE, invert=TRUE),]
  
  new_size <- nrow(x)
  
  if(new_size != original_size)
  {
      memo <- paste0("Multi Nucleotide codes are not currently supported, ", 
                     "removed: ", original_size - new_size,
                     " multi-nucleotides present in the data")
      warning(memo) 
  }
  
  return(x)
}