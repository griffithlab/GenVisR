#' Remove indels
#' 
#' Given a data frame with columns reference and variants remove all indels from data
#' @name rm.indel
#' @param x Object of class data frame containing columns 'reference', 'variant'
#' @return Object of class data frame with indels removed

rm.indel <- function(x)
{
  original_size <- nrow(x)
  
  # Find and remove insertions and deletions
  x <- x[grep('-|0', x$reference, perl=TRUE, invert=TRUE),]
  x <- x[grep('-|0', x$variant, perl=TRUE, invert=TRUE),]
  
  new_size <- nrow(x)
  
  # Print message if indels have been removed
  if(new_size != original_size)
  {
      message("Removed ", original_size - new_size, " indels present in data")
  }
  return(x)
}