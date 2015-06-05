#' Silent Mutation Removal
#' 
#' Subset a MAF file setting keeping only sample information if a mutation is silent
#' @name mutation_silent_rmv
#' @param x a data frame with columns 'sample', 'gene', 'trv_type'
#' @return a subset data frame

mutation_silent_rmv <- function(x)
{
  # Index and remove those rows which contain silent mutations
  x[which(toupper(x$trv_type) == toupper('silent')), c('gene')] <- NA
  x[which(toupper(x$trv_type) == toupper('silent')), c('trv_type')] <- NA
  
  return(x)
}