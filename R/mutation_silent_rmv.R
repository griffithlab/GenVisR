#' Silent Mutation Removal
#' 
#' Subset a MAF file keeping only samples/genes that have >= 1 non-silent mutations
#' @name mutation_silent_rmv
#' @param x a data frame with columns 'sample', 'gene', 'trv_type'
#' @return a subset data frame

mutation_silent_rmv <- function(x)
{
  # Index and remove those rows which contain silent mutations
  x <- x[which(toupper(x$trv_type) != toupper('silent')),]
  
  return(x)
}