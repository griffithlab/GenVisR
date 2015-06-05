#' Silent Mutation Removal
#' 
#' Subset a MAF file setting keeping only sample information if a mutation is silent
#' @name mutation_silent_rmv
#' @param x a data frame with columns 'sample', 'gene', 'trv_type'
#' @return a subset data frame

mutation_silent_rmv <- function(x)
{
  # Index and remove those rows which contain silent mutations
  x[which(toupper(x$trv_type) == toupper('silent') | toupper(x$trv_type) == toupper('3_prime_untranslated_region') | toupper(x$trv_type) == toupper('intronic') | toupper(x$trv_type) == toupper('5_prime_untranslated_region')), c('gene')] <- NA
  x[which(toupper(x$trv_type) == toupper('silent') | toupper(x$trv_type) == toupper('3_prime_untranslated_region') | toupper(x$trv_type) == toupper('intronic') | toupper(x$trv_type) == toupper('5_prime_untranslated_region')), c('trv_type')] <- NA
  
  return(x)
}