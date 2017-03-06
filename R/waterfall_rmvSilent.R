#' Silent Mutation Removal
#'
#' Subset a MAF file setting keeping only sample information if a mutation
#' is silent
#' @name waterfall_rmvSilent
#' @param x a data frame with columns 'sample', 'gene', 'trv_type'
#' @return a subset data frame

waterfall_rmvSilent <- function(x)
{
    message("Removing silent mutations...")
    
    # Index and remove those rows which contain silent mutations
    x[which(toupper(x$trv_type) == toupper('silent')), c('gene')] <- NA
    if("label" %in% colnames(x)){
        x[which(toupper(x$trv_type) == toupper('silent')), c('label')] <- NA
    }  
    x[which(toupper(x$trv_type) == toupper('silent')), c('trv_type')] <- NA

    return(x)
}
