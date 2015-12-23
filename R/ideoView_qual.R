#' Check input to ideoView
#' 
#' Check that input to ideoView is properly formatted
#' 
#' @name ideoView_qual
#' @param x a data frame with rows representing cytogenetic bands for a genome.
#' The data frame should have columns "chrom", "chromStart", "chromEnd", "name",
#' "gieStain'.
#' @return data frame

ideoView_qual <- function(x)
{
    # Check that input is a data frame
    if(!is.data.frame(x))
    {
        memo <- paste0("Input does not appear to be a data frame... ",
                       "attempting to coerce.")
        message(memo)
        x <- as.data.frame(x)
    }
    
    # Check that proper columns are present
    col <- c("chrom", "chromStart", "chromEnd", "name", "gieStain")
    if(!all(col %in% colnames(x)))
    {
        memo <- paste0("Did not detect correct columns in input to x, ",
                       "missing one of the following: ", toString(col))
        stop(memo)
    }
    
    # check that the giestain column contains the proper levels
    gie_stain <- c('gneg',
                   'stalk',
                   'acen',
                   'gpos',
                   'gvar',
                   paste0('gpos', 1:100))
    
    if(any(!x$gieStain %in% gie_stain))
    {
        memo <- paste0("Found variable in column gieStain not recognized,",
                       " please ensure variable supplied are one of the",
                       " following: ", toString(gie_stain))
        stop(memo)
    }
    
    return(x)
}