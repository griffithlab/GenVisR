#' subset based on chr
#'
#' given a data frame subset out specific a chromosome
#' @name multi_subsetChr
#' @param x a data frame with columns chromosome
#' @param chr character string specifying UCSC chromosome to subset on
#' @return object of class data frame
#' @noRd

multi_subsetChr <- function(x, chr)
{
    # if all is specified nothing to subset, return dataframe as is
    if(toupper(chr) == 'ALL')
    {
        return(x)
    }
    
    # if chromosome column is not a factor coerce to one
    if(class(x$chromosome) != "factor")
    {
        memo <- paste0("chromosome column supplied is not a factor, attempting",
                       " to coerce...")
        warning(memo)
        
        x$chromosome <- as.factor(x$chromosome)
    }
    
    # subset data frame based on value specified in chr argument
    if(any(chr == levels(x$chromosome)))
    {
        x <- x[x$chromosome == chr,]
        return(x)
    } else {
        memo <- paste0("Argument supplied to chr does not match levels found",
                       " in the chromosome column of argument supplied to x")
        stop(memo)
    }
}
