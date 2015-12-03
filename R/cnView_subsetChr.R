#' subset based on chr
#'
#' given a data frame subset out specific a chromosome
#' @name cnView_subsetChr
#' @param x a data frame with columns chromosome
#' @param chr character string specifying UCSC chromosome to subset on
#' @return object of class data frame

cnView_subsetChr <- function(x, chr)
{
    # subset data frame based on value specified in chr argument
    if(chr == 'all')
    {
        return(x)
    }
    else if(any(chr == levels(x$chromosome)))
    {
        x <- x[x$chromosome == chr,]
        return(x)
    } else {
        memo <- paste0("Argument supplied to chr does not match levels found",
                       " in the chromosome column of argument supplied to x")
        stop(memo)
    }
}
