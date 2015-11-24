#' subset based on chr
#'
#' given a data frame with cytogenetic band locations subset out specific
#' chromosome
#' @name cnView_subsetChr
#' @param x a data frame with columns Chr, Coord, Tumor, Normal, Diff
#' @param chr character string specifying UCSC chromosome to plot one of chr...
#' or all
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
