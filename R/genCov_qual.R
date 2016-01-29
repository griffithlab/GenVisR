#' Perform quality control on genCov data
#'
#' Ensure data input into genCov is of the proper type and format
#' @name genCov_qual
#' @param x named list containing data frames with columns end and cov
#' @param txdb A TxDb object for a genome
#' @param gr A Granges object specifying a region of interest
#' @param genome Object of class BSgenome specifying the genome
#' @return a list of objects passing basic quality control

genCov_qual <- function(x=x, txdb=txdb, gr=gr, genome=genome)
{
    # Check input to x
    if(class(x) != 'list')
    {
        warning("class of x does not appear to be a list, attempting to coerce")
        x <- as.data.frame(x)
        x <- list('Sample 1'=x)
    }

    if(!all(as.logical(lapply(x, is.data.frame))))
    {
        memo <- paste0("all arguments in  are not of type data frame", 
                       "... attempting to coerce")
        warning(memo)
        x <- lapply(x, as.data.frame)
    }

    if(!all(as.logical(lapply(x, function(x) all(c('end', 'cov') %in%
                                                     colnames(x))))))
    {
        memo <- paste0("one or more elements of x are missing column names,",
                       " one of end or cov")
        stop(memo)
    }

    # Check the TxDb object
    if(class(txdb)[1] != 'TxDb')
    {
        memo <- paste0("txdb does not appear to be an object of class TxDb")
        stop(memo)
    }

    # Check the Genomic ranges object
    if(class(gr)[1] != 'GRanges')
    {
        memo <- paste0("gr does not appear to be an object of class GRanges")
        stop(memo)
    }

    # Check the biostrings object
    if(class(genome)[1] != 'BSgenome')
    {
        memo <- paste0("genome does not appear to be an object",
                       "of class BSgenome")
        stop(memo)
    }

    return(list('x'=x, 'txdb'=txdb, 'gr'=gr, 'genome'=genome))
}
