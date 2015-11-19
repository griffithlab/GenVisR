#' Perform quality control on genCov data
#'
#' Ensure data input into genCov is of the proper type and format
#' @name genCov.qual
#' @param x named list containing data frames with columns end and cov
#' @param txdb A TxDb object for a genome
#' @param gr A Granges object specifying a region of interest
#' @param genome Object of class BSgenome specifying the genome
#' @return a list of objects passing basic quality control

genCov.qual <- function(x, txdb, gr, genome)
{
  # Check input to x
    if(!is.data.frame(x) && !is.list(x))
    {
        warning("x does not appear to be a list, attempting to coerce")
        x <- as.list(x)
    }else if(is.data.frame(x)){
      x<-list(c(x))
      names(x)<-c('Sample 1')
    }

    if(!all(as.logical(lapply(x, is.data.frame))))
    {
        warning("all arguments in x are not of type data frame, attempting to coerce")
        x <- lapply(x, as.data.frame)
    }

    if(!all(as.logical(lapply(x, function(x) all(c('end', 'cov') %in% colnames(x))))))
    {
        stop("one or more elements of x are missing the columns end or cov")
    }

  # Check the TxDb object
    if(class(txdb)[1] != 'TxDb')
    {
        message("txdb does not appear to be an object of class TxDb, if you believe this message is in error press [enter]")
        line <- readline()
    }

  # Check the Genomic ranges object
    if(class(gr)[1] != 'GRanges')
    {
        message("gr does not appear to be an object of class GRanges, if you believe this message is in error press [enter]")
        line <- readline()
    }

  # Check the biostrings object
    if(class(genome)[1] != 'BSgenome')
    {
        message("genome does not appear to be an object of class BSgenome, if you believe this message is in error press [enter]")
        line <- readline()
    }

    return(list(x, txdb, gr, genome))
}
