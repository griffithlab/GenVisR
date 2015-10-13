#' Construct CN cohort plot
#'
#' given a data frame construct a plot to display CN information for a group of samples
#' @name cnSpec.qual
#' @param x object of class data frame containing columns Chromosome, Start, Stop, SegMean, Sample.name
#' @param y object of class data frame containing user supplied chromosome locations
#' @param genome character string specifying a user supplied genome
#' @return character string specifying input passed quality check

cnSpec.qual <- function(x, y, genome)
{
  # Check for internet connection or for a value to y
    if(is.null(y))
    {
    # Check for internet connectivity if genome is not preloaded
        preloaded <- c("hg38", "hg19", "mm10", "mm9", "rn5")
        if(!any(genome == preloaded))
        {
            if(!is.character(getURL("www.google.com")))
            {
                stop("Did not detect an internet connection or a preloaded genome, check internet connectivity or supply a value to y")
            }
        }
    }

    if(!is.null(y))
    {
    # Check that y is a data frame
        if(!is.data.frame(y))
        {
            message("y is not a data frame, attempting to coerce")
            y <- as.data.frame(y)
        }

    # Check column names of y
        if(!all(c('chromosome', 'start', 'end') %in% colnames(y)))
        {
            stop("Did not detect correct columns in y, missing one of chromosome, start, end")
        }
    }

  # Check that x is a data frame
    if(!is.data.frame(x))
    {
        message("x is not a data frame, attempting to coerce")
        x <- as.data.frame(x)
    }

  # Check column names of x
    if(!all(c('chromosome', 'start', 'end', 'segmean', 'sample') %in% colnames(x)))
    {
        stop("Did not detect correct columns in x, missing one of chromosome, start, end, segmean, sample")
    }

    return(list(x, y))
}
