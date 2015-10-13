#' Check input to lolliplot
#'
#' Perform Basic quality checks for lolliplot input
#' @name lolliplot.qual
#' @param x object of class data frame containing columns transcript_name, gene,
#' and amino_acid_change and rows denoting mutations
#' @param y object of class data frame containing columns transcript_name, and
#' amino_acid_change and rows denoting mutations
#' @return character string specifying input passed quality check

lolliplot.qual <- function(x, y)
{
    # Check that x and y are a data frame, if not attempt to coerce
    if(!is.data.frame(x))
    {
        message(x, " is not a data frame, attempting to coerce")
        x <- as.data.frame(x)
        x <- droplevels(x)
    }

    if(!is.null(y))
    {
        if(!is.data.frame(y))
        {
            message(y, "is not a data frame, attempting to coerce")
            y <- as.data.frame(y)
            y <- droplevels(y)
        }
    }

    # Check for internet connectivity
    if(!is.character(getURL("www.google.com")))
    {
        stop("Did not detect an internet connection,
             check internet connectivity")
    }

    # Check for correct columns in x
    if(!all(c('transcript_name', 'gene', 'amino_acid_change') %in% colnames(x)))
    {
        stop("Did not detect correct columns in x,
             missing one of transcript_name, gene, amino_acid_change")
    }

    # Check that "transcript_name" in x contains only 1 transcript
    if(length(unique(x$transcript_name)) != 1)
    {
        stop("Detected more than 1 transcript in ", x)
    }

    # Check for correct columns in y
    if(!is.null(y))
    {
        if(!all(c('transcript_name', 'amino_acid_change') %in% colnames(y)))
        {
            stop("Did not detect correct columns in y, missing one of
                 transcript_name, amino_acid_change")
        }
    }

    return(list(x, y))
}
