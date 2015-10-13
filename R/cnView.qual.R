#' check input to cnView
#'
#' Perform a data quality check for inputs to cnView
#' @name cnView.qual
#' @param x a data frame with columns chromosome, coordinate, cn, p_value
#' @param y a data frame with columns "chrom", "chromStart", "chromEnd", "name", "gieStain"
#' @param genome character string specifying UCSC genome to use
#' @return a list of data frames passing quality checks

cnView.qual <- function(x, y, genome)
{
  # Check input to x
    if(!is.data.frame(x))
    {
        warning("did not detect a data frame in x, attempting to coerce")
        x <- as.data.frame(x)
        x <- droplevels(x)
    }

    if(!all(c('chromosome', 'coordinate', 'cn') %in% colnames(x)))
    {
        stop("Did not detect correct columns in x, missing one of chromosome, coordinate, cn, p_value")
    }

  # Check chromosome column in x
    if(!all(grepl("^chr", x$chromosome)))
    {
        message("did not detect the prefix chr in the chromosome column of x, adding prefix")
        x$chromosome <- paste0("chr", x$chromosome)
    } else if(all(grepl("^chr", x$chromosome)))
    {
        message("detected chr in the chromosome column of x, proceeding")
    } else {
        stop("Detected unknown or mixed prefixes in the chromosome column of x, should either have a chr prefix or none at all")
    }

  # make sure the chromosome column is of class factor
    x$chromosome <- as.factor(x$chromosome)

  # Check y data
    preloaded <- c("hg38", "hg19", "mm10", "mm9", "rn5")
    if(!is.null(y))
    {
        if(!is.data.frame(y))
        {
            warning("did not detect a data frame in y, attempting to coerce")
            y <- as.data.frame(y)
            y <- droplevels(y)
        }

        if(!all(c("chrom", "chromStart", "chromEnd", "name", "gieStain") %in% colnames(y)))
        {
            stop("Did not detect correct columns in y, missing one of chrom, chromStart, chromEnd, name, gieStain")
        }
    } else if(any(genome == preloaded)){
    # Do nothing here, just a control structure
    } else{
        if(!is.character(getURL("www.google.com")))
        {
            stop("Did not detect an internet connection, this is required if y is not specified and genome is not preloaded")
        }
    }

    return(list(x, y))
}
