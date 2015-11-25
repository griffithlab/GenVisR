#' check input to cnView
#'
#' Perform a data quality check for inputs to cnView
#' @name cnView_qual
#' @param x a data frame with columns chromosome, coordinate, cn
#' @param y a data frame with columns "chrom", "chromStart", "chromEnd", "name",
#' "gieStain"
#' @param genome character string specifying UCSC genome to use
#' @return a list of data frames passing quality checks

cnView_qual <- function(x, y, genome)
{
    # Check input to x
    if(!is.data.frame(x))
    {
        memo <- paste0("Argument supplied to x does not appear to be a",
                       " data frame... attempting to coerce.")
        warning(memo)
        x <- as.data.frame(x)
        x <- droplevels(x)
    }

    if(!all(c('chromosome', 'coordinate', 'cn') %in% colnames(x)))
    {
        memo <- paste0("Did not detect correct columns in argument supplied",
                       " to x. missing one of \"chromosome\", \"coordinate\", ",
                       "\"cn\", \"p_value\"")
        stop(memo)
    }

    # Check chromosome column in x
    if(!all(grepl("^chr", x$chromosome)))
    {
        memo <- paste0("Did not detect the prefix chr in the chromosome column",
                       " of x... adding prefix")
        message(memo)
        x$chromosome <- paste0("chr", x$chromosome)
    } else if(all(grepl("^chr", x$chromosome))) {
        memo <- paste0("detected chr in the chromosome column of x...",
                       "proceeding")
        message(memo)
    } else {
        memo <- paste0("Detected unknown or mixed prefixes in the chromosome",
                       " column of x... should either be chr or none i.e. ",
                       "chr1 or 1")
        stop(memo)
    }

    # make sure the chromosome column is of class factor
    x$chromosome <- as.factor(x$chromosome)

    # Check y data
    preloaded <- c("hg38", "hg19", "mm10", "mm9", "rn5")
    if(!is.null(y))
    {
        if(!is.data.frame(y))
        {
            memo <- paste0("Agrument supplied to y does not appear to be a",
                           "data frame... attempting to coerce")
            message(memo)
            y <- as.data.frame(y)
            y <- droplevels(y)
        }

        if(!all(c("chrom",
                  "chromStart",
                  "chromEnd",
                  "name",
                  "gieStain") %in% colnames(y)))
        {
            memo <- paste0("Did not detect correct columns in y... missing one",
                           " of \"chrom\", \"chromStart\", \"chromEnd\",",
                           " \"name\", \"gieStain\"")
            stop(memo)
        }
    } else if(!any(genome == preloaded)) {
        
        if(grepl("NCBI|GRC|RGSC|BROAD|BAYLOR|WUGSC", genome, ignore.case=TRUE))
        {
            memo <- paste0("Detected a genome that does not appear to be,",
                           "in UCSC terms, please specify a genome in UCSC",
                           " terms to attempt query to UCSC mySQL databae.",
                           "Alternativly supply a value to y.")
            warning(memo)
        }
    }

    return(list(x, y))
}
