#' Construct CN cohort plot
#'
#' given a data frame construct a plot to display CN information for a group of
#' samples
#' @name cnSpec_qual
#' @param x object of class data frame containing columns chromosome, start,
#' stop, segmean, sample
#' @param y object of class data frame containing user supplied chromosome
#' locations
#' @param CNscale Character string specifying if copy number calls supplied are
#' relative (i.e.copy neutral == 0) or absolute (i.e. copy neutral ==2). One of 
#' "relative" or "absolute"
#' @param genome character string specifying a user supplied genome
#' @return character string specifying input passed quality check
#' @noRd

cnSpec_qual <- function(x, y, genome, CNscale)
{
    # Check genome is acceptable name if y is not supplied
    if(is.null(y))
    {
        # Check that genome specified is not the ensembl name
        preloaded <- c("hg38", "hg19", "mm10", "mm9", "rn5")
        if(!any(genome == preloaded))
        {
            if(grepl("NCBI|GRC|RGSC|BROAD|BAYLOR|WUGSC",
                     genome, ignore.case=TRUE))
            {
                memo <- paste0("Detected a genome that does not appear to be,",
                               "in UCSC terms, please specify a genome in UCSC",
                               " terms to attempt query to UCSC mySQL databae.",
                               "Alternativly supply a value to y.")
                warning(memo)
            }
        }
    }

    # Check input to y
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
            memo <- paste0("Did not detect correct column names in y, missing",
                           "one of \"chromosome\", \"start\", \"end\"")
            stop(memo)
        }
        
        # Ensure that columns in data frame are of proper type
        y$chromosome <- as.character(y$chromosome)
        y$start <- as.integer(as.character(y$start))
        y$end <- as.integer(as.character(y$end))
    }

    # Check that x is a data frame
    if(!is.data.frame(x))
    {
        message("x is not a data frame, attempting to coerce")
        x <- as.data.frame(x)
    }

    # Check column names of x
    if(!all(c('chromosome',
              'start',
              'end',
              'segmean',
              'sample') %in% colnames(x)))
    {
        memo <- paste0("Did not detect correct columns in input to x, ",
                       "missing one of \"chromosome\", \"start\", \"end\",",
                       " \"segmean\", \"sample\"!")
        stop(memo)
    }
    
    if(CNscale == "absolute")
    {
        # if any cn values are negative something fishy is happening, report
        if(any(as.numeric(as.character(x$segmean)) < 0))
        {
            memo <- paste0("Detected negative values in the segmean",
                           " column but CNscale is set to \"absolute\"!")
            warning(memo)
        }
    } else if(CNscale == "relative") {
        # Do nothing
    } else {
        memo <- paste0("Did not recognize input to parameter CNscale",
                       " please specify one of \"relative\" or \"absolute\"!")
        stop(memo)
    }
    
    # Ensure that data frame columns are of proper type
    x$chromosome <- as.character(x$chromosome)
    x$start <- as.integer(as.character(x$start))
    x$end <- as.integer(as.character(x$end))
    x$segmean <- as.numeric(as.character(x$segmean))
    if(!is.factor(x$sample))
    {
        x$sample <- as.character(x$sample)
    }

    return(list(x, y))
}
