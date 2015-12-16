#' Check input to lohView
#'
#' Perform data quality checks on input supplied to lohView
#' @name lohView_qual
#' @param x object of class data frame with columns 'chromosome', 'position',
#' 'n_vaf', 't_vaf', 'sample'
#' @param y object of class data frame with columns 'chromosome', 'start',
#' 'end' specifying chromosomal boundaries for a genome assembly
#' (required if genome is not specified)
#' @param genome character string specifying the genome assembly from which
#' input data is based
#' @return list of inputs passing basic quality controls

lohView_qual <- function(x, y, genome)
{
    # check input data to x
    if(!is.data.frame(x))
    {
        stop("Did not detect a data frame for input to x")
    }

    # check that correct columns are supplied in x
    x_col <- c('chromosome', 'position', 'n_vaf', 't_vaf', 'sample')
    if(!all(x_col %in% colnames(x)))
    {
        stop('Did not detect required column names in x, required columns are: '
             , paste0(x_col, sep="\t"))
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
    } else {
        
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

    return(list(x, y))
}
