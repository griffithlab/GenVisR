#' Check input to lolliplot
#'
#' Perform Basic quality checks for lolliplot input
#' @name lolliplot_qual
#' @param x object of class data frame containing columns transcript_name, gene,
#' and amino_acid_change and rows denoting mutations
#' @param y object of class data frame containing columns transcript_name, and
#' amino_acid_change and rows denoting mutations
#' @param z Object of class data frame containing columns "description", "start",
#'  "stop" specifying gene regions to highlight
#' @return objects passing basic quality checks

lolliplot_qual <- function(x, y, z)
{
    # Check input to x
    if(!is.data.frame(x))
    {
        message("Input to x is not a data frame, attempting to coerce")
        x <- as.data.frame(x)
        x <- droplevels(x)
    }
    
    # Check for correct columns in x
    if(!all(c('transcript_name', 'gene', 'amino_acid_change') %in% colnames(x)))
    {
        stop("Did not detect correct columns in x,
             missing one of transcript_name, gene, amino_acid_change")
    }
    
    # Make sure columns in x are of the proper class
    x$transcript_name <- as.factor(x$transcript_name)
    x$gene <- as.factor(x$gene)
    x$amino_acid_change <- as.factor(x$amino_acid_change)
    
    # Check that "transcript_name" in x contains only 1 transcript
    if(length(unique(x$transcript_name)) != 1)
    {
        stop("Detected more than 1 transcript in input to x")
    }
    
    # Check input to y
    if(!is.null(y))
    {
        # is y a data frame?
        if(!is.data.frame(y))
        {
            message(y, "is not a data frame, attempting to coerce")
            y <- as.data.frame(y)
            y <- droplevels(y)
        }
        
        # does y have correct columns?
        if(!all(c('transcript_name', 'amino_acid_change') %in% colnames(y)))
        {
            stop("Did not detect correct columns in y, missing one of
                 transcript_name, amino_acid_change")
        }
        
        # make sure columns in y are of proper class
        y$transcript_name <- as.factor(y$transcript_name)
        y$amino_acid_change <- as.factor(y$amino_acid_change)
    }

    # Check input to z
    if(!is.null(z))
    {
        # is z a data frame?
        if(!is.data.frame(z))
        {
            memo <- paste0("Input to z is not a data frame",
                           ", attempting to coerce")
            warning(memo)
            z <- as.data.frame(z)
            z <- droplevels(z)
        }
        
        # does z contain correct columns
        if(!all(c("description", "start", "stop") %in% colnames(z)))
        {
            memo <- paste0("Did not detect correct columns in input to z, ",
                           "missing one of \"description\", \"start\",",
                           " \"stop\"")
            stop(memo)
        }
        
        # make sure column class of z are of proper type
        z$description <- as.factor(z$description)
        z$start <- as.numeric(as.character(z$start))
        z$stop <- as.numeric(as.character(z$stop))
    }

    return(list(x, y, z))
}
