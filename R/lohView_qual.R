#' check input to lohView
#'
#' Perform a data quality check for inputs to lohView
#' @name lohView_qual
#' @param x object of class data frame with rows representing germline calls.
#' The data frame must contain columns with the following names "chromosome",
#' "position", "n_vaf", "t_vaf", "sample".
#' @param y a data frame with columns "chrom", "chromStart", "chromEnd", "name",
#' "gieStain"
#' @param genome character string specifying UCSC genome to use
#' @return a list of data frames passing quality checks

lohView_qual <- function(x, y, genome)
{
    ############################## Check input to x #####################
    if(!is.data.frame(x))
    {
        memo <- paste0("Argument supplied to x does not appear to be a",
                       " data frame... attempting to coerce.")
        warning(memo)
        x <- as.data.frame(x)
        x <- droplevels(x)
    }
    
    if(!all(c('chromosome', 'position', 'n_vaf', 't_vaf', 'sample') %in% colnames(x)))
    {
        memo <- paste0("Did not detect correct columns in argument supplied",
                       " to x. missing one of \"chromosome\", \"position\", ",
                       "\"n_vaf\", \"t_vaf\",, \"sample\"")
        stop(memo)
    }
    
    # Check chromosome column in x
    if(!any(grepl("^chr", x$chromosome)))
    {
        memo <- paste0("Did not detect the prefix chr in the chromosome column",
                       " of x... adding prefix")
        message(memo)
        x$chromosome <- paste0("chr", x$chromosome)
    } else if(all(grepl("^chr", x$chromosome))) {
        memo <- paste0("Detected chr in the chromosome column of x...",
                       " proceeding")
        message(memo)
    } else {
        memo <- paste0("Detected mixed prefixes in the chromosome",
                       " column of x... should either be chr or none i.e. ",
                       "chr1 or 1")
        stop(memo)
    }
    
    # make sure the chromosome columns are of proper class
    x$chromosome <- as.factor(x$chromosome)
    x$position <- as.numeric(as.character(x$position))
    x$n_vaf <- as.numeric(as.character(x$n_vaf))
    x$t_vaf <- as.numeric(as.character(x$t_vaf))
    x$sample <- as.factor(x$sample)
    
    ################### Check input to y ###############################
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
                           " terms to attempt query to UCSC mySQL database.",
                           "Alternativly supply a value to y.")
            warning(memo)
        }
    }
    
    return(list(x, y))
}
