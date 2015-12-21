#' check input to cnView
#'
#' Perform a data quality check for inputs to cnView
#' @name cnView_qual
#' @param x a data frame with columns chromosome, coordinate, cn
#' @param y a data frame with columns "chrom", "chromStart", "chromEnd", "name",
#' "gieStain"
#' @param CNscale Character string specifying if copy number calls supplied are
#' relative (i.e.copy neutral == 0) or absolute (i.e. copy neutral ==2). One of 
#' "relative" or "absolute"
#' @param z a data frame with columns chromosome, start, end , segmean
#' @param genome character string specifying UCSC genome to use
#' @return a list of data frames passing quality checks

cnView_qual <- function(x, y, z, genome, CNscale)
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

    if(!all(c('chromosome', 'coordinate', 'cn') %in% colnames(x)))
    {
        memo <- paste0("Did not detect correct columns in argument supplied",
                       " to x. missing one of \"chromosome\", \"coordinate\", ",
                       "\"cn\"")
        stop(memo)
    }
    
    if(CNscale == "absolute")
    {
        # if any cn values are negative something fishy is happening, report
        if(any(x$cn < 0))
        {
            memo <- paste0("Detected negative values in the copy number",
                           " column but CNscale is set to \"absolute\"!")
            warning(memo)
        }
    } else if(CNscale == "relative") {
        next
    } else {
        memo <- paste0("Did not recognize input to parameter CNscale",
                       " please specify one of \"relative\" or \"absolute\"!")
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

    # make sure the chromosome column is of class factor
    x$chromosome <- as.factor(x$chromosome)

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
                           " terms to attempt query to UCSC mySQL databae.",
                           "Alternativly supply a value to y.")
            warning(memo)
        }
    }
    
    ###################### Check input to Z #########################
    if(!is.null(z))
    {
        if(!is.data.frame(z))
        {
            memo <- paste0("Agrument supplied to z does not appear to be a",
                           "data frame... attempting to coerce")
            message(memo)
            z <- as.data.frame(z)
            z <- droplevels(z)
        }
        
        if(!all(c("chromosome",
                  "start",
                  "end",
                  "segmean") %in% colnames(z)))
        {
            memo <- paste0("Did not detect correct columns in z... missing one",
                           " of \"chromosome\", \"start\", \"end\",",
                           " \"segmean\"")
            stop(memo)
        }
        
        # Check chromosome column in z
        if(!any(grepl("^chr", z$chromosome)))
        {
            memo <- paste0("Did not detect the prefix chr in the chromosome",
                           " column of z... adding prefix")
            message(memo)
            z$chromosome <- paste0("chr", z$chromosome)
        } else if(all(grepl("^chr", z$chromosome))) {
            memo <- paste0("detected chr in the chromosome column of z...",
                           "proceeding")
            message(memo)
        } else {
            memo <- paste0("Detected mixed prefixes in the",
                           " chromosome column of z... should either be chr or",
                           "none i.e. chr1 or 1")
            stop(memo)
        }
        
        # make sure the chromosome column is of class factor
        z$chromosome <- as.factor(z$chromosome)
    }

    return(list(x, y, z))
}
