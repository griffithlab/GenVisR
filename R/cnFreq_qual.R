#' check input to cnFreq
#'
#' Perform a data quality check for inputs to cnFreq
#' @name cnFreq_qual
#' @param x a data frame with columns chromosome, start, end, gain, and loss,
#' or chromosome, start, end, segmean, and sample
#' @return list containing data frame passing quality checks and the type
#' of plot (proportional losses/gain or frequency of losses/gains)

cnFreq_qual <- function(x)
{
    # Check that x is a data frame
    if(!is.data.frame(x))
    {
        memo <- paste0("Did not detect a data frame in argument supplied",
                       " to x... attempting to coerce")
        warning(memo)
        x <- as.data.frame(x)
        x <- droplevels(x)
    }

    # Check that x has at least 1 row
    if(nrow(x) < 1)
    {
        memo <- paste0("x needs at least one row")
        stop(memo)
    }

    # Check that columns have the correct headers
    plotType <- NULL
    if(all(c('chromosome',
             'start',
             'end',
             'gain',
             'loss') %in% colnames(x)))
    {
        plotType <- "prop"
        x$chromosome <- as.factor(x$chromosome)
        x$start <- as.integer(as.character(x$start))
        x$end <- as.integer(as.character(x$end))
        x$loss <- as.numeric(as.character(x$loss))
        x$gain <- as.numeric(as.character(x$gain))
    } else if(all(c('chromosome',
                    'start','end',
                    'segmean',
                    'sample') %in% colnames(x))) {
        plotType <- "freq"
        x$chromosome <- as.factor(x$chromosome)
        x$start <- as.integer(as.character(x$start))
        x$end <- as.integer(as.character(x$end))
        x$segmean <- as.numeric(as.character(x$segmean))
        x$sample <- as.factor(x$sample)
        
        # also make sure windows are consistent (this is temporary)
        tmp_vec <- x$end
        tmp <- split(x$sample, x$sample)
        if(!all(sapply(tmp, length) == length(tmp[[1]]))){
            memo <- paste0("Input to x must have consistent windows across samples",
                           " output may be incorrect!")
            warning(memo)
        }
        if(any(!sapply(tmp, function(x) x[,3] %in% tmp_vector))){
            memo <- paste0("Input to x must have consistent windows across samples",
                           " output may be incorrect!")
            warning(memo)            
        }
        rm(tmp)
        rm(tmp_vec)
    } else {
        memo <- paste0("Did not detect correct columns in argument supplied",
                       " to x!")
        stop(memo)
    }

    # Check chromosome column in x
    if(!all(grepl("^chr", x$chromosome)))
    {
        memo <- paste0("Did not detect the prefix \"chr\" in the chromosome",
                       " column of x... adding prefix")
        message(memo)
        x$chromosome <- paste0("chr", x$chromosome)
    } else if(all(grepl("^chr", x$chromosome))) {
        memo <- paste0("Detected \"chr\" in the chromosome column of x...",
                       " proceeding")
        message(memo)
    } else {
        memo <- paste0("Detected unknown or mixed prefixes in the chromosome ",
                       " column of x, should either have a chr prefix or ",
                       "none at all!")
        stop(memo)
    }

    # Make sure the chromosome column is of class factor
    x$chromosome <- as.factor(x$chromosome)

    if(plotType=="prop")
    {
        # Check that all proportions are between 0 and 1
        if(any(x$gain<0) | any(x$gain>1) | any(x$loss<0) | any(x$loss>1))
        {
            stop("Gain and loss columns must be in the range [0,1]")
        }

        # Check that no proportions add up to more than 1 for the same window
        tmpsum <- apply(x[,c("gain","loss")], 1, sum, na.rm=TRUE)
		if(any(round(tmpsum, digits=1) > 1))
        {
            memo <- paste0("The proportions of gain + loss sums to greater ",
                           "than 1 for ", sum(tmpsum>1), " elements!")
            warning(memo)
        }

    }

    # Make sure that columns are the correct data type
    if(!all(x$start == as.numeric(as.character(x$start))))
    {
        stop("The start column is not numeric")
    }

    if(!all(x$end == as.numeric(as.character(x$end))))
    {
        stop("The end column is not numeric")
    }

    if(plotType=="freq")
    {
        if(!all(x$segmean == as.numeric(as.character(x$segmean))))
        {
            stop("The segmean column is not numeric")
        }
    }

    return(list(x,plotType))
}
