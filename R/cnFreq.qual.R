#' check input to cnFreq
#' 
#' Perform a data quality check for inputs to cnFreq
#' @name cnFreq.qual
#' @param x a data frame with columns chromosome, start, end, gain, and loss, or chromosome, start, end, segmean, and sample
#' @return list containing data frame passing quality checks and the type of plot (proportional losses/gain or frequency of losses/gains)

cnFreq.qual <- function(x)
{
    # Check that x is a data frame
    if(!is.data.frame(x))
    {
        warning("Did not detect a data frame in x, attempting to coerce")
        x <- as.data.frame(x)
        x <- droplevels(x)
    }
  
    # Check that x has at least 1 row
    if(nrow(x) < 1)
        stop("x needs at least one row")
  
    # Check that columns have the correct headers
    plotType <- NULL
    if(all(c('chromosome','start','end','gain','loss') %in% colnames(x))) {
        plotType <- "prop"
    } else if(all(c('chromosome','start','end','segmean','sample') %in% colnames(x))) {
        plotType <- "freq"
    } else {
        stop("Did not detect correct columns in x")
    }

    # Check chromosome column in x
    if(!all(grepl("^chr", x$chromosome)))
    {
        message("did not detect the prefix chr in the chromosome column of x, adding prefix")
        x$chromosome <- paste0("chr", x$chromosome)
    } else if(all(grepl("^chr", x$chromosome)))
    {
        #message("Detected chr in the chromosome column of x, proceeding")
    } else {
        stop("Detected unknown or mixed prefixes in the chromosome column of x, should either have a chr prefix or none at all")
    }
  
    # Make sure the chromosome column is of class factor
    x$chromosome <- as.factor(x$chromosome)
  
    if(plotType=="prop") 
    {
        # Check that all proportions are between 0 and 1
    if(any(x$gain<0) | any(x$gain>1) | any(x$loss<0) | any(x$loss>1))
        stop("Gain and loss columns must be in the range [0,1]")
      
    # Check that no proportions add up to more than 1 for the same window
    tmpsum <- apply(x[,c("gain","loss")], 1, sum, na.rm=TRUE)
    if(any(tmpsum>1))
        warning(paste("The proportions of gain + loss sums to greater than 1 for",sum(tmpsum>1),"elements"))
    }
  
    # Make sure that columns are the correct data type
    if(!all(x$start == as.numeric(as.character(x$start))))
        stop("The start column is not numeric")
    if(!all(x$end == as.numeric(as.character(x$end))))
        stop("The end column is not numeric")
    if(plotType=="freq")
    {
        if(!all(x$segmean == as.numeric(as.character(x$segmean))))
            stop("The segmean column is not numeric")
    }
  
    return(list(x,plotType))
}