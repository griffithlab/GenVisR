#' check input to cnFreq
#'
#' Perform a data quality check for input to cnFreq
#' @name cnFreq_qual
#' @param x a data frame with columns chromosome, start, end, segmean, and sample
#' @return  data frame passing quality checks 
#' @noRd

cnFreq_qual <- function(x)
{
    # Check that x is a data frame
    if(!is.data.frame(x)){
        memo <- paste0("Did not detect a data frame in argument supplied",
                       " to x... attempting to coerce")
        warning(memo)
        x <- as.data.frame(x)
        x <- droplevels(x)
    }

    # Check that x has at least 1 row
    if(nrow(x) < 1){
        memo <- paste0("x needs at least one row")
        stop(memo)
    }
    
    # remove any NA values in the data
    if(any(is.na(x))){
        na_rows_removed <- nrow(x) - nrow(na.omit(x))
        memo <- paste0("Removing ", na_rows_removed, " rows containing NA values")
        message(memo)
        x <- na.omit(x)
    }

    if(all(c('chromosome', 'start','end', 'segmean', 'sample') %in% colnames(x))){
        
        # make sure columns are of the correct type
        x$chromosome <- as.factor(x$chromosome)
        x$start <- as.integer(as.character(x$start))
        x$end <- as.integer(as.character(x$end))
        x$segmean <- as.numeric(as.character(x$segmean))
        x$sample <- as.factor(x$sample)
        
        # make sure windows are consistent if not disjoin them
        tmp <- split(x, x$sample)
        tmp_vec <- tmp[[1]]$end
        if(any(!unlist(sapply(tmp, function(x) x[,"end"] %in% tmp_vec), use.names=F))){
            memo <- paste0("Did not detect identical genomic segments for all samples",
                           " ...Performing disjoin operation")
            message(memo) 
            
            # here we split the DF up in an attempt to avoid complaints that lists are to large
            x <- split(x, f=x$chromosome)
            x <- lapply(x, cnFreq_disjoin)

            
            x <- do.call(rbind, x)
        }
        rm(tmp)
        rm(tmp_vec)
    } else {
        memo <- paste0("Did not detect correct columns in argument supplied",
                       " to x!")
        stop(memo)
    }

    # Check chromosome column in x
    if(!all(grepl("^chr", x$chromosome))){
        memo <- paste0("Did not detect the prefix \"chr\" in the chromosome",
                       " column of x... adding prefix")
        message(memo)
        x$chromosome <- paste0("chr", x$chromosome)
        x$chromosome <- as.factor(x$chromosome)
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

    return(x)
}
