#' Construct coverage cohort plot
#'
#' given a matrix construct a plot to display coverage as percentage bars for a
#' group of samples
#' @name covBars_qual
#' @param x object of class matrix containing rows for the coverage and columns
#' the sample names
#' @return a list of data frame and color vector
#' @noRd

covBars_qual <- function(x)
{
    # Check that x is a matrix with at least 1 row
    if(!is.matrix(x))
    {
        memo <- paste0("Argument supplied to x is not a matrix... ",
                       "attempting to coerce")
        message(memo)
        x <- as.matrix(x)
    }
    
    if(nrow(x) < 1)
    {
        memo <- paste0("argument supplied to x needs at least one row")
        stop(memo)
    }

    # Check that rownames of x can be converted to integers
    if(is.null(rownames(x)))
    {
        memo <- paste0("all rownames of x are missing, they will be converted",
                       " to integers starting at 0")
        message(memo)
        rownames(x) = as.character(0:(nrow(x)-1))
    } else {
        naind <- which(is.na(as.integer(rownames(x))))
        if(length(naind)==nrow(x))
        {
            memo <- paste0("no rownames of x can be interpreted as integers, ",
                           "they will be converted to integers starting at 0")
            message(memo)
            rownames(x) = as.character(0:(nrow(x)-1))
        } else if(length(naind) > 0) {
            paste0("some rownames of x cannot be interpreted as integers, ",
                   "they will be removed")
            message(memo)
            x <- x[-naind,]
        }
    }

    return(list(x))
}
