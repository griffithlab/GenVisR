#' Construct coverage cohort plot
#'
#' given a matrix construct a plot to display coverage as percentage bars for a
#' group of samples
#' @name covBars_qual
#' @param x object of class matrix containing rows for the coverage and columns
#' the sample names
#' @param col vector of colors for the coverage bars
#' @importFrom grDevices rainbow
#' @return a list of data frame and color vector

covBars_qual <- function(x, col)
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

    # Check that col is the same length as the number of rows in x
    if(length(col)==0)
    {
        memo <- paste0("Argument supplied to col has 0 length... ",
                       "default colour scheme will be used")
        message(memo)
        col <- grDevices::rainbow(nrow(x),start=0,end=0.9)
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

    return(list(x, col))
}
