#' Construct coverage cohort plot
#'
#' given a matrix construct a plot to display coverage as percentage bars for a group of samples
#' @name covBars.qual
#' @param x object of class matrix containing rows for the coverage and columns the sample names
#' @param col vector of colors for the coverage bars
#' @return a list of data frame and color vector

covBars.qual <- function(x, col)
{
  # Check that x is a matrix with at least 1 row
  if(!is.matrix(x))
  {
    message("x is not a matrix, attempting to coerce")
    x <- as.matrix(x)
  }
  if(nrow(x) < 1)
  	stop("x needs at least one row")

  # Check that col is the same length as the number of rows in x
  if(length(col)==0)
  {
    message("col has a length of 0, default color scheme will be used")
    col <- rainbow(nrow(x),start=0,end=0.9)
  } #else if(length(col) < nrow(x))
  #{
  #  message("col is shorter than the number of rows in x, col will be repeated")
  #  col <- rep_len(col, nrow(x))
  #} else
  #{
  #  message("col is longer than the number of rows in x, col will be truncated")
  #  col <- col[1:nrow(x)]
  #}
  
  # Check that rownames of x can be converted to integers
  if(is.null(rownames(x))) 
  {
    message("all rownames of x are missing, they will be converted to integers starting at 0")
  	rownames(x) = as.character(0:(nrow(x)-1))
  } else
  {
    naind <- which(is.na(as.integer(rownames(x))))
    if(length(naind)==nrow(x)) {
      message("no rownames of x can be interpreted as integers, they will be converted to integers starting at 0")
      rownames(x) = as.character(0:(nrow(x)-1))
    } else if(length(naind) > 0)
    {
      message("some rownames of x cannot be interpreted as integers, they will be removed")
      x <- x[-naind,]
      #col <- col[-naind]
    }
  }
  
  return(list(x, col))
}