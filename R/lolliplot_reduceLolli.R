#' Reduce Lolli
#' 
#' Reduce lollis stacked ontop of each other to the amount specified
#' @name lolliplot_reduceLolli
#' @param x Data frame with column name mutation_coord to reduce lollis on
#' @param max Integer specifying the maximum number of lollis to allow
#' @return Object of class data frame taking the reduced form of x
#' @importFrom plyr count

lolliplot_reduceLolli <- function(x, max=NULL)
{
    # if max is null no reduction is required
    if(is.null(max))
    {
        return(x)
    }
    
    # Get a frequency of counts from x
    coordFreq <- plyr::count(x$mutation_coord)
    colnames(coordFreq) <- c("x", "freq")
    
    keep <- vector('numeric')
    # Loop through mutations keeping only what is under max
    for(i in 1:nrow(coordFreq))
    {
        index <- which(x$mutation_coord == coordFreq$x[i])
        index <- index[1:max]
        keep <- c(keep, index)
    }
    
    # subset the input on what we want to keep
    x <- x[keep,]
    
    return(x)
}