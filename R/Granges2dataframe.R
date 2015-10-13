#' Convert Granges object to dataframe
#'
#' Convert a Granges object with meta data GC content to a object of class data frame
#' @name Granges2dataframe
#' @param gr A Granges object to convert to data frame
#' @return Object of class data frame

Granges2dataframe <- function(gr)
{
  # extract the relevant data and bind together as a data frame
    range <- as.data.frame(ranges(gr))
    meta <- as.data.frame(mcols(gr))
    genomic_data <- cbind(range, meta)

  # subset and return the desired output
    genomic_data <- genomic_data[,c('start', 'end', 'GC', 'width', 'txname')]

    return(genomic_data)
}
