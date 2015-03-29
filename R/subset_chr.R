#' subset based on chr
#' 
#' given a data frame with cytogenetic band locations subset out specific chromosome
#' @name subset_chr
#' @param data_frame a data frame with columns Chr, Coord, Tumor, Normal, Diff, p_value
#' @param chr character string specifying UCSC chromosome to plot one of chr... or all
#' @return object of class data frame

subset_chr <- function(x, chr)
{
  # subset data frame based on value specified in chr argument
  if(chr == 'all')
  {
    return(x)
  }
  else if(any(chr == levels(x$Chr)))
  {
    x <- x[x$Chr == chr,]
    return(x)
  } else {
    stop("chr string does not match chr in data.frame")
  }
}