#' transform intronic space
#' 
#' transform intronic space
#' @name intronic_transform
#' @param x object of class data frame containing columns Type, end, start
#' @return Object of class data frame

intronic_transform <- function(x)
{
  x$diff <- x$end - x$start + 1
  x$diff <- ifelse(x$Type == 'Intron', log(x$diff), x$diff)
  x <- x[order(x$start),]
  return(x)
}