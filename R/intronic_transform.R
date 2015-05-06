#' transform intronic space
#' 
#' transform intronic space
#' @name intronic_transform
#' @param x object of class data frame containing columns Type, end, start
#' @return Object of class data frame

intronic_transform <- function(x)
{
  x$width <- ifelse(x$Type == 'Intron', log(x$width), x$width)
  x <- x[order(x$start),]
  return(x)
}