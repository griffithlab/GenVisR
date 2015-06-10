#' dodge coordinates
#' 
#' given amino acid position coordinates on, before grouping and dodge on x axis
#' @name dodge_coord_x
#' @param x numeric vector of position coordinates on x axis
#' @return numeric vector of dodged position coordinates on x axis
#' @import FField

dodge_coord_x <- function(x, rep.fact=5000, rep.dist.lmt=500, attr.fact=.1, adj.max=.1, adj.lmt=.5, iter.max=50000)
{
  # Format into data frame with columns as x and y
  x <- as.data.frame(cbind(x, 0))
  colnames(x) <- c('x', 'y')
  
  # take the data frame and apply a repulsive force to coordiantes
  x <- FFieldPtRep(x, rep.fact=rep.fact, rep.dist.lmt=rep.dist.lmt, attr.fact=attr.fact, adj.max=adj.max, adj.lmt=adj.lmt, iter.max=iter.max)
  
  return(x$x)
}