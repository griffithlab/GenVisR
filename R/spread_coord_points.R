#' spread coordinates
#' 
#' given a data frame with median and group spreadd coordinates within that group
#' @name spread_coord_points
#' @param x list 
#' @return numeric vector of spread points

spread_coord_points <- function(x)
{
  # spread mutation points based on the count of points in a group
  spread_points <- seq(from=x[3] - (2.5*15), to=x[3] + (2.5*15), length.out=x[4])
  return(spread_points)
}