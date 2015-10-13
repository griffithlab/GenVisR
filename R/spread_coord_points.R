#' spread coordinates
#'
#' given a data frame with median and group spreadd coordinates within that group (I think this function is deprecated)
#' @name spread_coord_points
#' @param x list
#' @param spread_degree degree of repulsion between points
#' @return numeric vector of spread points

spread_coord_points <- function(x, spread_degree)
{
  # spread mutation points based on the count of points in a group
    spread_points <- seq(from=x[3] - (2.5*spread_degree), to=x[3] + (2.5*spread_degree), length.out=x[4])
    return(spread_points)
}
