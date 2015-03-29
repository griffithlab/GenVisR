#' correct spread coordinates
#' 
#' correct skewing of coord points introduced by spread_coord_points for 
#' @name correct_spread_coord_points
#' @param x list 
#' @return numeric vector of spread points

correct_spread_coord_points <- function(x)
{
  
  # if frequency of group == 1 return the original coordinate else return the new spread coordinate
  if(x[4] == 1)
  {
    return(x[2])
  } else {
    return(x[5])
  }
}