#' dodge coordinates
#' 
#' given amino acid position coordinates on, before grouping and dodge on y axis
#' @name dodge_coord_y
#' @param x data frame containing columns coord_x_dodge and group
#' @return numeric vector of dodged position coordinates on y axis

dodge_coord_y <- function(x, track='top')
{
  
  for(i in 1:length(x$group))
  {
    if(track == 'top')
    {
      pos <- .7
      orig_pos <- .7
      pos_change <- .1
    } else if(track == 'bottom') {
      pos <- -.7
      orig_pos <- -.7
      pos_change <- -.1
    }
    
    if(i == 1)
    {
      y_axis_vec <- c(pos)
      next
    } else {
      group_a <- x$group[i-1]
      x_coord_a <- x$coord_x_dodge[i-1]
      group_b <- x$group[i]
      x_coord_b <- x$coord_x_dodge[i]
    }
    
    if(group_a == group_b & x_coord_b == x_coord_a)
    {
      new_pos <- pos + pos_change
      y_axis_pos <- new_pos
      y_axis_vec <- c(y_axis_vec, y_axis_pos)
      pos <- new_pos
      next
    } else {
      pos <- orig_pos
    }
    
    if(group_a != group_b & x_coord_a <= x_coord_b)
    {
      new_pos <- pos
      y_axis_pos <- new_pos
      y_axis_vec <- c(y_axis_vec, y_axis_pos)
      pos <- new_pos
      next
    } else if(group_a == group_b) {
      y_axis_pos <- pos
      y_axis_vec <- c(y_axis_vec, y_axis_pos)
      next
    } else {
      pos <- orig_pos
      y_axis_pos <- pos
      y_axis_vec <- c(y_axis_vec, y_axis_pos)
      next
    }
  }
  
  return(y_axis_vec)
}