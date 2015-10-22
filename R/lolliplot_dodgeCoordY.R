#' dodge coordinates
#'
#' given amino acid position coordinates on, before grouping and dodge on y axis
#' @name lolliplot_dodgeCoordY
#' @param x data frame containing columns coord_x_dodge and group
#' @param track character vector, one of "top", "bottom" specifying whether to
#' dodge in a positive or negative fashion
#' @return numeric vector of dodged position coordinates on y axis

lolliplot_dodgeCoordY <- function(x, track='top')
{
    for(i in 1:length(x$coord_x_dodge))
    {
        if(track == 'top' & i == 1)
        {
            pos <- 2
            orig_pos <- 2
            pos_change <- .5
        } else if(track == 'bottom' & i == 1) {
            pos <- -2
            orig_pos <- -2
            pos_change <- -.5
        }

        if(i == 1)
        {
            y_axis_vec <- c(pos)
            next
        } else {
            x_coord_a <- x$coord_x_dodge[i-1]
            x_coord_b <- x$coord_x_dodge[i]
        }

        if(x_coord_b == x_coord_a)
        {
            new_pos <- pos + pos_change
            y_axis_pos <- new_pos
            y_axis_vec <- c(y_axis_vec, y_axis_pos)
            pos <- new_pos
            next
        } else {
            pos <- orig_pos
            y_axis_vec <- c(y_axis_vec, pos)
        }
    }

    return(y_axis_vec)
}
