#' dodge coordinates
#'
#' given a data frame, dodge x coordinates ontop of each other
#' @name lolliplot_dodgeCoordY
#' @param x data frame containing columns coord_x_dodge
#' @param track character vector, one of "top", "bottom" specifying whether to
#' dodge in a positive or negative fashion
#' @return numeric vector of dodged position coordinates on y axis

lolliplot_dodgeCoordY <- function(x, track='top')
{
    for(i in 1:length(x$coord_x_dodge))
    {
        if(track == 'top' & i == 1)
        {
            pos <- .3
            orig_pos <- .3
            pos_change <- .1
        } else if(track == 'bottom' & i == 1) {
            pos <- -.3
            orig_pos <- -.3
            pos_change <- -.1
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
