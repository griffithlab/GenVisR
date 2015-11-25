#' group coordinates
#'
#' given amino acid position coordinates on, before grouping and dodge on x axis
#' @name lolliplot_groupMutcoord
#' @param x numeric vector of position coordinates on x axis
#' @return numeric vector of grouped position coordinates on x axis

lolliplot_groupMutcoord <- function(x)
{
    # Sort the vector, Then group coordinates meeting pre-defined threshold
    x <- sort(x)
    for(i in 1:length(x))
    {
        # define the first element as belonging to group 1
        if(i == 1)
        {
            group <- 1
            begin <- x[i]
            group_vec <- c(group)
            next
        }

        # assign subsequent coordinates to the same group if threshold is met,
        # else start new group
        if(x[i] - x[i - 1] <= 50 & x[i] - begin <= 200)
        {
            group_vec <- c(group_vec, group)
        } else if(x[i] - x[i - 1] > 50 | x[i] - begin > 200) {
            begin <- x[i]
            group <- group + 1
            group_vec <- c(group_vec, group)
        }
    }

    return(group_vec)
}
