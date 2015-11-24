#' extract ggplotGrob height
#'
#' extract plot height of ggplotGrob object
#' @name genCov_extr_ggplotGrob_height
#' @param x ggplotGrob object
#' @return ggplotGrob height parameters

genCov_extr_ggplotGrob_height <- function(x)
{
  # Extract ggplotGrob height
    x <- x$heights[2:5,]
    return(x)
}
