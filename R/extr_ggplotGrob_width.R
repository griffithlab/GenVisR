#' extract ggplotGrob width
#'
#' extract plot width of ggplotGrob object
#' @name extr_ggplotGrob_width
#' @param x ggplotGrob object
#' @return ggplotGrob width parameters

extr_ggplotGrob_width <- function(x)
{
  # Extract ggplotGrob width
    x <- x$widths[2:5,]
    return(x)
}
