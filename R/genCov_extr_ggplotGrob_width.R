#' extract ggplotGrob width
#'
#' extract plot width of ggplotGrob object
#' @name genCov_extr_ggplotGrob_width
#' @param x ggplotGrob object
#' @return ggplotGrob width parameters

genCov_extr_ggplotGrob_width <- function(x)
{
  # Extract ggplotGrob width
    x <- x$widths[2:5,]
    return(x)
}
