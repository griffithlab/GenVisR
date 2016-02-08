#' extract ggplotGrob height
#' 
#' extract plot height of ggplotGrob object
#' @name extr_ggplotGrob_height
#' @param x ggplotGrob object
#' @return ggplotGrob height parameters

extr_ggplotGrob_height <- function(x)
{
  # Extract ggplotGrob height
  x <- x$heights[2:5,]
  return(x)
}