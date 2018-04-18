#' assign ggplotGrob width
#'
#' assign width of ggplotGrob object
#' @name genCov_assign_ggplotGrob_width
#' @param x ggplotGrob object
#' @param max_width grob object specifying width to reassign ggplotGrob with
#' @return ggplotGrob
#' @noRd

genCov_assign_ggplotGrob_width <- function(x, max_width)
{
  # Assign a max width to ggplotGrob
    x$widths[2:5] <- as.list(max_width)
    return(x)
}
