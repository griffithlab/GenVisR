#' assign ggplotGrob height
#'
#' assign height of ggplotGrob object
#' @name genCov_assign_ggplotGrob_height
#' @param x ggplotGrob object
#' @param max_height grob object specifying width to reassign ggplotGrob with
#' @return ggplotGrob

genCov_assign_ggplotGrob_height <- function(x, max_height)
{
  # Assign a max width to ggplotGrob
    x$heights[2:5] <- as.list(max_height)
    return(x)
}
