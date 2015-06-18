#' align TvTi plots on y axis
#' 
#' align transition/transversion plots
#' @name align_y_TvTi
#' @param p1 main plot
#' @param p2 left subplot
#' @return ggplot object
#' @import gtable
#' @import gridExtra

align_y_TvTi <- function(p1, p2)
{
  # define the ggplot's as grobs and create a blank plot
  gA <- ggplotGrob(p1)
  gB <- ggplotGrob(p2)
  
  # Adjust the grob heights so p1, and p2 plots line up
  maxheight = grid::unit.pmax(gA$heights[2:5,], gB$heights[2:5,])
  gA$heights[2:5] <- as.list(maxheight)
  gB$heights[2:5] <- as.list(maxheight)
  
  # plot the grobs with grid.arrange
  p1 <- grid.arrange(gB, gA, ncol=2, nrow=1, widths=c(1,4))
  
  return(p1)
}