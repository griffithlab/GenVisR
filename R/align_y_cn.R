#' align cn plots on x axis
#' 
#' given a chromosome and cn plot align plot widths
#' @name align_y_cn
#' @param p1 ggplot object of chromosome
#' @param p2 ggplot object of CN
#' @return ggplot object

align_y_cn <- function(p1, p2)
{
  
  require(gridExtra)
  require(gtable)
  
  # define the ggplot's as grobs and create a blank plot
  gA <- ggplotGrob(p1)
  gB <- ggplotGrob(p2)
  
  # Adjust the grob heights so p1, and p2 plots line up
  maxwidth = grid::unit.pmax(gA$widths[2:5,], gB$widths[2:5,])
  gA$widths[2:5] <- as.list(maxwidth)
  gB$widths[2:5] <- as.list(maxwidth)
  
  # plot the grobs with grid.arrange
  p1 <- grid.arrange(gA, gB, ncol=1, nrow=2, heights=c(1,2))
  
  return(p1)
}