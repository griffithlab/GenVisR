#' align TvTi plots on y axis
#'
#' align transition/transversion plots
#' @name TvTi_alignPlot
#' @param p1 main plot
#' @param p2 left subplot
#' @return ggplot object

TvTi_alignPlot <- function(p1, p2)
{
    # define the ggplot's as grobs and create a blank plot
    gA <- ggplot2::ggplotGrob(p1)
    gB <- ggplot2::ggplotGrob(p2)

    # Adjust the grob heights so p1, and p2 plots line up
    maxheight = grid::unit.pmax(gA$heights[2:5,], gB$heights[2:5,])
    gA$heights[2:5] <- as.list(maxheight)
    gB$heights[2:5] <- as.list(maxheight)

    # plot the grobs with grid.arrange
    p1 <- gridExtra::arrangeGrob(gB, gA, ncol=2, nrow=1, widths=c(1,6))

    return(p1)
}
