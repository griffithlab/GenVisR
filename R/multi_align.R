#' align CN/LOH plots on x axis
#'
#' given a chromosome and CN/LOH plot align plot widths
#' @name multi_align
#' @param p1 ggplot object of chromosome
#' @param p2 ggplot object of CN or LOH
#' @return ggplot object
#' @importFrom grid unit.pmax
#' @importFrom gridExtra arrangeGrob

multi_align <- function(p1, p2)
{
    # define the ggplot's as grobs and create a blank plot
    gA <- ggplot2::ggplotGrob(p1)
    gB <- ggplot2::ggplotGrob(p2)

    # Adjust the grob heights so p1, and p2 plots line up
    maxwidth = grid::unit.pmax(gA$widths, gB$widths)
    gA$widths <- as.list(maxwidth)
    gB$widths <- as.list(maxwidth)

    # plot the grobs with grid.arrange
    p1 <- gridExtra::arrangeGrob(gA, gB, ncol=1, nrow=2, heights=c(1,2))

    return(p1)
}
