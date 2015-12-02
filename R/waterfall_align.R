#' align plots
#'
#' align mutation landscape, mutation burden on sample, and mutation burden on
#' gene plots
#' @name waterfall_align
#' @param p1 ggplot object displaying a mutation landscape
#' @param p2 ggplot object displaying mutation burden on gene
#' @param p3 ggplot object displaying mutation burden on sample
#' @param p4 ggplot object displaying clinical information "optional"
#' @return a grob object
#' @importFrom gridExtra arrangeGrob

waterfall_align <- function(p2, p1, p3, p4)
{
    # define the ggplot's as grobs and create a blank plot
    gA <- suppressWarnings(ggplot2::ggplotGrob(p2))
    gB <- ggplot2::ggplotGrob(p1)
    gC <- ggplot2::ggplotGrob(p3)
    blankPanel <- grid::grid.rect(gp=grid::gpar(col="white"))
    if(!missing(p4))
    {
        gD <- ggplot2::ggplotGrob(p4)
    }

    # Adjust the grob widths so p1 and p3 plots line up
    if(!missing(p4))
    {
        maxwidth = grid::unit.pmax(gB$widths[2:5,],
                                   gC$widths[2:5,],
                                   gD$widths[2:5,])
        gC$widths[2:5] <- as.list(maxwidth)
        gB$widths[2:5] <- as.list(maxwidth)
        gD$widths[2:5] <- as.list(maxwidth)
    } else {
        maxwidth = grid::unit.pmax(gB$widths[2:5,], gC$widths[2:5,])
        gC$widths[2:5] <- as.list(maxwidth)
        gB$widths[2:5] <- as.list(maxwidth)
    }

    # Adjust the grob heights so p1, and p2 plots line up
    maxheight = grid::unit.pmax(gA$heights[2:5,], gB$heights[2:5,])
    gA$heights[2:5] <- as.list(maxheight)
    gB$heights[2:5] <- as.list(maxheight)

    # plot the grobs with grid.arrange
    if(!missing(p4))
    {
        p1 <- gridExtra::arrangeGrob(blankPanel, gC, gA, gB, blankPanel, gD,
                                     ncol=2, nrow=3, widths=c(.8,4),
                                     heights=c(1,4,1.2))
    } else {
        p1 <- gridExtra::arrangeGrob(blankPanel, gC, gA, gB, ncol=2, nrow=2,
                                     widths=c(1,4), heights=c(1,4))
    }

    return(p1)
}
