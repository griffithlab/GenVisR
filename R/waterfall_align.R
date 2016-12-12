#' align plots
#'
#' align mutation landscape, mutation burden on sample, and mutation burden on
#' gene plots
#' @name waterfall_align
#' @param p1 ggplot object displaying a mutation landscape
#' @param p2 ggplot object displaying mutation burden on gene
#' @param p3 ggplot object displaying mutation burden on sample
#' @param p4 ggplot object displaying clinical information "optional"
#' @param p5 ggplot object displaying proportion of mutation types "optional"
#' @param section_heights Heights of each section (should sum to one)
#' @return a grob object
#' @importFrom gridExtra arrangeGrob

waterfall_align <- function(p2, p1, p3, p4, p5, 
    section_heights) {
    if (missing(section_heights)) {    
        if (missing(p4)) {
            if (is.null(p5)) {
                section_heights <- c(1, 4)
            } else {
                section_heights <- c(1, 4, 0.5)
            }
        } else {
            if (is.null(p5)) {
                section_heights <- c(1, 4, 0.5)
            } else {
                section_heights <- c(1, 4, 1.2, 0.5)
            }
        }
    }

    # define the ggplot's as grobs and create a blank plot
    gA <- suppressWarnings(ggplot2::ggplotGrob(p2))
    gB <- ggplot2::ggplotGrob(p1)
    gC <- ggplot2::ggplotGrob(p3)
    if (!is.null(p5)) {
        gE <- ggplot2::ggplotGrob(p5)
    }
    blankPanel <- grid::grid.rect(gp=grid::gpar(col="white"))
    if(!missing(p4))
    {
        gD <- ggplot2::ggplotGrob(p4)
    }

    # Adjust the grob widths so p1 and p3 plots line up
    if(!missing(p4))
    {
        if (!is.null(p5)) {
            maxwidth <- grid::unit.pmax(gB$widths,
                                       gC$widths,
                                       gD$widths,
                                       gE$widths)
            gE$widths <- as.list(maxwidth)
        } else {   
            maxwidth <- grid::unit.pmax(gB$widths,
                                       gC$widths,
                                       gD$widths
                                       )
        }
        gC$widths <- as.list(maxwidth)
        gB$widths <- as.list(maxwidth)
        gD$widths <- as.list(maxwidth)
    } else {
        if (!is.null(p5)) {
            maxwidth <- grid::unit.pmax(gB$widths,
                                       gC$widths,
                                       gE$widths)
            gE$widths <- as.list(maxwidth)
        } else {   
            maxwidth <- grid::unit.pmax(gB$widths,
                                       gC$widths
                                       )
        }
        gC$widths <- as.list(maxwidth)
        gB$widths <- as.list(maxwidth)
    }

    # Adjust the grob heights so p1, and p2 plots line up
    maxheight <- grid::unit.pmax(gA$heights, gB$heights)
    gA$heights <- as.list(maxheight)
    gB$heights <- as.list(maxheight)

    # plot the grobs with grid.arrange
    if(!missing(p4))
    {
        if (!is.null(p5)) {
            nrows <- 4
            grobs <- list(blankPanel, gC, gA, gB, blankPanel, gE, blankPanel, gD)
        } else {
            nrows <- 3
            grobs <- list(blankPanel, gC, gA, gB, blankPanel, gD)
        }
        p1 <- gridExtra::arrangeGrob(grobs = grobs,
                                     ncol=2, nrow=nrows, widths=c(.8, 4),
                                     heights=section_heights)
    } else {
        if (!is.null(p5)) {
            nrows <- 3
            grobs <- list(blankPanel, gC, gA, gB, blankPanel, gE)
        } else {
            grobs <- list(blankPanel, gC, gA, gB)
            nrows <- 2
        }
        p1 <- gridExtra::arrangeGrob(grobs = grobs,
                                     ncol=2, nrow=nrows,
                                     widths=c(0.8, 4), heights=section_heights)
    }

    return(p1)
}
