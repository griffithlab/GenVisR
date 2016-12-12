#' align TvTi plots on y axis
#'
#' align transition/transversion plots
#' @name TvTi_alignPlot
#' @param p1 main plot
#' @param p2 left expected value subplot
#' @param p3 bottom clinical subplot
#' @return ggplot object

TvTi_alignPlot <- function(p1=NULL, p2=NULL, p3=NULL)
{
    # define the ggplot's as grobs and create a blank plot
    gA <- ggplot2::ggplotGrob(p1)

    # Adjust the grob heights so p1, and p2 plots line up if p2 exists
    if(!is.null(p2))
    {
        # convert expected plot to grob
        gB <- ggplot2::ggplotGrob(p2)
        
        maxheight = grid::unit.pmax(gA$heights, gB$heights)
        gA$heights <- as.list(maxheight)
        gB$heights <- as.list(maxheight)
    }
    
    # adjust the grob widths so p1 and p3 line up if p3 exists
    if(!is.null(p3))
    {
        # Convert clinical plot to grob
        gC <- ggplot2::ggplotGrob(p3)
        gD <- grid::grid.rect(gp=grid::gpar(col="white"))
        
        maxwidth = grid::unit.pmax(gA$widths, gC$widths)
        gA$widths <- as.list(maxwidth)
        gC$widths <- as.list(maxwidth)        
    }

    # Build the final plot
    if(is.null(p3) & is.null(p2))
    {
        finalPlot <- gridExtra::arrangeGrob(gA, ncol=1, nrow=1)
    } else if(is.null(p3) & !is.null(p2)) {
        finalPlot <- gridExtra::arrangeGrob(gB, gA, ncol=2, nrow=1, widths=c(1,6))
    } else if(!is.null(p3) & is.null(p2)) {
        finalPlot <- gridExtra::arrangeGrob(gA, gC, ncol=1, nrow=2, heights=c(6,1))
    } else if(!is.null(p3) & !is.null(p2)) {
        finalPlot <- gridExtra::arrangeGrob(gB, gA, gD, gC, ncol=2, nrow=2, heights=c(6,1), widths=c(1,6))
    }

    return(finalPlot)
}
