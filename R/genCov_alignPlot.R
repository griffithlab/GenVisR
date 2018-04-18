#' align plots on an axis
#'
#' given a list of plots, align them on plotting space
#' @name genCov_alignPlot
#' @param plot_list list of ggplot objects
#' @param axis character string to specify the axis to align plotting space on,
#' one of both, width, height
#' @return ggplotGrob object
#' @importFrom gridExtra arrangeGrob
#' @noRd

genCov_alignPlot <- function(plot_list, axis='both')
{
    # convert all ggplot objects to grob obects
    plots <- lapply(plot_list, ggplot2::ggplotGrob)

    # if specified align the plot widths
    if(axis == 'width' | axis == 'both')
    {
        # Obtain the max width in all plots in the list
        maxwidth <- do.call(grid::unit.pmax,
                            lapply(plots, genCov_extr_ggplotGrob_width))

        # set the max width for all plots in the list
        plots <- lapply(plots, genCov_assign_ggplotGrob_width, maxwidth)
    }

    # if specified alter the plot heights
    if(axis == 'height' | axis == 'both')
    {
        # Obtain the max height in all plots in the list
        maxheight <- do.call(grid::unit.pmax,
                             lapply(plots, genCov_extr_ggplotGrob_height))

        # set the max height for all plots in the list
        plots <- lapply(plots, genCov_assign_ggplotGrob_width, maxheight)
    }

    # Combine the plots of now equal widths/heights
    plots <- do.call(gridExtra::arrangeGrob, c(plots, ncol=1, nrow=length(plots)))

    return(plots)
}
