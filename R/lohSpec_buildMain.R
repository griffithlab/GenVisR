#' Plot LOH data
#'
#' Build a ggplot2 object displaying calculated LOH data
#' @name lohSpec_buildMain
#' @param x object of class dataframe with loh difference
#' calculations and column names "window_start", "window_stop", "chromosome",
#' "sample", and "loh_diff"
#' @param dummyData object of class dataframe with column names "chromosome",
#' "start", "end" specifying chromosome boundaries
#' @param colourScheme Character vector specifying the colour scale to use from
#' the viridis package. One of "viridis", "magma", "plasma", or "inferno".
#' @param plotLayer Valid ggplot2 layer to be added to the plot.
#' for the legend's parameters
#' @return object of class ggplot2
#' @import ggplot2
#' @importFrom viridis scale_fill_viridis

lohSpec_buildMain <- function(x, dummyData, colourScheme="inferno",
                              plotLayer=NULL)
{
    # define dummy data which will be chromosome boundaries, these are plotted
    # but are transparent and will not appear in the plot
    dummy_data <- geom_rect(data=dummyData, aes_string(xmin='start', xmax='end',
                                                       ymin=-1, ymax=1),alpha=0)
    # Define the main plot
    data <- geom_rect(data=x, aes_string(xmin='window_start',
                                         xmax='window_stop',
                                         ymin=-1,
                                         ymax=1, fill='loh_diff_avg'))


    # Define additional plot parameters
    facet <- facet_grid(sample ~ chromosome, scales="free", space="free")

    x_scale <- scale_x_continuous(expand = c(0, 0))
    y_scale <- scale_y_continuous(expand = c(0,0))

    lab_x <- xlab("Chromosome")
    lab_y <- ylab("Sample")

    # Define plot aesthetics
    BWscheme <- theme_bw()
    plotTheme <- theme(axis.ticks.x=element_blank(),
                       axis.text.x=element_blank(),
                       axis.ticks.y=element_blank(),
                       axis.text.y=element_blank(),
                       panel.grid.major=element_blank(),
                       panel.grid.minor=element_blank())
    
    # plot an additional layer if specified
    if(!is.null(plotLayer))
    {
        plotLayer <- plotLayer
    } else {
        plotLayer <- geom_blank()
    }
    
    # Define colour rame
    LOHgradient <- viridis::scale_fill_viridis("Avg. VAF Difference",
                                               option=colourScheme)

    # Build the plot
    tmp <- data.frame(x=0, y=0)
    p1 <- ggplot(data=tmp, aes(y=0)) + dummy_data + data + facet + x_scale + y_scale + 
        lab_x + lab_y + BWscheme + LOHgradient + plotTheme + plotLayer
     
    return(p1)
}
