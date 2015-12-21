#' build coverage plot
#'
#' given data build a coverage plot to represent the data
#' @name genCov_buildCov
#' @param data_frame an object of class data frame containing columns
#' stop and cov
#' @param x_limits vector giving x-axis limits for plot, inferred from data
#' if not specified
#' @param colour character string specifying the color of the data in the plot
#' @param plot_type character string specifying one of line, area, bar, or point
#' for data display
#' @param display_x_axis boolean specifying whether to plot x-axis labels
#' @param layers additional ggplot2 layers to plot
#' @return ggplot object
#' @import ggplot2

genCov_buildCov <- function(data_frame, x_limits=NULL, display_x_axis=TRUE,
                            colour="blue", plot_type="line", layers=NULL)
{
     # Specify various parameters of the plot
    line <- geom_line(colour=colour)
    area <- geom_area(colour=colour, fill=colour)
    bar <- geom_bar(fill=colour, width=1, stat="identity")
    point <- geom_point(fill=colour, colour=colour)

    if(is.null(x_limits))
    {
        x_limits <- xlim(c(min(data_frame$end), max(data_frame$end)))
    } else {
        x_limits <- xlim(x_limits)
    }

    if(!is.null(layers))
    {
        layers <- layers
    } else {
        layers <- geom_blank()
    }

    # Define the theme
    if(display_x_axis == TRUE)
    {
        theme <- theme(axis.title.x=element_blank())
    } else {
        theme <- theme(axis.title.x=element_blank(),
                       axis.text.x=element_blank(),
                       axis.ticks.x=element_blank())
    }
    
    # Define y label
    y_label <- ylab("Coverage Depth")
    
    # Define the main plot
    cov_plot <- ggplot(data_frame, aes_string(x='end', y='cov')) + x_limits +
        theme_bw() + theme + y_label + layers

    # Define Control structure for plot type
    if(toupper(plot_type) == "LINE")
    {
        cov_plot <- cov_plot + line
    } else if(toupper(plot_type) == "AREA") {
        cov_plot <- cov_plot + area
    } else if(toupper(plot_type) == "BAR"){
        cov_plot <- cov_plot + bar
    } else if(toupper(plot_type) == "POINT") {
        cov_plot <- cov_plot + point
    } else {
        memo <- paste0("Do not recoginze: ", plot_type,
                       ", please specify 'line', 'area' or 'bar'")
        stop(memo)
    }

    return(cov_plot)
}
