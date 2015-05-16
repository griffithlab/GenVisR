#' build coverage plot
#' 
#' given data build a coverage plot to represent the data
#' @name build_coverage
#' @param data_frame an object of class data frame containing columns stop and cov
#' @param xlimits vector giving x-axis limits for plot, inferred from data if not specified
#' @param colour character string specifying the color of the data in the plot
#' @param plot_type character string specifying one of line, area for data display
#' @param display_x_axis boolean specifying whether to plot x-axis labels
#' @return ggplot object
#' @import ggplot2

build_coverage <- function(data_frame, x_limits=NULL, display_x_axis=TRUE, colour="blue", plot_type="line")
{ 
  # Specify various parameters of the plot
  line <- geom_line(colour=colour)
  area <- geom_area(colour=colour)
  
  if(is.null(x_limits))
  {
    x_limits <- xlim(c(min(data_frame$end), max(data_frame$end)))
  } else {
    x_limits <- xlim(x_limits)
  }
  
  # Define the theme
  if(display_x_axis == TRUE)
  {
    theme <- theme(axis.title.x=element_blank(), axis.title.y=element_blank())
  } else {
    theme <- theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank())
  }
  
  # Define the main plot
  cov_plot <- ggplot(data_frame, aes(x=end, y=cov)) + x_limits + theme
  
  # Define Control structure for plot type
  if(plot_type == "line")
  {
    cov_plot <- cov_plot + line
  } else if(plot_type == "area") {
    cov_plot <- cov_plot + area
  } else {
    message <- paste0("Do not recoginze: ", plot_type, ", please specify 'line' or 'area'")
    stop(message)
  }
  
  return(cov_plot)
}