#' build coverage plot
#' 
#' given data build a coverage plot to represent the data
#' @name build_coverage
#' @param data_frame an object of class data frame containing columns stop and cov
#' @param colour character string specifying the color of the data in the plot
#' @param plot_type character string specifying one of line, area for data display
#' @return ggplot object

build_coverage <- function(data_frame, colour="blue", plot_type="line")
{
  require(ggplot2)
  
  # Specify various parameters of the plot
  line <- geom_line(colour=colour)
  area <- geom_area(colour=colour)
  
  # Define the theme
  theme <- theme(axis.title.x=element_blank(), axis.title.y=element_blank())
  
  # Define the main plot
  cov_plot <- ggplot(data_frame, aes(x=end, y=cov)) + theme
  
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