#' build gene plot
#' 
#' given a data frame with gene feature information build the ggplot2 object
#' @name build_gene
#' @param data_frame an object of class data frame specifying gene feature information
#' @param display_axis Boolean specifying whether to display X axis coordinate values
#' @return ggplot object

build_gene <- function(data_frame, display_x_axis=T, x_limits=NULL)
{
  require(ggplot2)
  require(grid)
  
  # Define various parameters of plot
  gene_features <- geom_rect(data=data_frame, mapping=aes(xmin=start, xmax=end, ymin=Upper, ymax=Lower, fill=GC))
  gene_track <- geom_segment(data=data_frame, mapping=aes(x=segStart, xend=segEnd, y=Mid, yend=Mid))
  if(is.null(x_limits))
  {
    xlimits <- xlim(c(min(data_frame$start), max(data_frame$end)))
  } else {
    xlimits <- xlim(x_limits)
  }
  
  # Define the theme of the plot
  if(display_x_axis == TRUE)
  {
    theme <- theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position='top')
  } else {
    theme <- theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position='top')
  }
  
  # Define the main plot
  gene_plot <- ggplot() + gene_track + gene_features + theme + xlimits
  
  return(gene_plot)
}