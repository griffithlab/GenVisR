#' Construct coverage cohort plot
#'
#' given a data frame construct a plot to display coverage as percentage bars for a group of samples
#' @name build.covBars
#' @param data_frame object of class data frame containing columns depth, sample, bp
#' @param col vector of colors for the coverage bars
#' @param plot_title character string for title of plot
#' @param background character string specifying backround color of plot
#' @param x_lab_size integer specifying the size of the X label
#' @param y_lab_size integer specifying the size of the Y label
#' @param facet_lab_size integer specifying the size of the faceted labels
#' @param layers Additional layers to be plotted, can be a theme but must be a ggplot layer
#' @return ggplot object
#' @import ggplot2

build.covBars <- function(data_frame, col, plot_title=NULL, background='grey90', x_lab_size=12, y_lab_size=12, facet_lab_size=10, layers=NULL)
{
  #xmelt <- na.omit(data_frame)
  #dummy_data <- data_frame
  xmelt <- data_frame

  # Define Theme of plot
  theme <- theme(strip.text.y=element_text(angle=0, size=facet_lab_size), axis.text.y=element_blank(), axis.ticks.y=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_rect(fill=background), legend.position='right', axis.title.x=element_text(size=x_lab_size, face='bold'), axis.title.y=element_text(size=y_lab_size, face='bold'))

  # Define parameters of plot
  facet <- facet_grid(sample ~ ., scales='free', space='free')
  fill_gradient <- scale_fill_gradientn(colours=col)
  ylabel <- ylab('Sample')
  xlabel <- xlab('Coverage')
  
  # Define main plot using boundaries and then plot actual data
  p1 <- ggplot(data=xmelt, mapping=aes_string(xmin='xmin', xmax='bp', ymin=0, ymax=1)) + geom_rect(alpha=0) + scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0))
  p1 <- p1 + geom_rect(data=xmelt, mapping=aes_string(xmin='xmin', xmax='bp', ymin=0, ymax=1, fill='depth'))
  
  # build the plot
  p1 <- p1 + fill_gradient + ylabel + xlabel + facet + theme
  
  # if there are other layers, add them
  if(!is.null(layers))
    p1 <- p1 + layers
  
  # if title is supplied plot it
  if(!is.null(plot_title))
    p1 <- p1 + ggtitle(plot_title)
  
  return(p1)
}