#' Construct CN frequency plot
#'
#' given a data frame construct a plot to display proportions of losses and gains across the genome
#' @name build.cnFreq
#' @param data_frame object of class data frame containing columns chromosome, start, end, gain, and loss
#' @param plotType character string to determine whether to plot proportions or frequencies
#' @param plot_title character string for title of plot
#' @param background character string specifying backround color of plot
#' @param CN_low_colour character string specifying low value of colour gradient 
#' @param CN_high_colour character string specifying high value of colour gradient
#' @param x_lab_size integer specifying the size of the X label
#' @param y_lab_size integer specifying the size of the Y label
#' @param facet_lab_size integer specifying the size of the faceted labels
#' @param layers Additional layers to be plotted, can be a theme but must be a ggplot layer
#' @return ggplot object
#' @import ggplot2

build.cnFreq <- function(data_frame, plotType, plot_title=NULL, background='grey90', CN_low_colour='#002EB8', CN_high_colour='#A30000', x_lab_size=12, y_lab_size=12, facet_lab_size=10, layers=NULL)
{
  #x <- na.omit(data_frame) this causes problems converting factors to numbers
  x <- data_frame

  # Define Theme of plot
  #theme <- theme(strip.text.y=element_text(angle=0, size=facet_lab_size), strip.text.x=element_text(size=facet_lab_size), axis.text.y=element_blank(), axis.ticks.y=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), panel.background=element_rect(fill=background), legend.position='right', axis.title.x=element_text(size=x_lab_size, face='bold'), axis.title.y=element_text(size=y_lab_size, face='bold'))
  theme <- theme(strip.text.x=element_text(size=facet_lab_size), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), panel.background=element_rect(fill=background), legend.position='right', axis.title.x=element_text(size=x_lab_size, face='bold'), axis.title.y=element_text(size=y_lab_size, face='bold'))

  # Define parameters of plot
  facet <- facet_grid(. ~ chromosome, scales='free', space='free')
  ylabel <- ylab(ifelse(plotType=="prop","Proportion of Copy Number Gains/Losses","Frequency of Copy Number Gains/Losses"))
  xlabel <- xlab('Chromosomes')
  
  # Define main plot using boundaries in dummy data and then plot actual data
  #ymin <- ifelse(plotType=="prop", 1, min(as.numeric(as.character(x$loss)),na.rm=T))
  ymax <- ifelse(plotType=="prop", 1, max(as.numeric(as.character(x$obs)),na.rm=T))
  p1 <- ggplot(data=x, mapping=aes_string(xmin='start', xmax='end', ymin=-1*ymax, ymax=ymax)) + geom_rect(alpha=0) + scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0))
  p1 <- p1 + geom_rect(data=x, mapping=aes_string(xmin='start', xmax='end', ymin='loss', ymax=0), fill=CN_low_colour)
  p1 <- p1 + geom_rect(data=x, mapping=aes_string(xmin='start', xmax='end', ymin=0, ymax='gain'), fill=CN_high_colour)
  p1 <- p1 + geom_hline(aes(yintercept=0), linetype="dotted")
  
  # build the plot
  p1 <- p1 + ylabel + xlabel + facet + theme
  
  # if there are other layers, add them
  if(!is.null(layers))
    p1 <- p1 + layers
  
  # if title is supplied plot it
  if(!is.null(plot_title))
    p1 <- p1 + ggtitle(plot_title)
  
  return(p1)
}