#' construct CN plot
#' 
#' given a CN data frame plot points in ggplot
#' @name build_CN_plot
#' @param data_frame a data frame with columns Chr, Coord, Tumor, Normal, Diff, p_value
#' @param chr a character string specifying chromosome
#' @return ggplot2 object
#' @import ggplot2

build_CN_plot <- function(x, chr)
{
  # Define a THEME
  theme <- theme(axis.text.x=element_text(angle=30, hjust=1))
  
  # Define points to plot for the main plot and apply shading function
  cnpoints <- ggplot(data=x, mapping=aes(x=Coord, y=Diff)) + geom_point(data=x, mapping=aes(colour=Diff, alpha=1-p_value))
  shade_cn <- scale_color_gradient2(midpoint=0, low='#009ACD', mid='#646082', high='#C82536', space='Lab')
  
  # add additional labels and theme
  ylabel <- ylab('CNV Difference')
  xlabel <- xlab('Coordinate')
  
  # apply a smoothing function to the data and overlay on top of the raw data
  smooth_cn <- geom_smooth(fill='green', alpha=1, level=.999)
  
  # build the plot
  p1 <- cnpoints + shade_cn + smooth_cn + ylabel + xlabel + theme
  
  if(chr == 'all')
  {
    facet <- facet_wrap(~Chr, scales='free')
    p1 <- p1 + facet
  }
  
  return(p1)
}