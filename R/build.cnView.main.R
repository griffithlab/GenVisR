#' construct CN plot
#' 
#' given a CN data frame plot points in ggplot
#' @name build.cnView.main
#' @param x a data frame with columns chromosome, coordinate, cn, p_value
#' @param y a data frame with columns chromosome, coordinate for plotting boundaries
#' @param z a data frame with columns chromsome, start, end, segmean specifying segments called from copy number (optional)
#' @param chr a character string specifying chromosome
#' @param cnDiff Boolean specifying whether values in cn are copy number differences or actual copy number
#' @param layers additional ggplot2 layers to add
#' @return ggplot2 object
#' @import ggplot2

build.cnView.main <- function(x, y, z=NULL, chr, cnDiff=FALSE, layers=NULL)
{
  # Define various parameters of the plot
  dummy_data <- geom_point(data=y, mapping=aes_string(x='coordinate', y=2), alpha=0)
  
  theme <- theme(axis.text.x=element_text(angle=30, hjust=1))
  if(cnDiff == TRUE)
  {
    # cn fill colors
    shade_cn <- scale_color_gradient2(midpoint=0, low='#009ACD', mid='#646082', high='#C82536', space='Lab')
    # y label
    ylabel <- ylab('CN Difference')
  } else if(cnDiff == FALSE)
  {
    # cn fill colors
    shade_cn <- scale_color_gradient2(midpoint=2, low='#009ACD', mid='#646082', high='#C82536', space='Lab')
    # y label
    ylabel <- ylab('Copy Number')
  }
  xlabel <- xlab('Coordinate')
  
  if(!is.null(layers))
  {
    layers <- layers
  } else {
    layers <- geom_blank()
  }
  
  # if x contains a p_value column set an alpha for it and plot points
  if(any('p_value' %in% colnames(x)))
  {
    x$transparency <- 1-x$p_value
    cnpoints <- geom_point(data=x, mapping=aes_string(x='coordinate', y='cn', colour='cn', alpha='transparency'))
  } else {
    cnpoints <- geom_point(data=x, mapping=aes_string(x='coordinate', y='cn', colour='cn'))
  }
  
  # Define segments for main plot
  if(!is.null(z))
  {
    cnseg <- geom_segment(data=z, mapping=aes_string(x='start', xend='end', y='segmean', yend='segmean'), colour='green', size=2)
  } else {
    cnseg <- geom_blank()
  }
  
  # build the plot
  p1 <- ggplot() + cnpoints + shade_cn + ylabel + xlabel + theme_bw() + theme + cnseg + dummy_data + layers
  
  if(chr == 'all')
  {
    facet <- facet_wrap(~chromosome, scales='free')
    p1 <- p1 + facet
  }
  
  return(p1)
}