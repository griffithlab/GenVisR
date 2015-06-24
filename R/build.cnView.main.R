#' construct CN plot
#' 
#' given a CN data frame plot points in ggplot
#' @name build.cnView.main
#' @param data_frame a data frame with columns chromosome, coordinate, cn, p_value
#' @param chr a character string specifying chromosome
#' @param cnDiff Boolean specifying whether values in cn are copy number differences or actual copy number
#' @return ggplot2 object
#' @import ggplot2

build.cnView.main <- function(x, chr, cnDiff=FALSE)
{
  # Define various parameters of the plot
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
  
  # Define points to plot for the main plot and apply shading function
  cnpoints <- geom_point(data=x, mapping=aes(x=coordinate, y=cn, colour=cn, alpha=1-p_value))

  
  # build the plot
  p1 <- ggplot() + cnpoints + shade_cn + ylabel + xlabel + theme
  
  if(chr == 'all')
  {
    facet <- facet_wrap(~chromosome, scales='free')
    p1 <- p1 + facet
  }
  
  return(p1)
}