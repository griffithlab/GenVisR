#' plot clinical information
#' 
#' given a data frame with columns names sample, variable, and value create a ggplot2 object
#' @name plot_clinical
#' @param x a data frame in "long" format giving additional information to be plotted, requires columns "sample", "variable", and "value"
#' @param clin.legend.col an integer specifying the number of columns to plot in the legend
#' @param clin.var.colour a named character vector specifying the mapping between colors and variables
#' @param clin.var.order a character vector of variables to order the legend by
#' @return a grob object
#' @import ggplot2

plot_clinical <- function(x, clin.legend.col=1, clin.var.colour=NULL, clin.var.order=NULL)
{
  # Define parameters
  x_label <- xlab(paste0("Sample n=", length(unique(x$sample))))
  leg_guide <- guides(fill=guide_legend(ncol=clin.legend.col))
  
  if(!is.null(clin.var.colour) & is.null(clin.var.order))
  {
    clin_fill_colour <- scale_fill_manual(values=clin.var.colour)
  } else if(!is.null(clin.var.colour) & !is.null(clin.var.order))
  {
    clin_fill_colour <- scale_fill_manual(breaks=clin.var.order, values=clin.var.colour)
  } else if(is.null(clin.var.colour) & !is.null(clin.var.order))
  {
    clin_fill_colour <- scale_fill_manual(breaks=clin.var.order)
  }
  
  # Define the theme
  theme <- theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(), panel.background=element_rect(fill='white', colour='white'), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.text.x=element_blank(), axis.title.x=element_text(size=16), legend.title=element_blank(), axis.title.y=element_blank(), axis.text.y=element_text(size=14, colour='black'), legend.position='right', legend.title=element_text(size=14), legend.text=element_text(size=14))   
  
  # Define the main plot
  if(!is.null(clin.var.colour))
  {
    p1 <- ggplot(x, aes(x=sample, y=variable, fill=value)) + geom_tile() + theme + x_label + leg_guide + clin_fill_colour
  } else {
    p1 <- ggplot(x, aes(x=sample, y=variable, fill=value)) + geom_tile() + theme + x_label + leg_guide
  }
  
  
  return(p1)
}
