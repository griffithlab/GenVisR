#' plot clinical information
#' 
#' given a data frame with columns names sample, variable, and value create a ggplot2 object
#' @name plot_clinical
#' @param x a data frame in "long" format giving additional information to be plotted, requires columns "sample", "variable", and "value"
#' @return a grob object
#' @import ggplot2

plot_clinical <- function(x)
{
  # Define the theme
  theme <- theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(), panel.background=element_rect(fill='white', colour='white'), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.text.x=element_blank(), axis.title.x=element_text(size=16), legend.title=element_blank(), axis.title.y=element_blank(), axis.text.y=element_text(size=14, colour='black'), legend.position='right')
  
  
  # Define the main plot
  p1 <- ggplot(x, aes(x=sample, y=variable, fill=value)) + geom_tile() + theme
  
  return(p1)
}
