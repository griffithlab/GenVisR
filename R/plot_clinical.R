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
  theme <- theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(), panel.background=element_rect(fill='white', colour='white'), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.text.x=element_blank(), axis.title.x=element_text(size=16), legend.title=element_blank(), axis.title.y=element_blank(), axis.text.y=element_text(size=14, colour='black'), legend.position='right', plot.margin=unit(c(-1,1,.5,.7), 'lines'))
  
  # Define various parameters of plot
  fill <- scale_fill_manual(breaks=c('any relapse', 'no relapse', 'death from breast cancer (event)', 'alive (uncensored)', 'lobular', 'ductal', 'Her2', 'Basal', 'LumA', 'LumB', 'no relapse', 'Normal', 'unknown'), values=c('alive (uncensored)'='#225533', 'any relapse'='#44bbcc', 'lobular'='#CE1836', 'death from breast cancer (event)'='#F85931', 'ductal'='#C7FCD7', 'Her2'='hotpink2', 'Basal'='firebrick2', 'LumA'='blue4', 'LumB'='deepskyblue', 'no relapse'='#4169E1', 'Normal'='green4', 'unknown'='#000000'))
  leg_guide <- guides(fill=guide_legend(ncol=2))
  
  # Define the main plot
  p1 <- ggplot(x, aes(x=sample, y=variable, fill=value)) + geom_tile() + theme + fill + leg_guide + xlab('Sample n=625')
  
  return(p1)
}
