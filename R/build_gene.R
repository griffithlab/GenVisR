#' build gene plot
#' 
#' given a data frame with gene feature information build the ggplot2 object
#' @name build_gene
#' @param data_frame an object of class data frame specifying gene feature information
#' @return ggplot object

build_gene <- function(data_frame)
{
  require(ggplot2)
  require(grid)
  
  # Define various parameters of plot
  gene_features <- geom_rect(data=data_frame, mapping=aes(xmin=start, xmax=end, ymin=Upper, ymax=Lower, fill=GC))
  gene_track <- geom_segment(data=data_frame, mapping=aes(x=segStart, xend=segEnd, y=Mid, yend=Mid), arrow = arrow(length=unit(0.1,"cm")))
  
  # Define the theme of the plot
  theme <- theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position='top')
  
  # Define the main plot
  gene_plot <- ggplot() + gene_track + gene_features + theme
  
  return(gene_plot)
}