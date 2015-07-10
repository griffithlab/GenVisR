#' plot mutation burden
#' 
#' plot a barchart showing mutations per MB
#' @name build.mutRecurA.mutSpec
#' @param data_frame a data frame in MAF format
#' @param coverage_space an integer specifying the coverage space in base pairs from which a mutation could occur
#' @param layers Additional ggplot2 layers to plot
#' @return a ggplot object

build.mutRecurA.mutSpec <- function(data_frame, coverage_space, layers=NULL)
{
  #############################################################################################################################
  #################################### Function to plot the top margin barplot ################################################
  #############################################################################################################################
  
  # Add in mutations per MB calculation
  data_frame$mutation_per_MB <- data_frame$mutation_total/coverage_space * 1000000
  
  #print(data_frame)
  
  # Alter GGplot2 Theme 
  theme <- theme(axis.ticks.x=element_blank(), axis.text.x=element_blank(), axis.title.x=element_blank(), legend.title=element_text(size=14))
  
  # Add Legend
  legend <- scale_fill_manual(name="Translational Effect", values=c("red", "blue"), breaks=c("Synonymous", "Non Synonymous"), drop=FALSE)
  
  # add y label
  y_label <- ylab('Mutations per MB')
  
  # additional parameters
  if(!is.null(layers))
  {
    layers <- layers
  } else {
    layers <- geom_blank()
  }
  
  # ggplot2 call
  p1 <- ggplot(data_frame, aes(x=sample, y=mutation_per_MB, fill=trv_type)) + geom_bar(stat='identity', alpha=.75, width=1) + theme_bw() + theme + y_label + legend + layers
  
  return(p1)
}