#' plot mutation burden
#' 
#' plot a barchart showing mutations per MB
#' @name plot_mutation_recurrence
#' @param data_frame a data frame in MAF format
#' @param coverage_space an integer specifying the coverage space in base pairs from which a mutation could occur
#' @return a ggplot object

plot_mutation_recurrence <- function(data_frame, coverage_space)
{
  #############################################################################################################################
  #################################### Function to plot the top margin barplot ################################################
  #############################################################################################################################
  
  # Add in mutations per MB calculation
  data_frame$mutation_per_MB <- data_frame$mutation_total/coverage_space * 1000000
  
  #print(data_frame)
  
  # Alter GGplot2 Theme 
  theme <- theme(panel.border = element_blank(), axis.line = element_line(), panel.background=element_rect(fill='white'), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_blank(), axis.title.x=element_blank(), legend.title=element_text(size=14), axis.title.y=element_text(size=16))   
  
  # Add Legend
  legend <- scale_fill_manual(name="Translational Effect", values=c("red", "blue"), breaks=c("Synonymous", "Non Synonymous"), drop=FALSE)
  
  # add y label
  y_label <- ylab('Mutations per MB')
  
  # ggplot2 call
  p1 <- ggplot(data_frame, aes(x=sample, y=mutation_per_MB, fill=trv_type)) + geom_bar(stat='identity', alpha=.75, width=1) + theme_bw() + theme + y_label + legend
  
  return(p1)
}