#' plot mutation burden
#' 
#' plot a barchart showing mutation burden given by data frame
#' @name build_mutation_recurrence_b
#' @param x a data frame containing columns sample, mut_burden
#' @return a ggplot object

build_mutation_recurrence_b <- function(x)
{  
  # Define Theme
  theme <- theme(panel.border =  element_blank(), axis.line =  element_line(), panel.background=element_rect(fill='white'), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_blank(), axis.title.x=element_blank())
  
  # Define additional parameters
  y_label <- ylab("Mutation Burden")
  
  # ggplot2 call
  p1 <- ggplot(x, aes(x=sample, y=mut_burden)) + geom_bar(stat='identity', alpha=.75, width=1) + theme + y_label + legend
  
  return(p1)
}