#' plot mutation burden
#' 
#' plot a barchart showing mutation burden given by data frame
#' @name build_mutation_recurrence_b
#' @param x a data frame containing columns sample, mut_burden
#' @return a ggplot object
#' @import ggplot2

build_mutation_recurrence_b <- function(x)
{  
  # add in fake column for legend (necessary to have legend for proper plot alignment)
  # make everything white to hide legend
  x$Type <- "Non Synonymous"
  
  # Define Theme
  theme <- theme(panel.border =  element_blank(), axis.line =  element_line(), panel.background=element_rect(fill='white'), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_blank(), axis.title.x=element_blank(), legend.text=element_text(colour='white'), legend.title=element_text(colour='white'), legend.key=element_rect(colour='white', fill="white"))
  
  # Define additional parameters
  y_label <- ylab("Mutation Burden")
  legend <- scale_fill_manual(name="Translational Effect", values=c("Non Synonymous"="blue"))
  guide <- guides(fill=guide_legend(override.aes=list(fill="white")))
  
  # ggplot2 call
  p1 <- ggplot(x, aes(x=sample, y=mut_burden, fill=Type)) + geom_bar(stat='identity', alpha=.75, width=1) + theme + y_label + legend + guide
  
  return(p1)
}