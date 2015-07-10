#' plot mutation burden
#' 
#' plot a barchart showing mutation burden given by data frame
#' @name build.mutRecurB.mutSpec
#' @param x a data frame containing columns sample, mut_burden
#' @param layers additional ggplot2 layers to plot
#' @return a ggplot object
#' @import ggplot2

build.mutRecurB.mutSpec <- function(x, layers=NULL)
{  
  # add in fake column for legend (necessary to have legend for proper plot alignment)
  # make everything white to hide legend
  x$Type <- "Non Synonymous"
  
  # Define Theme
  theme <- theme(axis.ticks.x=element_blank(), axis.text.x=element_blank(), axis.title.x=element_blank(), legend.text=element_text(colour='white'), legend.title=element_text(colour='white'), legend.key=element_rect(colour='white', fill="white"))
  
  # Define additional parameters
  y_label <- ylab("Mutation Burden")
  legend <- scale_fill_manual(name="Translational Effect", values=c("Non Synonymous"="blue"))
  guide <- guides(fill=guide_legend(override.aes=list(fill="white")))
  
  if(!is.null(layers))
  {
    layers <- layers
  } else {
    layers <- geom_blank()
  }
  
  # ggplot2 call
  p1 <- ggplot(x, aes(x=sample, y=mut_burden, fill=Type)) + geom_bar(stat='identity', alpha=.75, width=1) + theme_bw() + theme + y_label + legend + guide + layers
  
  return(p1)
}