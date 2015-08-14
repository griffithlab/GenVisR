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
    # add in fake column for legend
    # (necessary to have legend for proper plot alignment)
    x$Type <- c("Undefined")
    x$Type <- factor(x$Type,
                     levels=c("Synonymous", "Non Synonymous", "Undefined"))
    
    # Define Theme
    theme <- theme(axis.ticks.x=element_blank(),
                   axis.text.x=element_blank(),
                   axis.title.x=element_blank(),
                   legend.title=element_text(size=14))
    
    # Define additional parameters
    y_label <- ylab("Mutation Burden")
    legend <- scale_fill_manual(name="Translational Effect",
                                values=c("Non Synonymous"="blue",
                                         "Synonymous"="red",
                                         "Undefined"="black"),
                                drop=FALSE)
    
    if(!is.null(layers))
    {
        layers <- layers
    } else {
        layers <- geom_blank()
    }
    
    bar <- geom_bar(stat='identity', alpha=.75, width=1)
    
    # ggplot2 call
    p1 <- ggplot(x, aes_string(x='sample', y='mut_burden', fill='Type')) + bar +
        theme_bw() + theme + y_label + legend + layers
    
    return(p1)
}