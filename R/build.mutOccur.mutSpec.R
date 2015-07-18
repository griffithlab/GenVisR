#' plot mutation recurrence in genes
#' 
#' plot a bar graph displaying the percentage of samples with a mutation
#' @name build.mutOccur.mutSpec
#' @param data_frame a data frame in MAF format
#' @param layers additional ggplot2 layers
#' @return a ggplot object
#' @import plyr

build.mutOccur.mutSpec <- function(data_frame, layers=NULL)
{
    # Convert all silent mutations to Synonymous, and all else to non-synonymous
    data_frame$trv_type <- as.character(data_frame$trv_type)
    data_frame$trv_type[toupper(data_frame$trv_type) != toupper('silent')] <- 'Non Synonymous'
    data_frame$trv_type[toupper(data_frame$trv_type) == toupper('silent')] <- 'Synonymous'
    data_frame$trv_type <- factor(data_frame$trv_type, levels=c('Synonymous', 'Non Synonymous'))
  
    # Define the number of samples for the Percentage calculation (Note, to pass a variable outside of aes into aes it needs to be defined again)
    total_number_sample <- nlevels(data_frame$sample)
  
    # Make a count column
    data_frame <- count(data_frame, c("gene", "trv_type"))
    data_frame$prop <- data_frame$freq/total_number_sample * 100
  
    # Define Theme and various other layers to be passed to ggplot
    theme <- theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank(), legend.position=('none'))
    y_limits <- ylim(100, 0)
    y_label <- ylab('% Samples With Mutation')	
    legend <- scale_fill_manual(name="Translational Effect", values=c("red", "blue"), breaks=c('Synonymous', 'Non Synonymous'), drop=FALSE)
    if(!is.null(layers))
    {
        layers <- layers
    } else {
        layers <- geom_blank()
    }
  
    # Plotting call
    p1 <- ggplot(na.omit(data_frame), aes_string(x='gene', y='prop', fill='trv_type')) + geom_bar(position='stack', alpha=.75, width=1, stat='identity') + theme_bw() + coord_flip() + theme + y_label + scale_y_reverse() + legend + layers	
  
    return(p1)
}