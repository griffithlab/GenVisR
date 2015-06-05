#' plot mutation recurrence in genes
#' 
#' plot a bar graph displaying the percentage of samples with a mutation
#' @name plot_bar
#' @param data_frame a data frame in MAF format
#' @return a ggplot object
#' @import scales

plot_bar <- function(data_frame)
{
  #####################################################################################################################
  ####### Function to create the left margin bar plot, function plots the percentage of samples with a mutation #######
  #####################################################################################################################
  
  # Convert all silent mutations to Synonymous, and all else to non-synonymous
  data_frame$trv_type <- as.character(data_frame$trv_type)
  data_frame$trv_type[toupper(data_frame$trv_type) != toupper('silent')] <- 'Non Synonymous'
  data_frame$trv_type[toupper(data_frame$trv_type) == toupper('silent')] <- 'Synonymous'
  data_frame$trv_type <- factor(data_frame$trv_type, levels=c('Synonymous', 'Non Synonymous'))
  
  # Define the number of samples for the Percentage calculation (Note, to pass a variable outside of aes into aes it needs to be defined again)
  total_number_sample <- nlevels(data_frame$sample)
  
  # Define Theme and various other layers to be passed to ggplot
  theme <- theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank(), legend.position=('none'), axis.title.x=element_text(size=18))
  y_limits <- ylim(100, 0)
  y_label <- ylab('% Samples With Mutation')	
  legend <- scale_fill_manual(name="Translational Effect", values=c("red", "blue"), breaks=c('Synonymous', 'Non Synonymous'), drop=FALSE)
  
  # Plotting call
  p1 <- ggplot(na.omit(data_frame), aes(x=gene, total_number_sample=total_number_sample, y=(..count..)/total_number_sample * 100, fill=trv_type), environment = environment()) + geom_bar(position='stack', alpha=.75) + coord_flip() + theme + y_label + scale_y_reverse() + legend	
  
  return(p1)
}