#' build transitions/transversions
#' 
#' Given a data frame with columns 'trans_tranv', 'sample', 'Freq', and 'Prop', build a transition/transversion plot
#' @name build_trans_tranv
#' @param x Object of class data frame containing columns 'trans_tranv', 'sample', 'Freq', and 'Prop'
#' @param type Object of class character specifying whether to plot the Proportion or Frequency, one of "Prop"
#' @param label_x_axis boolean specifying wheter to label x axis
#' @param x_axis_text_angle Integer specifying the angle to labels on x_axis
#' @param palette Character vector of length 6 specifying colors for trans/tranv type
#' @return GGplot Object
#' @import ggplot2

build_trans_tranv <- function(x, type='Proportion', label_x_axis=TRUE, x_axis_text_angle=45, palette=c('#7BC374', '#EFCD8D', '#8763A0', '#6677A0', '#EDEE8D', '#EF8D8D'))
{
  # Define various parameters of plot
  if(type == 'Proportion')
  {
    bar <- geom_bar(data=x, mapping=aes(x=sample, y=Prop, fill=trans_tranv), stat='identity', width=1)
  } else if(type == 'Frequency')
  {
    bar <- geom_bar(data=x, mapping=aes(x=sample, y=Freq, fill=trans_tranv), stat='identity', width=1)
  }
  ylabel <- ylab(type)
  xlabel <- xlab(paste0("Sample: n=", length(unique(x$sample))))
  fill_palette <- scale_fill_manual(name='Transistion/Transverstion', values=palette)
  
  # Define theme of plot
  if(label_x_axis == TRUE)
  {
    theme <- theme(axis.text.x=element_text(angle=x_axis_text_angle, hjust=1, vjust=1))
  } else {
    theme <- theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
  }
  
  
  # Define plot
  p1 <- ggplot() + bar + xlabel + ylabel + theme + fill_palette
  
  return(p1)
}