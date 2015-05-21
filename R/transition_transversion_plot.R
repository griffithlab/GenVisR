#' plot transitions/transversions
#' 
#' Given a data frame with columns reference, variant, and sample construct a transition/transversion plot
#' @name transition_transversion_plot
#' @param x Object of class data frame containing columns 'reference', 'variant'
#' @param type Object of class character specifying whether to plot the Proportion or Frequency, one of "Prop"
#' @param label_x_axis boolean specifying wheter to label x axis
#' @param x_axis_text_angle Integer specifying the angle to labels on x_axis
#' @param palette Character vector of length 6 specifying colors for trans/tranv type
#' @return Object of class data frame with indels removed
#' @import plyr
#' @export

transition_transversion_plot <- function(x, type='Proportion', label_x_axis=TRUE, x_axis_text_angle=45, palette=c('#7BC374', '#EFCD8D', '#8763A0', '#6677A0', '#EDEE8D', '#EF8D8D'))
{
  # Check that columns are named appropriatley, if not print error
  if(any(grepl('^reference$', colnames(x))) && any(grepl('^variant$', colnames(x))) && any(grepl('^sample$', colnames(x))))
  {
    message("Found appropriate columns")
  } else {
    stop("Could not find all columns requested, missing one of reference, variant, sample")
  } 
  
  # Remove any indels present in the data
  x <- rm.indel(x)
  
  # TODO Convert IUPAC codes
  
  # add transition/transversion info
  x <- adply(x, 1, anno_trans_tranv)
  
  # Calculate the proportion of transitions/transversions
  x <- calc_trans_tranv_freq(x)
  
  #re-level based on proportion values
  sample_order <- x[which(x$trans_tranv %in% "G<-T or C<-A"),]
  sample_order <- sample_order[order(sample_order$Prop),]
  x$sample <- factor(x$sample, levels=sample_order$sample)
  
  # Build the Transition/Transversion Plot
  p1 <- build_trans_tranv(x, type=type, x_axis_text_angle=x_axis_text_angle, palette=palette, label_x_axis=label_x_axis)
  
  return(p1)
}