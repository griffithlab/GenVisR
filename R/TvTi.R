#' plot transitions/transversions
#' 
#' Given a data frame with columns reference, variant, and sample construct a transition/transversion plot
#' @name TvTi
#' @param x Object of class data frame containing columns 'sample', reference', 'variant'
#' @param y named vector containing expected transition/transversion rates
#' @param type Object of class character specifying whether to plot the Proportion or Frequency, one of "Prop"
#' @param label_x_axis boolean specifying wheter to label x axis
#' @param x_axis_text_angle Integer specifying the angle to labels on x_axis
#' @param palette Character vector of length 6 specifying colors for trans/tranv type
#' @param file_type Character string specifying the format the input is in, one of 'MAF', 'MGI'
#' @return Object of class data frame with indels removed
#' @import plyr
#' @export

TvTi <- function(x, y=NULL, type='Proportion', label_x_axis=TRUE, x_axis_text_angle=45, palette=c('#7BC374', '#EFCD8D', '#8763A0', '#6677A0', '#EDEE8D', '#EF8D8D'), file_type='MAF')
{
  # Perform quality checks
  if(!is.null(y))
  {
    out <- TvTi.qual(x, y, file_type=file_type)
    x <- out$input1
    y <- out$input2
  } else {
    x <- TvTi.qual(x, file_type=file_type)
  }
  
  # add transition/transversion info
  x <- adply(x, 1, anno_trans_tranv)
  
  # Calculate the proportion of transitions/transversions
  x <- calc_trans_tranv_freq(x)
  
  #re-level based on proportion values
  sample_order <- x[order(x$trans_tranv, x$Prop),]
  sample_order <- unique(sample_order$sample)
  x$sample <- factor(x$sample, levels=sample_order)
  
  # Build the Transition/Transversion Plot
  p1 <- build_trans_tranv(x, type=type, x_axis_text_angle=x_axis_text_angle, palette=palette, label_x_axis=label_x_axis)
  
  return(p1)
}