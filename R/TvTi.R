#' plot transitions/transversions
#' 
#' Given a data frame with columns reference, variant, and sample construct a transition/transversion plot
#' @name TvTi
#' @param x Object of class data frame containing columns 'sample', reference', 'variant'
#' @param y named vector containing expected transition/transversion proportions with names "A->C or T->G", "A->G or T->C", "A->T or T->A", "G->A or C->T", "G->C or C->G", "G->T or C->A" or a data frame with column names "Prop", "trans_tranv" and levels of trans_tranv matching "A->C or T->G", "A->G or T->C", "A->T or T->A", "G->A or C->T", "G->C or C->G", "G->T or C->A"
#' @param type Object of class character specifying whether to plot the Proportion or Frequency, one of "Prop"
#' @param label_x_axis boolean specifying wheter to label x axis
#' @param x_axis_text_angle Integer specifying the angle to labels on x_axis
#' @param palette Character vector of length 6 specifying colors for trans/tranv type
#' @param file_type Character string specifying the format the input is in, one of 'MAF', 'MGI'
#' @param layers Additional ggplot2 layers to add
#' @return Object of class data frame with indels removed
#' @import plyr
#' @import ggplot2
#' @export

TvTi <- function(x, y=NULL, type='Proportion', label_x_axis=TRUE, x_axis_text_angle=45, palette=c('#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e', '#e6ab02'), file_type='MAF', layers=NULL)
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
  
  # re-level based on proportion values
  sample_order <- x[order(x$trans_tranv, -x$Prop),]
  sample_order <- sample_order[sample_order$Prop != 0,]
  sample_order <- unique(sample_order$sample)
  x$sample <- factor(x$sample, levels=sample_order)
  
  # Perform a quality control on y to ensure fill levels match x
  if(!is.null(y))
  {
    y$trans_tranv <- factor(y$trans_tranv, levels=levels(x$trans_tranv))
  }
  
  # Build the Transition/Transversion Plot
  p1 <- build.TvTi(x, y, type=type, x_axis_text_angle=x_axis_text_angle, palette=palette, label_x_axis=label_x_axis, layers=layers)
  
  if(!is.null(y))
  {
    # If y is input plot the expected values
    p2 <- build.TvTi(y, y, type=type, x_axis_text_angle=x_axis_text_angle, palette=palette, label_x_axis=label_x_axis, plot_expected=TRUE, layers=layers)
    
    # Align the plots
    p3 <- align_y_TvTi(p1, p2)
    
    return(p3)
  }

  
  return(p1)
}