#' plot transitions/transversions
#' 
#' Given a data frame with columns reference, variant, and sample construct a transition/transversion plot
#' @name TvTi
#' @param x Object of class data frame containing columns 'sample', reference', 'variant'
#' @param y named vector containing expected transition/transversion proportions with names "A->C or T->G", "A->G or T->C", "A->T or T->A", "G->A or C->T", "G->C or C->G", "G->T or C->A" or a data frame with column names "Prop", "trans_tranv" and levels of trans_tranv matching "A->C or T->G", "A->G or T->C", "A->T or T->A", "G->A or C->T", "G->C or C->G", "G->T or C->A"
#' @param type Object of class character specifying whether to plot the Proportion or Frequency, one of "Proportion", "Frequency"
#' @param label_x_axis boolean specifying wheter to label x axis
#' @param x_axis_text_angle Integer specifying the angle of labels on the x_axis
#' @param palette Character vector of length 6 specifying colors for each trans/tranv type
#' @param file_type Character string specifying the format the input is in, one of 'MAF', 'MGI'
#' @param tvti.layers Additional ggplot2 layers to add to the main panel plot
#' @param expec.layers Additional ggplot2 layers to add to the expected values plot
#' @param name_sort Boolean specifying whether to sort samples by name
#' @examples
#' TvTi(brcaMAF, type='Frequency', 
#' palette=c("#77C55D", "#A461B4", "#C1524B", "#93B5BB", "#4F433F", "#BFA753"), x_axis_text_angle=60)
#' @return ggplot2 object
#' @import plyr
#' @import gridExtra
#' @importFrom "gtools" mixedsort
#' @export


TvTi <- function(x, y=NULL, type='Proportion', label_x_axis=TRUE, x_axis_text_angle=45, palette=c("#2A2650", "#95fE52", "#F9C59B", "#0F722C", "#D7C7E9", "#FFB93F"), file_type='MAF', tvti.layers=NULL, expec.layers=NULL, name_sort=FALSE)
{ 
  # Perform quality checks
    out <- TvTi.qual(x, y, file_type=file_type)
    x <- out$input1
    y <- out$input2
  
  # add transition/transversion info
  message("annotating transitions and transversions")
  x <- adply(x, 1, anno_trans_tranv, .progress='text')
  
  # Calculate the proportion of transitions/transversions
  x <- calc_trans_tranv_freq(x)
  
  # re-level based on proportion values or via a smart sort if requested
  if(name_sort == TRUE)
  {
      sample_order <- unique(x$sample)
      sample_order <- mixedsort(sample_order)
      x$sample <- factor(x$sample, levels=sample_order)
  } else {
      sample_order <- x[order(x$trans_tranv, -x$Prop),]
      sample_order <- sample_order[sample_order$Prop != 0,]
      sample_order <- unique(sample_order$sample)
      x$sample <- factor(x$sample, levels=sample_order)
  }
  
  # Perform a quality control on y to ensure fill levels match x
  if(!is.null(y))
  {
    y$trans_tranv <- factor(y$trans_tranv, levels=levels(x$trans_tranv))
  }
  
  # Build the Transition/Transversion Plot
  p1 <- build.TvTi(x, y, type=type, x_axis_text_angle=x_axis_text_angle, palette=palette, label_x_axis=label_x_axis, tvti.layers=tvti.layers, expec.layers=NULL)
  
  if(!is.null(y))
  {
    # If y is input plot the expected values
    p2 <- build.TvTi(y, y, type=type, x_axis_text_angle=x_axis_text_angle, palette=palette, label_x_axis=label_x_axis, plot_expected=TRUE, tvti.layers=NULL, expec.layers=expec.layers)
    
    # Align the plots
    p3 <- align_y_TvTi(p1, p2)
    
    return(grid.arrange(p3))
  }

  
  return(p1)
}