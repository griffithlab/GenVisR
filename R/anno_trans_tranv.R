#' Annotate Transitions and Transversions
#' 
#' Given a data frame with columns reference and variant annotate the base change occurring
#' @name anno_trans_tranv
#' @param x Object of class data frame containing columns 'reference', 'variant'
#' @return Object of class data frame with transition/transversion annotations appended

anno_trans_tranv <- function(x)
{
  # add an extra column with the reference to variant
  x$base_change <- paste0(toupper(x$reference), "2", toupper(x$variant))
  
  # annotate the grouping of the base change
  x$trans_tranv <- switch(x$base_change, A2C="A<-C or T<-G", T2G="A<-C or T<-G", A2G="A<-G or T<-C", T2C="A<-G or T<-C", A2T="A<-T or T<-A", T2A="A<-T or T<-A", G2A="G<-A or C<-T", C2T="G<-A or C<-T", G2C="G<-C or C<-G", C2G="G<-C or C<-G", G2T="G<-T or C<-A", C2A="G<-T or C<-A")
  
  # remove the temp base change column
  x$base_change <- NULL
  return(x)
}