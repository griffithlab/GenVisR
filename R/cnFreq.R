#' Construct copy number frequency plot
#'
#' given a matrix construct a plot to display copy number changes across the genome for a group of samples
#' @name cnFreq
#' @param x object of class data frame with rows representing the proportion of CN losses/gains across the genome (default), or actual CN values. The formor option contains columns chromosome, start, end, gain, and loss, and the latter option contains columns chromosome, start, end, segmean, and sample
#' @param CN_low_cutoff numeric value representing the cutoff for copy number losses. Only used if x represents CN values.
#' @param CN_high_cutoff numeric value representing the cutoff for copy number gains. Only used if x represents CN values.
#' @param plot_title character string for title of plot
#' @param CN_low_colour character string specifying low value of colour gradient 
#' @param CN_high_colour character string specifying high value of colour gradient
#' @param x_lab_size integer specifying the size of the x label
#' @param y_lab_size integer specifying the size of the y label
#' @param facet_lab_size integer specifying the size of the faceted labels
#' @param layers a valid ggplot layer to over-ride default parameters
#' @return ggplot object
#' @examples
#' # Create data
#' xstart <- seq(0,4990000,length.out=500)
#' xloss <- rep(runif(10,0,0.6),rep(50,10))/1.5
#' xloss <- xloss + jitter(xloss,amount=0.002)
#' x <- data.frame(chromosome=rep(paste0("chr",1:5),rep(500,5)), start=xstart, 
#' end=xstart+10000, loss=xloss, gain=(1-xloss))
#' 
#' cnFreq(x)
#' @export
#' @import ggplot2
#' @importFrom "gtools" mixedsort

cnFreq <- function(x, CN_low_cutoff=1.5, CN_high_cutoff=2.5, plot_title=NULL, CN_low_colour='#002EB8', CN_high_colour='#A30000', x_lab_size=12, y_lab_size=12, facet_lab_size=10, layers=NULL)
{
  # Perform quality check on input data
  data <- cnFreq.qual(x)
  x <- data[[1]]
  plotType <- data[[2]]
  
  # If x contains actual CN values, transform into frequencies
  if(plotType=="freq")
  {
    xuniq <- unique(x[,c("chromosome","start","end")])
    gain = loss = obs = rep(NA, nrow(xuniq))
    for(i in 1:nrow(xuniq)) {
      tmpind = Reduce(intersect, list(which(x$chromosome==xuniq[i,1]),which(x$start==xuniq[i,2]),which(x$end==xuniq[i,3])))
      gain[i] = sum(x[tmpind,"segmean"] >= CN_high_cutoff, na.rm=TRUE)
      loss[i] = sum(x[tmpind,"segmean"] <= CN_low_cutoff, na.rm=TRUE)
      obs[i] = length(tmpind)
    }
    x <- data.frame(chromosome=xuniq$chromosome, start=xuniq$start, end=xuniq$end, gain=gain, loss=loss, obs=obs)
  }
  
  # Transform losses to be negative values
  x$loss <- -1*x$loss
  x <- as.data.frame(x)
 
  # Maintain the order of chromosomes
  x$chromosome <- factor(x$chromosome, levels=mixedsort(unique(x$chromosome)))
  
  # Construct the plot
  p1 <- build.cnFreq(x, plotType, plot_title=plot_title, CN_low_colour=CN_low_colour, CN_high_colour=CN_high_colour, x_lab_size=x_lab_size, y_lab_size=y_lab_size, facet_lab_size=facet_lab_size, layers=layers)
  
  return(p1)
}