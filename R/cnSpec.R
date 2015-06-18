#' Construct CN cohort plot
#' 
#' given a data frame construct a plot to display CN information for a group of samples
#' @name cnSpec
#' @param x object of class data frame containing columns chromosome, start, end, segmean, sample
#' @param y object of class data frame containing columns chromosome start end
#' @param genome character string specifying UCSC genome from which data is derived
#' @param plot_title character string for title of plot
#' @param background character string specifying backround color of plot
#' @param CN_low_colour character string specifying low value of colour gradient 
#' @param CN_high_colour character string specifying high value of colour gradient
#' @param x_lab_size integer specifying the size of the X label
#' @param Y_lab_size integer specifying the size of the Y label
#' @param facet_lab_size integer specifying the size of the faceted labels
#' @return ggplot object
#' @export
#' @import plyr
#' @import reshape2
#' @import gtools

cnSpec <- function(x, y=NULL, genome='hg19', plot_title=NULL, background='grey90', CN_low_colour='#002EB8', CN_high_colour='#A30000', x_lab_size=12, y_lab_size=12, facet_lab_size=10)
{
  # Perform quality check on input data
  cnSpec.qual(x, y)
  
  # Get dummy data for genome
  if(!is.null(y))
  {
    y
  } else {
    UCSC_Chr_pos <- CN_dummy_data(genome=genome)
  }
  
  
  # Dcast the input data into a recognizable format
  CN_data <- dcast(x, chromosome + start + end ~ sample, value.var = "segmean")
  
  # Rbind fill the dummy and CN data and remove any chr appendices
  CN_data <- rbind.fill(CN_data, UCSC_Chr_pos)
  CN_data[is.na(CN_data)] <- NA
  CN_data$chromosome <- gsub('chr', '', CN_data$chromosome)
  
  # melt the data for ggplot2 call
  CN_data <- melt(CN_data, id.vars=c('chromosome', 'start', 'end'))
  colnames(CN_data) <- c('chromosome', 'start', 'end', 'sample', 'cn')
  
  # Change the order of chromosomes and samples (natural sort order)
  chromosome_sorted <- unique(mixedsort(CN_data$chromosome))
  CN_data$chromosome <- factor(CN_data$chromosome, levels=chromosome_sorted)
  sample_sorted <- unique(mixedsort(CN_data$sample))
  CN_data$sample <- factor(CN_data$sample, levels=sample_sorted)
  
  # Construct the plot
  p1 <- buildCN_cohort(CN_data, plot_title=plot_title, background=background, CN_low_colour=CN_low_colour, CN_high_colour=CN_high_colour, x_lab_size=x_lab_size, y_lab_size=y_lab_size, facet_lab_size=facet_lab_size)
  
  return(p1)
}