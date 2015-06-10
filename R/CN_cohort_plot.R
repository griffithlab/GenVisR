#' Construct CN cohort plot
#' 
#' given a data frame construct a plot to display CN information for a group of samples
#' @name CN_cohort_plot
#' @param data_frame object of class data frame containing columns Chromosome, Start, Stop, SegMean, Sample.name
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

CN_cohort_plot <- function(data_frame, genome='hg19', plot_title=NULL, background='grey90', CN_low_colour='#002EB8', CN_high_colour='#A30000', x_lab_size=12, y_lab_size=12, facet_lab_size=10)
{
  # Get dummy data for genome
  UCSC_Chr_pos <- CN_dummy_data(genome=genome)
  
  # Dcast the input data into a recognizable format
  CN_data <- dcast(data_frame, Chromosome + Start + Stop ~ Sample.name, value.var = "SegMean")
  
  # Rbind fill the dummy and CN data and remove any chr appendices
  CN_data <- rbind.fill(CN_data, UCSC_Chr_pos)
  CN_data[is.na(CN_data)] <- NA
  CN_data$Chromosome <- gsub('chr', '', CN_data$Chromosome)
  
  # melt the data for ggplot2 call
  CN_data <- melt(CN_data, id.vars=c('Chromosome', 'Start', 'Stop'))
  colnames(CN_data) <- c('Chromosome', 'Start', 'Stop', 'Sample', 'CN')
  
  # Change the order of chromosomes and samples (natural sort order)
  chromosome_sorted <- unique(mixedsort(CN_data$Chromosome))
  CN_data$Chromosome <- factor(CN_data$Chromosome, levels=chromosome_sorted)
  sample_sorted <- unique(mixedsort(CN_data$Sample))
  CN_data$Sample <- factor(CN_data$Sample, levels=sample_sorted)
  
  # Construct the plot
  p1 <- buildCN_cohort(CN_data, plot_title=plot_title, background=background, CN_low_colour=CN_low_colour, CN_high_colour=CN_high_colour, x_lab_size=x_lab_size, y_lab_size=y_lab_size, facet_lab_size=facet_lab_size)
  
  return(p1)
}