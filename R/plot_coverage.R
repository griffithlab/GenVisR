#' produce a coverage plot
#' 
#' produce a coverage plot displaying gene and coverage information
#' @name plot_coverage
#' @param coverage_data named list containing data frames with columns stop and cov
#' @param txdb A TxDb object for a genome
#' @param gr A Granges object specifying a region of interest
#' @param genome Object of class BSgenome specifying the genome
#' @param reduce Boolean specifying whether to collapse isoforms in the ROI
#' @param gene_name character string specifying the name of the gene or ROI
#' @param bg_fill character string giving the colour to fill the label
#' @param text_fill character string giving the colour to fill the text
#' @param border character string specifying the colour to fill the border of the label
#' @param size integer specifying the size of the text within the label
#' @param width_ratio vector of length 2 giving the ratio of track labels to plot
#' @param colour character string specifying the color of the data in the plot
#' @param plot_type character string specifying one of line, area for data display
#' @return ggplot object
#' @export

plot_coverage <- function(coverage_data, txdb, gr, genome, reduce=F, gene_name='test', bg_fill="black", text_fill="white", border="black", size=10, width_ratio=c(1, 10), colour="blue", plot_type="line")
{
  # Obtain a plot for the gene overlapping the Granges object and covert to a named list
  gene <- gene_plot(txdb, gr, genome, reduce=reduce)
  gene_list <- list()
  gene_list[[gene_name]] <- gene
  
  # Obtain coverage plots for the data input as a list
  coverage_plot <- lapply(coverage_data, build_coverage, colour=colour, plot_type=plot_type)
  
  # Combine both gene and coverage plot lists
  merged_data <- c(gene_list, coverage_plot)
  
  # Plot the data on a track
  track_coverage_plot <- plot_track(merged_data, gene_name=gene_name, bg_fill=bg_fill, text_fill=text_fill, border=border, size=size, axis_align='width', width_ratio=width_ratio, nested_list=T)
  
  return(track_coverage_plot)
}