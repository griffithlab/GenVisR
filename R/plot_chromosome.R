#' plot chromosome
#' 
#' plot chromosome in ggplot using UCSC database
#' @name plot_chromosome
#' @param genome character string specifying UCSC genome
#' @param chromosome character string specifying UCSC chromosome to plot
#' @param chr_txt_angle integer specifying angle of text when plotting band text
#' @param chr_txt_size integer specifying size of text when plotting band text
#' @return ggplot object

plot_chromosome <- function(genome='hg19', chromosome='chr1', chr_txt_angle=45, chr_txt_size=5)
{
  # obtain the cytogenetic band information for the requested reference
  cytobands <- get_cytobands(genome=genome)
  
  # format the data obtained from ucsc into something ggplot can understand
  cytobands <- format_cytobands(cytobands, chromosome=chromosome)
  
  # plot the chromosome
  chr_plot <- build_chromosome(cytobands, chromosome=chromosome, chr_txt_angle=chr_txt_angle, chr_txt_size=chr_txt_size)
  
  return(chr_plot)
}