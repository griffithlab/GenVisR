#' plot chromosome
#' 
#' plot chromosome in ggplot using UCSC database
#' @name ideoView
#' @param x data frame specifying cytogenetic information for a genome, should contain columns "chrom", "chromStart", "chromEnd", "name", "gieStain"
#' @param chromosome character string specifying UCSC chromosome to plot
#' @param chr_txt_angle integer specifying angle of text when plotting band text
#' @param chr_txt_size integer specifying size of text when plotting band text
#' @param layers additional ggplot2 layers for the ideogram
#' @return ggplot object
#' @export

ideoView <- function(x, chromosome='chr1', chr_txt_angle=45, chr_txt_size=5, layers=NULL)
{
  cytobands <- x
  
  # format the data obtained from ucsc into something ggplot can understand
  cytobands <- format_cytobands(cytobands, chromosome=chromosome)
  
  # plot the chromosome
  chr_plot <- build.ideogram(cytobands, chromosome=chromosome, chr_txt_angle=chr_txt_angle, chr_txt_size=chr_txt_size, layers=layers)
  
  return(chr_plot)
}