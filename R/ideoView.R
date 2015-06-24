#' plot chromosome
#' 
#' plot chromosome in ggplot using UCSC database
#' @name ideoView
#' @param genome character string specifying UCSC genome, not required if y is specified
#' @param y data frame specifying cytogenetic information for a genome, should contain columns "chrom", "chromStart", "chromEnd", "name", "gieStain"
#' @param chromosome character string specifying UCSC chromosome to plot
#' @param chr_txt_angle integer specifying angle of text when plotting band text
#' @param chr_txt_size integer specifying size of text when plotting band text
#' @return ggplot object
#' @export

ideoView <- function(genome='hg19', y=NULL, chromosome='chr1', chr_txt_angle=45, chr_txt_size=5)
{
  # use y input or query UCSC for the data
  if(is.null(y))
  {
    # obtain the cytogenetic band information for the requested reference
    cytobands <- get_cytobands(genome=genome)
  } else {
    cytobands <- y
  }

  # format the data obtained from ucsc into something ggplot can understand
  cytobands <- format_cytobands(cytobands, chromosome=chromosome)
  
  # plot the chromosome
  chr_plot <- build.ideogram(cytobands, chromosome=chromosome, chr_txt_angle=chr_txt_angle, chr_txt_size=chr_txt_size)
  
  return(chr_plot)
}