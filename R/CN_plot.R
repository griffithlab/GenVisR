#' plot CN values
#' 
#' given data frame with columns Chr, Coord, Tumor, Normal, Diff, p_value plot CN
#' @name CN_plot
#' @param x a data frame with columns Chr, Coord, Tumor, Normal, Diff, p_value
#' @param chr character string specifying UCSC chromosome to plot one of chr... or all
#' @param chr_txt_angle integer specifying angle of text when plotting band text
#' @param chr_txt_size integer specifying size of text when plotting band text
#' @return ggplot object
#' @export

CN_plot <- function(x, genome='hg19', chr='chr1', chr_txt_angle=45, chr_txt_size=5)
{
  
  if(chr == 'all')
  {
    x <- reformat_cn(x)
    p1 <- build_CN_plot(x, chr=chr)
    return(p1)
  }
  
  # plot chromosome 
  chromosome_plot <- plot_chromosome(genome=genome, chromosome=chr, chr_txt_angle=chr_txt_angle, chr_txt_size=chr_txt_size)
  
  # reformat cn data frame for ggplot
  x <- reformat_cn(x)
  
  # if requested plot only selected chromosome
  x <- subset_chr(x, chr)
  
  # build the cn plot
  CN_plot <- build_CN_plot(x, chr=chr)
  
  p1 <- align_y_cn(chromosome_plot, CN_plot)
  
  return(p1)
}