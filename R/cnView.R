#' plot raw copy number calls for a sample
#' 
#' given data frame with columns Chr, Coord, Tumor, Normal, Diff, p_value plot raw copy number calls
#' @name cnView
#' @param x object of class data frame with columns "chromosome", "coordinate", "cn", "p_value" consisting of raw cn calls from a single sample
#' @param y object of class data frame with columns "chrom", "chromStart", "chromEnd", "name", "gieStain" for plotting the ideogram (optional)
#' @param z object of class data frame with columns "chromosome", "start", "end", "segmean" consisting of cn segment calls (optional)
#' @param genome character string specifying UCSC genome to use, uneccessary if y is supplied, defaults to "hg19"
#' @param chr character string specifying UCSC chromosome to plot one of "chr..." or "all"
#' @param main.cnDiff Boolean specifying whether values in cn are copy number differences or actual copy numbers
#' @param ideo.chr_txt_angle integer specifying angle of text when plotting chromosome band text
#' @param ideo.chr_txt_size integer specifying size of text when plotting chromosome band text
#' @param main.layers additional ggplot2 layers for the main cn plot
#' @param ideo.layers additional ggplot2 layers for the ideogram
#' @return ggplot object
#' @examples
#' cnView(Luc2CNraw, chr='chr14', genome='hg19', ideo.chr_txt_size=4)
#' @import gridExtra
#' @export

cnView <- function(x, y=NULL, z=NULL, genome='hg19', chr='chr1', main.cnDiff=FALSE, ideo.chr_txt_angle=45, ideo.chr_txt_size=5, main.layers=NULL, ideo.layers=NULL)
{
  # Perform a basic quality check
  input <- cnView.qual(x, y, genome)
  x <- input[[1]]
  y <- input[[2]]
  
  # Obtain Cytogenetic Band information
  # use y input or query UCSC for the data if it's not preloaded
  preloaded <- c("hg38", "hg19", "mm10", "mm9", "rn5")
  if(is.null(y) && any(genome == preloaded))
  {
    message("genome specified is preloaded, retrieving data...")
    cytobands <- GenVisR::cytoGeno[GenVisR::cytoGeno$genome == genome,]
    cytobands <- cytobands[,-which(colnames(cytobands) == "genome")]
  } else if(is.null(y)){
    # obtain the cytogenetic band information for the requested reference
    message("attempting to query UCSC sql database for cytogenic band information")
    cytobands <- suppressWarnings(get_cytobands(genome=genome))
  } else {
    message("using y for cytogenic band information...")
    cytobands <- y
  }
  
  # Create Dummy data and add to x for proper plot dimensions
  fakeStart <- aggregate(data=cytobands, FUN=min, chromStart~chrom)
  colnames(fakeStart) <- c("chromosome", "coordinate")
  fakeEnd <- aggregate(data=cytobands, FUN=max, chromEnd~chrom)
  colnames(fakeEnd) <- c("chromosome", "coordinate")
  dummyData <- rbind(fakeStart, fakeEnd)
  dummyData$chromosome <- as.factor(dummyData$chromosome)
  dummyData <- subset_chr(dummyData, chr)
  
  # Plot all chromosomes at once if specified
  if(chr == 'all')
  {
    p1 <- build.cnView.main(x, z=z, dummyData, chr=chr)
    return(p1)
  }
  
  # plot chromosome 
  chromosome_plot <- ideoView(cytobands, chromosome=chr, chr_txt_angle=ideo.chr_txt_angle, chr_txt_size=ideo.chr_txt_size, layers=ideo.layers)
  
  # if requested plot only selected chromosome
  x <- subset_chr(x, chr)
  
  # build the cn plot
  CN_plot <- build.cnView.main(x, dummyData, z=z, chr=chr, cnDiff=main.cnDiff, layers=main.layers)
  
  p1 <- align_y_cn(chromosome_plot, CN_plot)
  
  return(grid.arrange(p1))
}