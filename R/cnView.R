#' plot CN values
#' 
#' given data frame with columns Chr, Coord, Tumor, Normal, Diff, p_value plot CN
#' @name cnView
#' @param x a data frame with columns chromosome, coordinate, cn, p_value
#' @param y a data frame with columns "chrom", "chromStart", "chromEnd", "name", "gieStain"
#' @param z a data frame with columns "chromosome", "start", "end", "segmean"
#' @param genome character string specifying UCSC genome to use
#' @param chr character string specifying UCSC chromosome to plot one of chr... or all
#' @param main.cnDiff Boolean specifying whether values in cn are copy number differences or actual copy number
#' @param chr_txt_angle integer specifying angle of text when plotting band text
#' @param chr_txt_size integer specifying size of text when plotting band text
#' @return ggplot object
#' @export

cnView <- function(x, y=NULL, z=NULL, genome='hg19', chr='chr1', main.cnDiff=FALSE, chr_txt_angle=45, chr_txt_size=5)
{
  # Perform a basic quality check
  input <- cnView.qual(x, y)
  x <- input[[1]]
  y <- input[[2]]
  
  # Obtain Cytogenetic Band information
  # use y input or query UCSC for the data
  if(is.null(y))
  {
    # obtain the cytogenetic band information for the requested reference
    cytobands <- suppressWarnings(get_cytobands(genome=genome))
  } else {
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
  chromosome_plot <- ideoView(cytobands, chromosome=chr, chr_txt_angle=chr_txt_angle, chr_txt_size=chr_txt_size)
  
  # if requested plot only selected chromosome
  x <- subset_chr(x, chr)
  
  # build the cn plot
  CN_plot <- build.cnView.main(x, dummyData, z=z, chr=chr, cnDiff=main.cnDiff)
  
  p1 <- align_y_cn(chromosome_plot, CN_plot)
  
  return(p1)
}