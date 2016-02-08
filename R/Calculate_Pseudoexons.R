#' calculate pseudoexons
#' 
#' calculate pseudoexons i.e. an exonic region in one isoform but not another
#' @name Calculate_Pseudoexons
#' @param x object of class dataframe with column names start end 
#' @param chr character string specifying the chromosome for pseudoexon calculation
#' @return Object of class data frame

Calculate_Pseudoexons <- function(x, chr)
{
  require("plyr")
  
  # Create GRanges object
  gr <- GRanges(seqnames=chr, ranges=IRanges(start=x$start, end=x$end))
  # Add in meta data
  mcols(gr, use.names=TRUE) <- x[,3:ncol(x)]
  
  # Calculate gaps
  pseudo_exon <- gaps(gr)
  
  # Convert back to a data frame
  range <- as.data.frame(ranges(pseudo_exon))
  meta <- as.data.frame(mcols(pseudo_exon))
  pseudo_exon <- cbind(range, meta)
  pseudo_exon <- pseudo_exon[,c('start', 'end')]
  
  # Specify that the calculated gap is a pseudo_exon
  if(nrow(pseudo_exon) != 0)
  {
    pseudo_exon$Type <- 'pseudo_exon'
  }
  
  # Bind the pseudo exons with the genomic features and sort by start
  genomic_data <- rbind.fill(x, pseudo_exon)
  genomic_data <- genomic_data[order(genomic_data$start),]
  
  return(genomic_data)
}
