#' format 5' UTR
#' 
#' given a Granges object specifying a region of interest, format into a form recognizable by ggplot2
#' @name formatfiveUTR
#' @param txdb A TxDb object for a genome
#' @param gr A Granges object to format
#' @param genome Object of class BSgenome specifying the genome for GC content calculation
#' @param reduce Boolean specifying whether to collapse isoforms in the Granges object ROI
#' @return Object of class data frame

formatfiveUTR <- function(txdb, gr, genome, reduce=FALSE)
{
  # Extract the CDS for each isoform overlapping GRanges object
  fiveUTR <- extrfiveUTR(txdb, gr, reduce=reduce)
  
  # Calculate GC content for retrieved data
  fiveUTR <- sapply(fiveUTR, calcGC, genome=genome)
  
  # Coerce the relevant data in the Granges object to a data frame
  fiveUTR <- lapply(fiveUTR, Granges2dataframe)
  
  # if the data frame has a size format, else return a list of NA dataframes
  if(nrow(fiveUTR[[1]]) != 0)
  {
    # Format the CDS list
    fiveUTR <- lapply(fiveUTR, function(x){cbind(x, Type = c("UTR"))})
    fiveUTR <- lapply(fiveUTR, function(x){cbind(x, Upper = c(.5))})
    fiveUTR <- lapply(fiveUTR, function(x){cbind(x, Lower = c(-.5))})
    fiveUTR <- lapply(fiveUTR, function(x){cbind(x, Mid = c(0))})
    fiveUTR <- lapply(fiveUTR, function(x){cbind(x, segStart = c(min(x$start)))})
    fiveUTR <- lapply(fiveUTR, function(x){cbind(x, segEnd = c(max(x$end)))})
  } else {
    fiveUTR <- mapply(rbind, fiveUTR, NA)
  }
  
  return(fiveUTR)
}