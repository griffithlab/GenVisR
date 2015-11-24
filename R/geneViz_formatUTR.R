#' format UTR
#'
#' given a Granges object specifying a region of interest, format into a form
#' recognizable by ggplot2
#' @name geneViz_formatUTR
#' @param txdb A TxDb object for a genome
#' @param gr A Granges object to format
#' @param genome Object of class BSgenome specifying the genome for GC content
#' calculation
#' @param reduce Boolean specifying whether to collapse isoforms in the Granges
#' object ROI
#' @return Object of class data frame

geneViz_formatUTR <- function(txdb=NULL, gr=NULL, genome=NULL, reduce=FALSE)
{
  # Extract the CDS for each isoform overlapping GRanges object
    UTR <- geneViz_extrUTR(txdb=txdb, gr=gr, reduce=reduce)

    if(is.null(UTR)){
        return(NA)
    }
  # Calculate GC content for retrieved data
    UTR <- sapply(UTR, geneViz_calcGC, genome=genome)

  # Coerce the relevant data in the Granges object to a data frame
    UTR <- lapply(UTR, geneViz_Granges2dataframe)

  # if the data frame has a size format, else return a list of NA dataframes
    if(nrow(UTR[[1]]) != 0)
    {
    # Format the CDS list
    UTR <- lapply(UTR, function(x){cbind(x, Type = c("UTR"))})
    UTR <- lapply(UTR, function(x){cbind(x, Upper = c(.5))})
    UTR <- lapply(UTR, function(x){cbind(x, Lower = c(-.5))})
    UTR <- lapply(UTR, function(x){cbind(x, Mid = c(0))})
    UTR <- lapply(UTR, function(x){cbind(x, segStart = c(min(x$start)))})
    UTR <- lapply(UTR, function(x){cbind(x, segEnd = c(max(x$end)))})
    } else {
        UTR <- mapply(rbind, UTR, NA)
    }

    return(UTR)
}
