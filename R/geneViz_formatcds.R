#' format cds
#'
#' given a Granges object specifying a region of interest, format into a form
#' recognizable by ggplot2
#' @name geneViz_formatCDS
#' @param txdb A TxDb object for a genome
#' @param gr A Granges object to format
#' @param genome Object of class BSgenome specifying the genome for GC content
#' calculation
#' @param reduce Boolean specifying whether to collapse isoforms in the Granges
#' object ROI
#' @return Object of class data frame

geneViz_formatCDS <- function(txdb=NULL, gr=NULL, genome=NULL, reduce=FALSE)
{
  # Extract the CDS for each isoform overlapping GRanges object
    cds <- geneViz_extrCDS(txdb, gr, reduce=reduce)

    if(is.null(cds)){
        return(NA)
    }
    
  # Calculate GC content for retrieved data
    cds <- sapply(cds, calcGC, genome=genome)

  # Coerce the relevant data in the Granges object to a data frame
    cds <- lapply(cds, Granges2dataframe)

  # if the data frame has a size format, else return a list of NA dataframes
    if(nrow(cds[[1]]) != 0)
    {
    # Format the CDS list
    cds <- lapply(cds, function(x){cbind(x, Type = c("CDS"))})
    cds <- lapply(cds, function(x){cbind(x, Upper = c(1))})
    cds <- lapply(cds, function(x){cbind(x, Lower = c(-1))})
    cds <- lapply(cds, function(x){cbind(x, Mid = c(0))})
    cds <- lapply(cds, function(x){cbind(x, segStart = c(min(x$start)))})
    cds <- lapply(cds, function(x){cbind(x, segEnd = c(max(x$end)))})

    } else {
        cds <- mapply(rbind, cds, NA)
    }

    return(cds)
}
