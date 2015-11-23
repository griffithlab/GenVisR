#' Calculate GC content
#'
#' Calculate GC content for elements in a GRanges object
#' @name geneViz_calcGC
#' @param gr A Granges object to calculate GC content for
#' @param genome Object of class BSgenome specifying the genome to calculate GC
#' content from
#' @return Object of class GRanges

geneViz_calcGC <- function(gr, genome)
{
  # Calculate the GC content for the given Granges object
    GC <- as.data.frame(Biostrings::alphabetFrequency(Biostrings::getSeq(genome,
                                                                         gr),
                                                      as.prob=TRUE))
    GC <- GC$G + GC$C

  # append this information to the Granges object
    GenomicRanges::values(gr) <- cbind(as.data.frame(GC),
                                       as.data.frame(GenomicRanges::mcols(gr)))

    return(gr)
}
