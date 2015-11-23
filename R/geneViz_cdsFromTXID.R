#' cdsFromTXID
#'
#' Return CDS coordinates as a GRanges object given transcript IDs
#' @name geneViz_cdsFromTXID
#' @param txdb A TxDb object for a genome
#' @param txid A list of TXIDs
#' @return Object of class Granges

geneViz_cdsFromTXID <- function(txdb, txid)
{
    # Declare Anonymous functions
    f <- function(r)
    {
        g <- GenomicRanges::GRanges(seqnames = GenomicRanges::unlist(r['CDSCHROM']),
                     ranges=IRanges::IRanges(start = as.numeric(GenomicRanges::unlist(r['CDSSTART'])), end = as.numeric(GenomicRanges::unlist(r['CDSEND']))),
                     strand=GenomicRanges::unlist(r['CDSSTRAND']),
                     txname=GenomicRanges::unlist(r['TXNAME']),
                     exonrank=GenomicRanges::unlist(r['EXONRANK']))
        return(g)
    }    
    
    r <- AnnotationDbi::select(txdb, as.character(txid), c("CDSSTART","CDSEND",
                                                           "CDSCHROM",
                                                           "CDSSTRAND","CDSID",
                                                           "TXID", "EXONRANK",
                                                           "TXNAME"), "TXID")
    
    idx <- as.vector(!is.na(r['CDSSTART']))
    r <- r[idx,]
    
    if (!any(idx))
    {
        return(NA)
    }
    
    r <- split(r, r['TXID'])

    cds <- GenomicRanges::GRangesList(GenomicRanges::unlist(lapply(r, f)))
    
    return(cds)
}
