#' cdsFromTXID
#'
#' Return CDS coordinates as a GRanges object given transcript IDs
#' @name cdsFromTXID
#' @param txdb A TxDb object for a genome
#' @param txid A list of TXIDs
#' @return Object of class Granges
#' @importFrom "IRanges" IRanges

cdsFromTXID <- function(txdb, txid){
    r <- select(txdb, as.character(txid),
                c("CDSSTART","CDSEND","CDSCHROM","CDSSTRAND","CDSID","TXID","EXONRANK","TXNAME"), "TXID")
    idx <- as.vector(!is.na(r['CDSSTART']))
    r <- r[idx,]
    if (!any(idx)){
        return(NA)
    }
    r <- split(r, r['TXID'])
    f <- function(r){
        g <- GRanges(seqnames = unlist(r['CDSCHROM']),
                     ranges = IRanges(start = as.numeric(unlist(r['CDSSTART'])), end = as.numeric(unlist(r['CDSEND']))),
                     strand = unlist(r['CDSSTRAND']),
                     txname=unlist(r['TXNAME']),
                     exonrank=unlist(r['EXONRANK']))
        return(g)
    }
    cds <- GRangesList(unlist(lapply(r, f)))
    return(cds)
}
