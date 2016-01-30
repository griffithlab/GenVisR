#' Extract UTR
#'
#' Extract UTR coordinates within a GRanges object given a transcription database
#' @name geneViz_extrUTR
#' @param txdb A TxDb object for a genome
#' @param gr A Granges object specifying the region of interest
#' @param reduce Boolean specifying whether to collapse isoforms
#' @param gaps Boolean specifying whether to report space between UTR instead of UTR
#' @return Object of class Granges list
#' @importFrom GenomicRanges GRanges
#' @importFrom BiocGenerics unlist
#' @importFrom IRanges IRanges
#' @importFrom GenomicRanges mcols
#' @importFrom GenomicFeatures transcriptsByOverlaps
#' @importFrom AnnotationDbi select
#' @importFrom GenomicRanges GRangesList
#' @importFrom GenomicRanges reduce
#' @importFrom GenomicRanges strand
#' @importFrom GenomicRanges intersect
#' @importFrom GenomicRanges gaps
#' @importFrom IRanges relist

geneViz_extrUTR <- function(txdb=txdb, gr=gr, reduce=FALSE, gaps=FALSE)
{
    # Declare internal annonymous fuctions
    n <- function(x)
    {
        return(as.numeric(x))
    }

    f1 <- function(r)
    {
        if (is.na(r['CDSSTART']))
        {
            r['UTRSTART'] <- n(r['EXONSTART'])
            r['UTREND'] <- n(r['EXONEND'])
        }else if(r['CDSSTART'] != r['EXONSTART']){
            r['UTRSTART'] <- n(r['EXONSTART'])
            r['UTREND'] <- n(r['CDSSTART']) - 1
        }else{
            r['UTRSTART'] <- n(r['CDSEND']) + 1
            r['UTREND'] <- n(r['EXONEND'])
        }

        return(r)
    }

    f2 <- function(r)
    {
        g <- GenomicRanges::GRanges(seqnames=BiocGenerics::unlist(r['TXCHROM']),
                                    ranges=IRanges::IRanges(start = as.numeric(BiocGenerics::unlist(r['UTRSTART'])),
                                                            end = as.numeric(BiocGenerics::unlist(r['UTREND']))),
                                    strand=BiocGenerics::unlist(r['TXSTRAND']),
                                    txname=BiocGenerics::unlist(r['TXNAME']),
                                    exonrank=BiocGenerics::unlist(r['EXONRANK']))
        return(g)
    }

    f3 <- function(x)
    {
        x$txname[[1]]
    }

    f4 <- function(x)
    {
        return(length(x) > 0)
    }

    f5 <- function(gr, name)
    {
        GenomicRanges::mcols(gr)$txname <- name
        return(gr)
    }

    message("Obtaining UTR Coordinates")

    # get a list of transcript id's overlapping the Granges object
    transcripts <- GenomicFeatures::transcriptsByOverlaps(txdb, gr)
    map <- IRanges::relist(BiocGenerics::unlist(transcripts, use.names=FALSE)$tx_id,
                           transcripts)
    txid <- BiocGenerics::unlist(map, use.names=FALSE)

    # extract UTR from transcript database given transcript ID
    r <- AnnotationDbi::select(txdb, as.character(txid),
                               c("CDSSTART","CDSEND","TXCHROM","TXSTRAND",
                                 "CDSID","TXID", "EXONRANK","TXNAME",
                                 "EXONSTART","EXONEND"),
                               "TXID")

    idx <- as.vector(r['CDSSTART'] != r['EXONSTART'] | r['CDSEND'] != r['EXONEND'] | is.na(r['CDSSTART']))

    if (!any(idx))
    {
        return(NA)
    }

    r <- r[idx,]

    df <- as.data.frame(t(apply(r,1,f1)))
    r['UTRSTART'] <- as.numeric(as.character(df$UTRSTART))
    r['UTREND'] <- as.numeric(as.character(df$UTREND))
    r <- split(r, r['TXID'])

    UTR <- GenomicRanges::GRangesList(BiocGenerics::unlist(lapply(r, f2)))

    txnames <- lapply(UTR, f3)

    # reduce isoforms into one if set to true and convert to GRanges list
    if(reduce==TRUE)
    {
        UTR <- GenomicRanges::reduce(BiocGenerics::unlist(UTR))
        UTR <- GenomicRanges::GRangesList(UTR)
    }

    # If Granges object is not an entire gene zoom in to the region specified
    if(as.character(GenomicRanges::strand(gr)) == '*')
    {
        GenomicRanges::strand(gr) <- '+'
        UTR_forward <- lapply(UTR, GenomicRanges::intersect, gr)

        GenomicRanges::strand(gr) <- '-'
        UTR_reverse <- lapply(UTR, GenomicRanges::intersect, gr)

        UTR <- c(UTR_forward, UTR_reverse)
    } else {
        UTR <- lapply(UTR, GenomicRanges::intersect, gr)
    }

    idx <- as.vector(BiocGenerics::unlist(lapply(UTR, f4)))
    UTR <- UTR[idx]

    if(length(UTR) == 0)
    {
        return(NA)
    }

    # If gaps is set to true report report gaps between cds in lieu of cds
    if(gaps == TRUE)
    {
        message("Calculating UTR Gaps")

        # Calcluate gaps between cds regions for each isoform
        UTR <- lapply(UTR, GenomicRanges::gaps)

        # Limit the caluclated gaps to just the chromosomal region of interest
        UTR <- lapply(UTR, function(x) x[GenomicRanges::seqnames(x) == as.character(GenomicRanges::seqnames(gr))])

        # Limit the calculated gaps to just the strand of interest
        UTR <- lapply(UTR, function(x) x[GenomicRanges::strand(x) == as.character(GenomicRanges::strand(gr))])
    }

    if (reduce==FALSE)
    {
        keys <- names(UTR)
        UTR <- mapply(f5, UTR[keys], txnames[keys], SIMPLIFY=FALSE)
    }else{
        GenomicRanges::mcols(UTR[[1]])$txname <- 'merged'
        names(UTR) <- 'merged'
    }

    return(UTR)
}
