#' disjoin copynumber segments
#'
#' split genomic segments so that none are overlapping
#' @name cnFreq_disjoin
#' @noRd
#' @param x Object of class data frame with columns chromosome, start, end, segmean, and sample
#' @return Object of class data frame with disjoint genomic segments
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomicRanges mcols
#' @importFrom GenomicRanges disjoin
#' @importFrom IRanges IRanges
#' @importFrom IRanges extractList

cnFreq_disjoin <- function(x){
    # create the Granges object for the data
    x <- GenomicRanges::GRanges(seqnames=x$chromosome,
                                ranges=IRanges::IRanges(start=x$start, end=x$end),
                                "sample"=x$sample, "segmean"=x$segmean)
    
    # disjoin with grange, get a mapping of meta columns and expand it
    disJoint_x <- GenomicRanges::disjoin(x, with.revmap=TRUE)
    revmap <- GenomicRanges::mcols(disJoint_x)$revmap
    disJoint_x <- rep(disJoint_x, lengths(revmap))

    
    # exract the meta columns and map them back to the disJoint GRanges object
    sample <- unlist(IRanges::extractList(GenomicRanges::mcols(x)$sample, revmap))
    segmean <- unlist(IRanges::extractList(GenomicRanges::mcols(x)$segmean, revmap))
    GenomicRanges::mcols(disJoint_x)$sample <- sample
    GenomicRanges::mcols(disJoint_x)$segmean <- segmean
    
    # convert the GRanges Object back to a data frame
    disJoint_x <- as.data.frame(disJoint_x)[,c("seqnames", "start", "end", "width",
                                               "sample", "segmean")]
    colnames(disJoint_x) <- c("chromosome", "start", "end", "width", "sample", "segmean")
    return(disJoint_x)
}
