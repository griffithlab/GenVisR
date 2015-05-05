#' Extract 3' UTR
#' 
#' Extract 3' UTR coordinates within a GRanges object given a transcription database
#' @name extrthreeUTR
#' @param txdb A TxDb object for a genome
#' @param gr A Granges object specifying the region of interest
#' @param reduce Boolean specifying whether to collapse isoforms
#' @param gaps Boolean specifying whether to report space between 3'UTR instead of 3'UTR
#' @return Object of class Granges list

extrthreeUTR <- function(txdb, gr, reduce=FALSE, gaps=FALSE)
{
  # get a list of transcript id's overlapping the Granges object
  transcripts <- transcriptsByOverlaps(txdb, gr)
  map <- relist(unlist(transcripts, use.names=FALSE)$tx_id, transcripts)
  txid <- unlist(map, use.names=FALSE)
  
  # extract 3' UTR from transcript database given transcript ID
  threeUTR <- threeUTRsByTranscript(txdb)
  threeUTR <- threeUTR[names(threeUTR) %in% txid]
  
  # reduce isoforms into one if set to true and convert to GRanges list
  if(reduce==TRUE)
  {
    threeUTR <- reduce(unlist(threeUTR))
    threeUTR <- GRangesList(threeUTR)
  }
  
  # If Granges object is not an entire gene zoom in to the region specified
  threeUTR <- lapply(threeUTR, intersect, gr)
  
  # If gaps is set to true report report gaps between cds in lieu of cds
  if(gaps == TRUE)
  {
    # Calcluate gaps between cds regions for each isoform
    threeUTR <- lapply(threeUTR, gaps)
    
    # Limit the caluclated gaps to just the chromosomal region of interest
    threeUTR <- lapply(threeUTR, function(x) x[seqnames(x) == as.character(seqnames(gr))])
    
    # Limit the calculated gaps to just the strand of interest
    threeUTR <- lapply(threeUTR, function(x) x[strand(x) == as.character(strand(gr))])
  }
  
  return(threeUTR)
}