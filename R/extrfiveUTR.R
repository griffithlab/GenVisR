#' Extract 5' UTR
#' 
#' Extract 5' UTR coordinates within a GRanges object given a transcription database
#' @name extrfiveUTR
#' @param txdb A TxDb object for a genome
#' @param gr A Granges object specifying the region of interest
#' @param reduce Boolean specifying whether to collapse isoforms
#' @param gaps Boolean specifying whether to report space between 5'UTR instead of 5'UTR
#' @return Object of class Granges list

extrfiveUTR <- function(txdb, gr, reduce=FALSE, gaps=FALSE)
{
  message("Obtaining 5' UTR Coordinates")
  
  # get a list of transcript id's overlapping the GRanges object
  transcripts <- transcriptsByOverlaps(txdb, gr)
  map <- relist(unlist(transcripts, use.names=FALSE)$tx_id, transcripts)
  txid <- unlist(map, use.names=FALSE)
  
  # extract 5' UTR from transcript database given transcript ID
  fiveUTR <- fiveUTRsByTranscript(txdb)
  fiveUTR <- fiveUTR[names(fiveUTR) %in% txid]
  
  # reduce isoforms into one if set to true and convert to GRanges list
  if(reduce==TRUE)
  {
    fiveUTR <- reduce(unlist(fiveUTR))
    fiveUTR <- GRangesList(fiveUTR)
  }
  
  # If Granges object is not an entire gene zoom in to the region specified
  fiveUTR <- lapply(fiveUTR, intersect, gr)
  
  # If gaps is set to true report report gaps between cds in lieu of cds
  if(gaps == TRUE)
  {
    message("Calculating 5' UTR Gaps")
    
    # Calcluate gaps between cds regions for each isoform
    fiveUTR <- lapply(fiveUTR, gaps)
    
    # Limit the caluclated gaps to just the chromosomal region of interest
    fiveUTR <- lapply(fiveUTR, function(x) x[seqnames(x) == as.character(seqnames(gr))])
    
    # Limit the calculated gaps to just the strand of interest
    fiveUTR <- lapply(fiveUTR, function(x) x[strand(x) == as.character(strand(gr))])
  }
  
  return(fiveUTR)
}