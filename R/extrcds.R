#' Extract CDS
#' 
#' Extract CDS coordinates within a GRanges object given a transcription database
#' @name extrcds
#' @param txdb A TxDb object for a genome
#' @param gr A Granges object specifying the region of interest
#' @param reduce Boolean specifying whether to collapse isoforms
#' @param gaps Boolean specifying whether to report space between CDS instead of CDS
#' @return Object of class Granges list

extrcds <- function(txdb, gr, reduce=FALSE, gaps=FALSE)
{
  # get a list of transcript id's overlapping the Granges object
  transcripts <- transcriptsByOverlaps(txdb, gr)
  map <- relist(unlist(transcripts, use.names=FALSE)$tx_id, transcripts)
  txid <- unlist(map, use.names=FALSE)
  
  # extract CDS from transcript database given transcript ID
  cds <- cdsBy(txdb, "tx")
  cds <- cds[names(cds) %in% txid]
  
  # reduce isoforms into one if set to true and convert to GRanges list
  if(reduce==TRUE)
  {
    cds <- reduce(unlist(cds))
    cds <- GRangesList(cds)
  }
  
  # If Granges object is not an entire gene zoom in to the region specified
  cds <- lapply(cds, intersect, gr)
  
  # If gaps is set to true report report gaps between cds in lieu of cds
  if(gaps == TRUE)
  {
    # Calcluate gaps between cds regions for each isoform
    cds <- lapply(cds, gaps)
    
    # Limit the caluclated gaps to just the chromosomal region of interest
    cds <- lapply(cds, function(x) x[seqnames(x) == as.character(seqnames(gr))])
    
    # Limit the calculated gaps to just the strand of interest
    cds <- lapply(cds, function(x) x[strand(x) == as.character(strand(gr))])
  }
  
  return(cds)
}
