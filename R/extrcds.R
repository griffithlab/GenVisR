#' Extract CDS
#' 
#' Extract CDS coordinates within a GRanges object given a transcription database
#' @name extrCDS
#' @param txdb A TxDb object for a genome
#' @param gr A Granges object specifying the region of interest
#' @param reduce Boolean specifying whether to collapse isoforms
#' @param gaps Boolean specifying whether to report space between CDS instead of CDS
#' @return Object of class Granges list

extrCDS <- function(txdb, gr, reduce=FALSE, gaps=FALSE)
{
  message("Obtaining CDS Coordinates")
  
  # get a list of transcript id's overlapping the Granges object
  transcripts <- transcriptsByOverlaps(txdb, gr)
  map <- relist(unlist(transcripts, use.names=FALSE)$tx_id, transcripts)
  txid <- unlist(map, use.names=FALSE)
  
  # extract CDS from transcript database given transcript ID
  cds <- cdsFromTXID(txdb, txid)
  f1 <- function(x){x$txname[[1]]}
  txnames <- lapply(cds, f1)
  
  if(typeof(cds) != 'S4'){
    return(NA)
  }
  
  # reduce isoforms into one if set to true and convert to GRanges list
  if(reduce==TRUE)
  {
    cds <- reduce(unlist(cds))
    cds <- GRangesList(cds)
  }
  
  # If Granges object is not an entire gene zoom in to the region specified
  cds <- lapply(cds, intersect, gr)
  f2 <- function(x){
    return(length(x) > 0)
  }
  idx <- as.vector(unlist(lapply(cds, f2)))
  cds <- cds[idx]
  if(length(cds) == 0){
    return(NA)
  }
  
  # If gaps is set to true report report gaps between cds in lieu of cds
  if(gaps == TRUE)
  {
    message("Calculating CDS Gaps")
    
    # Calcluate gaps between cds regions for each isoform
    cds <- lapply(cds, gaps)
    
    # Limit the caluclated gaps to just the chromosomal region of interest
    cds <- lapply(cds, function(x) x[seqnames(x) == as.character(seqnames(gr))])
    
    # Limit the calculated gaps to just the strand of interest
    cds <- lapply(cds, function(x) x[strand(x) == as.character(strand(gr))])
  }
  
  keys <- names(cds)
  f3 <- function(gr, name){
    mcols(gr)$txname <- name
    return(gr)
  }
  cds <- mapply(f3, cds[keys], txnames[keys], SIMPLIFY=F)
  
  return(cds)
}
