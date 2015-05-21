#' Extract UTR
#' 
#' Extract UTR coordinates within a GRanges object given a transcription database
#' @name extrUTR
#' @param txdb A TxDb object for a genome
#' @param gr A Granges object specifying the region of interest
#' @param reduce Boolean specifying whether to collapse isoforms
#' @param gaps Boolean specifying whether to report space between UTR instead of UTR
#' @return Object of class Granges list

extrUTR <- function(txdb, gr, reduce=FALSE, gaps=FALSE)
{
  message("Obtaining UTR Coordinates")
  
  # get a list of transcript id's overlapping the Granges object
  transcripts <- transcriptsByOverlaps(txdb, gr)
  map <- relist(unlist(transcripts, use.names=FALSE)$tx_id, transcripts)
  txid <- unlist(map, use.names=FALSE)
  
  # extract UTR from transcript database given transcript ID
  r <- select(txdb, as.character(txid), 
                     c("CDSSTART","CDSEND","TXCHROM","TXSTRAND","CDSID","TXID","EXONRANK","TXNAME","EXONSTART","EXONEND"), "TXID")
  idx <- as.vector(r['CDSSTART'] != r['EXONSTART'] | r['CDSEND'] != r['EXONEND'] | is.na(r['CDSSTART']))
  if (!any(idx)){
    return(NA)
  }
  r <- r[idx,]
  n <- function(x){return(as.numeric(x))}
  f <- function(r){
    if (is.na(r['CDSSTART'])){
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
  df <- as.data.frame(t(apply(r,1,f)))
  r['UTRSTART'] <- as.numeric(as.character(df$UTRSTART))
  r['UTREND'] <- as.numeric(as.character(df$UTREND))
  r <- split(r, r['TXID'])
  f <- function(r){
    g <- GRanges(seqnames = unlist(r['TXCHROM']), 
                 ranges = IRanges(start = as.numeric(unlist(r['UTRSTART'])), end = as.numeric(unlist(r['UTREND']))),
                 strand = unlist(r['TXSTRAND']),
                 txname=unlist(r['TXNAME']),
                 exonrank=unlist(r['EXONRANK']))
    return(g)
  }
  UTR <- GRangesList(unlist(lapply(r, f)))
  
  # reduce isoforms into one if set to true and convert to GRanges list
  if(reduce==TRUE)
  {
    UTR <- reduce(unlist(UTR))
    UTR <- GRangesList(UTR)
  }
  
  # If Granges object is not an entire gene zoom in to the region specified
  UTR <- lapply(UTR, intersect, gr)
  
  f <- function(x){
    return(length(x) > 0)
  }
  idx <- as.vector(unlist(lapply(UTR, f)))
  UTR <- UTR[idx]
  if(length(UTR) == 0){
    return(NA)
  }
  # If gaps is set to true report report gaps between cds in lieu of cds
  if(gaps == TRUE)
  {
    message("Calculating UTR Gaps")
    
    # Calcluate gaps between cds regions for each isoform
    UTR <- lapply(UTR, gaps)
    
    # Limit the caluclated gaps to just the chromosomal region of interest
    UTR <- lapply(UTR, function(x) x[seqnames(x) == as.character(seqnames(gr))])
    
    # Limit the calculated gaps to just the strand of interest
    UTR <- lapply(UTR, function(x) x[strand(x) == as.character(strand(gr))])
  }
  
  return(UTR)
}