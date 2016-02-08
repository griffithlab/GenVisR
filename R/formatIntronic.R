#' format intronic region
#' 
#' given a Granges object obtain and format the intronic region
#' @name formatIntronic
#' @param txdb A TxDb object for a genome
#' @param gr object of class GRanges specifying the ROI
#' @return list of data frames

formatIntronic <- function(txdb, gr)
{
  # Obtain Intronic regions for 5'UTR, 3'UTR, and cds
  threeUTR <- extrthreeUTR(txdb, gr, gaps=T)
  fiveUTR <- extrfiveUTR(txdb, gr, gaps=T)
  cds <- extrcds(txdb, gr, gaps=T)
  
  # Intersect these coordinates giving intronic regions of entire gene
  intron <- mapply(intersect, cds, fiveUTR)
  intron <- mapply(intersect, intron, threeUTR)
  
  # Create Dummy grange object and then intersect isoforms
  #!!!!!!!!!! NOTE: Do this without loop !!!!!!!!!!!!!!!!!
  isoform_intersect <- GRanges(seqnames=c("chr1"), ranges=IRanges(start=c(1), end=c(2)), strand=strand(c("+")))
  for(i in 1:length(intron)-1)
  {
    if(i == 1)
    {
      isoform_intersect <- intersect(intron[[i]], intron[[i + 1]])
    } else {
      isoform_intersect <- intersect(isoform_intersect, intron[[i+1]])
    }
  }
  
  # Calculate GC and perform general formatting
  isoform_intersect <- calcGC(isoform_intersect, genome)
  isoform_intersect <- Granges2dataframe(isoform_intersect)
  isoform_intersect$Type <- c("Intron")
  isoform_intersect$Upper <- c(.5)
  isoform_intersect$Lower <- c(-.5)
  isoform_intersect$Mid <- c(0)
  isoform_intersect$segStart <- c(0)
  isoform_intersect$segEnd <- c(0)
  
  return(isoform_intersect)
}