#' Map regions to transformed space
#' 
#' Reference a master genomic region file to map original positions to a transformed space
#' @name map_coord_space
#' @param master an object of class data frame containing columns start, end, width, type, trans_start, trans_end representing a master genomic region with features from isoforms merged
#' @param coord an object of class data frame containing columns start and end to map to transformed space
#' @return an object of class data frame identical to coord but with extra columns for transformed coord

map_coord_space <- function(master, coord)
{
  #Convert master table into GRanges object to find overlaps in coord
  master_gr <- GRanges(seqnames=c("chr1"), ranges=IRanges(start=master$start, end=master$end))
  mcols(master_gr) <- master[,c('Type', 'trans_start', 'trans_end')]
  
  # Convert input coordinate to GRanges object then find the overlap between master and coord
  coord_gr <- GRanges(seqnames=c("chr1"), ranges=IRanges(start=coord$start, end=coord$end))
  overlap <- subsetByOverlaps(master_gr, coord_gr)
    
  # Convert the grange object specfying the overlap coord back to a data frame
  master_range <- as.data.frame(ranges(overlap))
  master_meta <- as.data.frame(mcols(overlap))
  master_overlap <- cbind(master_range, master_meta)
  master_overlap <- master_overlap[,c('start', 'end', 'width', 'Type', 'trans_start', 'trans_end')]
  master_overlap$width <- master_overlap$trans_end - master_overlap$trans_start
    
  # Convert the grange object containing the coord to transform back to a data frame
  coord_range <- as.data.frame(ranges(coord_gr))
  coord_range <- coord_range[,c('start', 'end', 'width')]
    
  # Check that there is only 1 row in coord_range
  if(nrow(coord_range) != 1)
  {
    stop("Expect data frame coord_range map_coord_space.R to be 1 row")
  }
  
  # Map the original coordinates to the transformed space
  n <- nrow(master_overlap)
  fraction.overlap <- (master_overlap[1,]$end - coord_range$start + 1) / (master_overlap[1,]$end - master_overlap[1,]$start + 1)
  trans_start <- master_overlap[1,]$trans_start + (1 - fraction.overlap) * master_overlap[1,]$width
  fraction.overlap <- (coord_range$end - master_overlap[n,]$start + 1) / (master_overlap[n,]$end - master_overlap[n,]$start + 1)
  trans_end <- master_overlap[n,]$trans_start + fraction.overlap * master_overlap[n,]$width
  
  # Add in the new transformed coordinates to the coord object
  coord$trans_start <- trans_start
  coord$trans_end <- trans_end
  
  return(coord)
}