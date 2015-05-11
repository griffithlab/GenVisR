#' Map regions to transformed space
#' 
#' Reference a master genomic region file to map original positions to a transformed space
#' @name mergeRegions
#' @param master an object of class data frame containing columns start, end, width, type, trans_start, trans_end representing a master genomic region with features from isoforms merged
#' @param coord an object of class data frame containing columns start and end to map to transformed space
#' @return an object of class data frame identical to coord but with extra columns for transformed coord

map_coord_space <- function(master, coord, base=exp(1))
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
    
  # Convert the grange object containing the coord to transform back to a data frame
  coord_range <- as.data.frame(ranges(coord_gr))
  coord_range <- coord_range[,c('start', 'end', 'width')]
    
  # Check that there is only 1 row in both coord_range and master_overlap
  if(nrow(coord_range) != 1 || nrow(master_overlap) != 1)
  {
    stop("Expect data frame coord_range and master_overlap in map_coord_space.R to be 1 row")
  }
    
  # Map the original coordinates to the transformed space
  if( master_overlap$Type != 'Intron' ){
    trans_start <- (coord_range$start - master_overlap$start) + master_overlap$trans_start
    trans_end <- trans_start + coord_range$width
  }else{
    trans_start <- log(coord_range$start - master_overlap$start, base=base) + master_overlap$trans_start
    trans_end <- log(coord_range$start - master_overlap$start + coord_range$width, base=base) + master_overlap$trans_start
  }
  
  # Add in the new transformed coordinates to the coord object
  coord$trans_start <- trans_start
  coord$trans_end <- trans_end
  
  return(coord)
}