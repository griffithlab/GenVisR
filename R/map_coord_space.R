#' Map regions to transformed space
#'
#' Reference a master genomic region file to map original positions to a transformed space
#' @name map_coord_space
#' @param master an object of class data frame containing columns start, end, width, type, trans_start, trans_end representing a master genomic region with features from isoforms merged
#' @param coord an object of class data frame containing columns start and end to map to transformed space
#' @return an object of class data frame identical to coord but with extra columns for transformed coord

map_coord_space <- function(master, coord)
{
  # Map the original coordinates to the transformed space
    lt <- coord$start <= master$end
    gt <- coord$end >= master$start

    idx = which(lt & gt)
    n <- length(idx)

    start.row <- master[idx[1],]
    end.row <- master[idx[n],]

    overlap <- (start.row$end - coord$start + 1) / (start.row$end - start.row$start + 1)
    trans_start <- start.row$trans_start + (1 - overlap) * start.row$width
    overlap <- (coord$end - end.row$start + 1) / (end.row$end - end.row$start + 1)
    trans_end <- end.row$trans_start + overlap * end.row$width

  # Add in the new transformed coordinates to the coord object
    coord$trans_start <- trans_start
    coord$trans_end <- trans_end

    return(coord)
}
