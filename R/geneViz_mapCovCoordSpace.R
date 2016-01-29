#' Map coverage track regions to transformed space
#'
#' Reference a master genomic region file to map original positions of a
#' coverage track to a transformed space
#' @name geneViz_mapCovCoordSpace
#' @param master an object of class data frame containing columns start, end,
#' width, type, trans_start, trans_end representing a master genomic region with
#' features from isoforms merged
#' @param cov.coords an object of class data frame containing columns start and
#' end to map to transformed space, with rows demarking single nucleotide
#' coverages
#' @return an object of class data frame identical to coord but with extra
#' columns for transformed coord
#' @importFrom plyr adply

geneViz_mapCovCoordSpace <- function(cov.coords, master)
{ 
    # Check that cov.coords is only single nucleotide coverages
    all(cov.coords$start==cov.coords$end)
    cov.coords$trans_start <- 0
    cov.coords$trans_end <- 0
  
    # Do the following for each row of master:
    transform.coordinates <- function(master, cov.coords){
        lt <- cov.coords$start <= master$end
        gt <- cov.coords$end >= master$start
        
        idx = which(lt & gt)
        
        cov.region.coords <- cov.coords[idx,]
        fraction.start <- (cov.region.coords$start - master$start) / (master$end - master$start + 1)
        fraction.end <- (cov.region.coords$start - master$start + 1) / (master$end - master$start + 1)
        cov.region.coords$trans_start <- (fraction.start * master$width) + master$trans_start
        cov.region.coords$trans_end <- fraction.end * master$width + master$trans_start
        return(cov.region.coords)
    }
    
    x <- plyr::adply(master, 1, transform.coordinates, cov.coords=cov.coords)

    return(x)
}
