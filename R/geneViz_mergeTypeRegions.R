#' Create Typed Region Table
#'
#' Create a master region table by merging isoforms
#' @name geneViz_mergeTypeRegions
#' @param type.master A dataframe of all elements of a certain type, such as CDS
#' @return type.master A dataframe of merged elements of a certain type

geneViz_mergeTypeRegions <- function(type.master)
{
    i <- 1
    while (i < nrow(type.master))
    {
        r.start <- type.master[i,2] + 1
        r.stop <- type.master[i+1,1] - 1
        width <- r.stop - r.start + 1
        
        if (width < 1)
        {
            new.start <- type.master[i,1]
            new.stop <- type.master[i+1,2]
            new.width <- new.stop - new.start + 1
            type.master[i,1:3] <- c(new.start,new.stop,new.width)
            type.master <- type.master[-(i+1),]
        }else{
            i <- i + 1
        }
    }
    
    return(type.master)
}
