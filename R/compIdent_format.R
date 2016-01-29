#' Format readcount tables from compIdent
#' 
#' Format readcount tables from compIdent for input into compIdent_buildMain
#' @name compIdent_format
#' @param x Named list of data frames with rows of the data frame corresponding
#' to target locations.
#' @return Formated data frame

compIdent_format <- function(x)
{
    # Function that calculates frequencies from readcounts.
    annonymous_A <- function(v)
    {
        v[,c("A", "C", "G", "T")] <- v[,c("A", "C", "G", "T")]/v[,c("total_reads")]
        return(v)
    }
    x <- lapply(x, annonymous_A)
    
    # Define VAF based upon variant allele 
    annonymous_B <- function(x)
    {
        vaf <- apply(x,1,function(y){return(as.numeric(y[as.character(y["var"])]))})
        x$vaf <- vaf
        return(x)
    }
    x <- lapply(x, annonymous_B)
    
    # add column specifying the name of the lists and convert to single data
    # frame
    x <- Map(cbind, x, sample=names(x))
    x <- do.call(rbind, x)
    
    return(x)
}