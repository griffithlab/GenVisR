#' retrieve and format CN_cohort plot supplemental data
#'
#' given a genome obtain Start and Stop positions for all chromosomes in the
#' genome
#' @name multi_chrBound
#' @param x data frame containg columns chromosome, start, end
#' @return object of class data frame formatted to internal specifications
#' @importFrom stats aggregate

multi_chrBound <- function(x)
{
    # Check that input has size
    if(nrow(x) < 1)
    {
        memo <- paste0("input has 0 rows, it is possible that the UCSC",
                       " MySQL query has failed")
        stop(memo)
    }
    
    # Extract the columns needed
    data <- x[,c('chrom' ,'chromStart' , 'chromEnd')]

    # Obtain max for each chromosome
    maxChrom <- stats::aggregate(chromEnd ~ chrom, data=data, max)
    maxChrom <- cbind(maxChrom, maxChrom[,2])
    colnames(maxChrom) <- c('chromosome', 'start', 'end')

    # Obtain min for each chromosome
    minChrom <- stats::aggregate(chromStart ~ chrom, data=data, min)
    minChrom <- cbind(minChrom, minChrom[,2])
    colnames(minChrom) <- c('chromosome', 'start', 'end')

    # bind all the data together
    data <- rbind(maxChrom, minChrom)

    return(data)
}
