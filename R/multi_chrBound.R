#' retrieve and format CN_cohort plot supplemental data
#'
#' given a genome obtain Start and Stop positions for all chromosomes in the
#' genome
#' @name multi_chrBound
#' @param x data frame containg columns chromosome, start, end
#' @return object of class data frame formatted to internal specifications

multi_chrBound <- function(x)
{
    # Extract the columns needed, should be chromosome start stop in that order
    data <- x[,c(1,2,3)]

    # Obtain max for each chromosome
    maxChrom <- aggregate(chromStart ~ chrom, data=data, max)
    maxChrom <- cbind(maxChrom, maxChrom[,2])
    colnames(maxChrom) <- c('chromosome', 'start', 'end')

    # Obtain max for each chromosome
    minChrom <- aggregate(chromStart ~ chrom, data=data, min)
    minChrom <- cbind(minChrom, minChrom[,2])
    colnames(minChrom) <- c('chromosome', 'start', 'end')

    # bind all the data together
    data <- rbind(maxChrom, minChrom)

    return(data)
}
