#' Construct identity snp comparison plot
#'
#' Given the bam file path, count the number of reads at the 24 SNP locations
#' @name compIdent
#' @param x data frame with column names sample_name, bamfile
#' @param genome Object of class BSgenome specifying the genome
#' @param debug Boolean specifying if test datasets should be used for
#' debugging.
#' @details
#' TODO
#' @return graphical object
#' @examples
#' #TODO
#' @export

compIdent <- function(x, genome, debug=FALSE)
{
    # Run with the Debug data set if specified
    if(isTRUE(debug))
    {
        bams <- list(normal=GenVisR::HCC1395_N, tumor=GenVisR::HCC1395_T)
        samplenames <- c('normal','tumor')
        count_tables <- lapply(bams,
                               compIdent_bamRcnt,
                               genome=genome,
                               debug=debug)
    } else { 
        # Grab the bam files and samples from the data frame
        bams <- as.character(x$bamfile)
        samplenames <- as.character(x$sample_name)

        # Readcount the supplied bam files
        count_tables <- lapply(bams, compIdent_bamRcnt, genome=genome)
    }

    # make an sample identity plot
    plot <- compIdent_buildMain(count_tables, samplenames)

    return(plot)
}
