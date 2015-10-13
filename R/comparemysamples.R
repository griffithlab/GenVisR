#' Compare the identities of multiple samples
#'
#' Given the bam file path, count the number of reads at the 24 SNP locations
#' @name comparemysamples
#' @param x data frame with column names sample_name, bamfile
#' @param genome Object of class BSgenome specifying the genome
#' @return grid object
#' @export

comparemysamples <- function(x, genome)
{
    # Grab the bam files and samples from the data frame
    bams <- as.character(x$bamfile)
    samplenames <- as.character(x$sample_name)

    # Readcount the supplied bam files
    count_tables <- lapply(bams, bamreadcount, genome=genome)

    # make an sample identity plot
    plot <- compare_identities_plot(count_tables,samplenames)

    return(plot)
}
