#' Construct identity snp comparison plot
#'
#' Given the bam file path, count the number of reads at the 24 SNP locations
#' @name compIdent
#' @param x data frame with column names sample_name, bamfile
#' @param genome Object of class BSgenome specifying the genome
#' @param targetbed Object of class data frame containing target locations in
#' 1-base format and containing columns "chr", "start", "end", "var", "name"
#' @param debug Boolean specifying if test datasets should be used for
#' debugging.
#' @details
#' compIdent is a function designed to compa
#' @return graphical object
#' @examples
#' # Read in BSgenome object (hg19)
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' hg19 <- BSgenome.Hsapiens.UCSC.hg19
#'
#' # Generate plot
#' compIdent(genome=hg19, debug=TRUE)
#' @export

compIdent <- function(x, genome, target=NULL, debug=FALSE)
{
    # Run with the Debug data set if specified
    if(isTRUE(debug))
    {
        bams <- list(normal=GenVisR::HCC1395_N, tumor=GenVisR::HCC1395_T)
        samplenames <- c('normal','tumor')
        count_tables <- lapply(bams,
                               compIdent_bamRcnt,
                               genome=genome,
                               debug=debug,
                               target=target)
    } else { 
        # Grab the bam files and samples from the data frame
        bams <- as.character(x$bamfile)
        samplenames <- as.character(x$sample_name)

        # Readcount the supplied bam files
        count_tables <- lapply(bams,
                               compIdent_bamRcnt,
                               genome=genome,
                               target=target)
        annonymous <- function(x)
        {
            return(x)
            
        }
        lapply(count_tables, annonymous)
    }

    # make an sample identity plot
    plot <- compIdent_buildMain(count_tables, samplenames)

    return(plot)
}
