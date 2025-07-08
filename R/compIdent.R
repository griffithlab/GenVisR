#' Construct identity snp comparison plot
#'
#' Given the bam file path, count the number of reads at the 24 SNP locations
#' @name compIdent
#' @param x data frame with rows representing samples and column names 
#' "sample_name", "bamfile". Columns should correspond to a sample name and a
#' bam file path.
#' @param genome Object of class BSgenome specifying the genome.
#' @param target Object of class data frame containing target locations in
#' 1-base format and containing columns names "chr", "start", "end", "var",
#' "name". Columns should correspond to chromosome, start, end, variant allele, 
#' name of location.
#' @param mainLayer Valid ggplot2 layer for altering the main plot.
#' @param covLayer Valid ggplot2 layer for altering the coverage plot.
#' @param out Character vector specifying the the object to output, one of
#' "data", "grob", or "plot", defaults to "plot" (see returns).
#' @param debug Boolean specifying if test datasets should be used for
#' debugging.
#' @details
#' compIdent is a function designed to comppare samples via variant allele
#' frequencies (VAF) at specific sites. By default these sites correspond to 24
#' identity snps originating from the hg19 assembly however the user can specify
#' alternate sites via the target paramter. To view the 24 identity snp
#' locations use GenVisR::SNPloci.
#' 
#' Samples from the same origin are expected to have similar VAF values however
#' results can skew based on copy number alterations (CNA). The user is expected 
#' to  ensure no CNA occur at the 24 identity snp sites.
#' 
#' For display and debugging purposes a debug parameter is available which will
#' use predefined data instead of reading in bam files. Note that data in the
#' debug parameter is only available at the afore mentioned 24 sites.
#' @return One of the following, a list of dataframes containing data to be
#' plotted, a grob object, or a plot.
#' @examples
#' # Read in BSgenome object (hg19)
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' hg19 <- BSgenome.Hsapiens.UCSC.hg19
#'
#' # Generate plot
#' compIdent(genome=hg19, debug=TRUE)
#' @export

compIdent <- function(x, genome, target=NULL, debug=FALSE, mainLayer=NULL,
                      covLayer=NULL, out="plot")
{

    # Warn if target is null
    if(is.null(target))
    {
        memo <- paste0("Argument not supplied to target, defaulting to",
                       " predefined identity SNPs from hg19 assembly!")
        message(memo)
    }
    
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
        names(bams) <- samplenames

        # Readcount the supplied bam files
        count_tables <- lapply(bams,
                               compIdent_bamRcnt,
                               genome=genome,
                               target=target)
    }
    
    # Format the readcount tables into a form acceptable by compIdent_buildMain
    count_tables <- compIdent_format(count_tables)
    
    # make an sample identity plot
    plot <- compIdent_buildMain(count_tables, mainLayer=mainLayer,
                                covLayer=covLayer)
    
    # Decide what to plot
    output <- multi_selectOut(data=count_tables, plot=plot, draw=TRUE, out=out)
    return(output)
}
