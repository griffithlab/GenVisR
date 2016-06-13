#' Construct LOH chromosome plot
#'
#' Given a data frame construct a plot to display Loss of Heterozygosity for
#' specific chromosomes.
#' @name lohView
#' @param x object of class data frame with rows representing Heterozygous
#' Germline calls. The data frame must contain columns with the following names 
#' "chromosome", "position", "n_vaf", "t_vaf", "sample".
#' @param y Object of class data frame with rows representing cytogenetic bands
#' for a chromosome. The data frame must contain columns with the following
#' names "chrom", "chromStart", "chromEnd", "name", "gieStain" for plotting the
#' ideogram (optional: see details).
#' @param genome Character string specifying a valid UCSC genome (see details).
#' @param chr Character string specifying which chromosome to plot one of
#' "chr..." or "all"
#' @param ideogram_txtAngle Integer specifying the angle of cytogenetic labels
#' on the ideogram subplot.
#' @param ideogram_txtSize Integer specifying the size of cytogenetic labels on
#' the ideogram subplot.
#' @param plotLayer Valid ggplot2 layer to be added to the copy number plot.
#' @param ideogramLayer Valid ggplot2 layer to be added to the ideogram
#' sub-plot.
#' @param out Character vector specifying the the object to output, one of
#' "data", "grob", or "plot", defaults to "plot" (see returns).
#' @details lohView is able to plot in two modes specified via the `chr`
#' parameter, these modes are single chromosome view in which an ideogram is
#' displayed and genome view where chromosomes are faceted. For the single
#' chromosome view cytogenetic band information is required giving the
#' coordinate, stain, and name of each band. As a convenience GenVisR stores
#' this information for the following genomes "hg19", "hg38", "mm9", "mm10", and
#' "rn5". If the genome assembly supplied to the `genome` parameter is not one
#' of the 5 afore mentioned genome assemblies GenVisR will attempt to query the
#' UCSC MySQL database to retrieve this information. Alternatively the user can
#' manually supply this information as a data frame to the `y` parameter, input
#' to the `y` parameter take precedence of input to `genome`.
#' 
#' A word of caution, users are advised to only use heterozygous germline calls
#' in input to `x`, failure to do so may result in a misleading visual!
#' @examples
#' # Plot loh for chromosome 5
#' lohView(HCC1395_Germline, chr='chr5', genome='hg19', ideogram_txtSize=4)
#' @return One of the following, a list of dataframes containing data to be
#' plotted, a grob object, or a plot.
#' @importFrom stats aggregate
#' @export

lohView <- function(x, y=NULL, genome='hg19', chr='chr1',
                   ideogram_txtAngle=45, ideogram_txtSize=5, plotLayer=NULL,
                   ideogramLayer=NULL, out="plot")
{
    # Perform a basic quality check
    input <- lohView_qual(x, y, genome)
    x <- input[[1]]
    y <- input[[2]]
    
    # Obtain Cytogenetic Band information
    # use y input or query UCSC for the data if it's not preloaded
    preloaded <- c("hg38", "hg19", "mm10", "mm9", "rn5")
    if(is.null(y) && any(genome == preloaded))
    {
        message("genome specified is preloaded, retrieving data...")
        cytobands <- GenVisR::cytoGeno[GenVisR::cytoGeno$genome == genome,]
        cytobands <- cytobands[,-which(colnames(cytobands) == "genome")]
    } else if(is.null(y)) {
        # Obtain data for UCSC genome and extract relevant columns
        memo <- paste0("attempting to query UCSC mySQL database for chromosome",
                       " positions and cytogenetic information")
        message(memo)
        cytobands <- suppressWarnings(multi_cytobandRet(genome=genome))
    } else {
        memo <- paste0("Detected argument supplied to y.. using y for", 
                       "position and cytogenetic information")
        message(memo)
        cytobands <- y
    }
    
    # Create Dummy data and add to x for proper plot dimensions
    fakeStart <- stats::aggregate(data=cytobands, FUN=min, chromStart~chrom)
    colnames(fakeStart) <- c("chromosome", "coordinate")
    fakeEnd <- stats::aggregate(data=cytobands, FUN=max, chromEnd~chrom)
    colnames(fakeEnd) <- c("chromosome", "coordinate")
    dummyData <- rbind(fakeStart, fakeEnd)
    dummyData$chromosome <- as.factor(dummyData$chromosome)
    dummyData <- multi_subsetChr(dummyData, chr)
    
    # Format the main data in x
    x <- reshape2::melt(x, id.vars=c("chromosome", "position", "sample"))
    colnames(x) <- c("chromosome", "position", "sample", "Tissue", "vaf")
    x$Tissue <- sapply(as.character(x$Tissue),
                       function(x) switch(x, "n_vaf"="Normal", "t_vaf"="Tumor"))
    x$Tissue <- as.factor(x$Tissue)
    
    # Plot all chromosomes at once if specified
    if(chr == 'all')
    {
        # plot the graphic
        p1 <- lohView_buildMain(x, dummyData, chr=chr)
    } else {
        # plot chromosome
        chromosome_plot <- ideoView(cytobands, chromosome=chr,
                                    txtAngle=ideogram_txtAngle,
                                    txtSize=ideogram_txtSize,
                                    plotLayer=ideogramLayer)
        
        # if requested plot only selected chromosome
        x <- multi_subsetChr(x, chr)
        
        # build the plot
        LOH_plot <- lohView_buildMain(x, dummyData, chr=chr, layers=plotLayer)
    }
    
    # Decide what to output
    dataOut <- list(main=x, dummyData=dummyData, cytobands=cytobands)
    if(exists("LOH_plot", inherits=FALSE))
    {
        p1 <- multi_align(chromosome_plot, LOH_plot)
        output <- multi_selectOut(data=dataOut, plot=p1, draw=TRUE, out=out)
    } else {
        output <- multi_selectOut(data=dataOut, plot=p1, draw=FALSE, out=out)
    }
    
    return(output)
}