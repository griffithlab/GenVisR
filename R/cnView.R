#' Construct copy-number single sample plot
#'
#' Given a data frame construct a plot to display raw copy number calls for a
#' single sample.
#' @name cnView
#' @param x Object of class data frame with rows representing copy number calls
#' from a single sample. The data frame must contain columns with the following
#' names "chromosome", "coordinate", "cn", "p_value".
#' @param y Object of class data frame with rows representing cytogenetic bands
#' for a chromosome. The data frame must contain columns with the following
#' names "chrom", "chromStart", "chromEnd", "name", "gieStain" for plotting the
#' ideogram (optional: see details).
#' @param z Object of class data frame with row representing copy number segment
#' calls. The data frame must contain columns with the following names
#' "chromosome", "start", "end", "segmean" (optional: see details)
#' @param genome character string specifying UCSC genome to use, uneccessary if
#' y is supplied, defaults to "hg19"
#' @param chr character string specifying UCSC chromosome to plot one of
#' "chr..." or "all"
#' @param main.cnDiff Boolean specifying whether values in cn are copy number
#' differences or actual copy numbers
#' @param ideo.chr_txt_angle integer specifying angle of text when plotting
#' chromosome band text
#' @param ideo.chr_txt_size integer specifying size of text when plotting
#' chromosome band text
#' @param main.layers additional ggplot2 layers for the main cn plot
#' @param ideo.layers additional ggplot2 layers for the ideogram
#' @return ggplot object
#' @examples
#' cnView(Luc2CNraw, chr='chr14', genome='hg19', ideo.chr_txt_size=4)
#' @export

cnView <- function(x, y=NULL, z=NULL, genome='hg19', chr='chr1',
                   main.cnDiff=FALSE, ideo.chr_txt_angle=45,
                   ideo.chr_txt_size=5, main.layers=NULL, ideo.layers=NULL)
{
    # Perform a basic quality check
    input <- cnView_qual(x, y, z, genome)
    x <- input[[1]]
    y <- input[[2]]
    z <- input[[3]]

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
    fakeStart <- aggregate(data=cytobands, FUN=min, chromStart~chrom)
    colnames(fakeStart) <- c("chromosome", "coordinate")
    fakeEnd <- aggregate(data=cytobands, FUN=max, chromEnd~chrom)
    colnames(fakeEnd) <- c("chromosome", "coordinate")
    dummyData <- rbind(fakeStart, fakeEnd)
    dummyData$chromosome <- as.factor(dummyData$chromosome)
    dummyData <- cnView_subsetChr(dummyData, chr)

    # Plot all chromosomes at once if specified
    if(chr == 'all')
    {
        p1 <- cnView_buildMain(x, z=z, dummyData, chr=chr)
        return(p1)
    }

    # plot chromosome
    chromosome_plot <- ideoView(cytobands, chromosome=chr,
                                chr_txt_angle=ideo.chr_txt_angle,
                                chr_txt_size=ideo.chr_txt_size,
                                layers=ideo.layers)

    # if requested plot only selected chromosome
    x <- cnView_subsetChr(x, chr)
    if(!is.null(z))
    {
        z <- cnView_subsetChr(z, chr)
    }

    # build the cn plot
    CN_plot <- cnView_buildMain(x, dummyData, z=z, chr=chr, cnDiff=main.cnDiff,
                                layers=main.layers)

    p1 <- cnView_align(chromosome_plot, CN_plot)

    return(grid::grid.draw(p1))
}
