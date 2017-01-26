#' Construct copy-number frequency plot
#'
#' Given a data frame construct a plot to display copy number changes across the
#' genome for a group of samples.
#' @name cnFreq
#' @param x Object of class data frame with rows representing genomic segments.
#' The data frame must contain columns with the following names "chromosome",
#' "start", "end", "segmean", and "sample".
#' @param CN_low_cutoff Numeric value representing the point at or below which
#' copy number alterations are considered losses. Only used if x represents CN
#' values.
#' @param CN_high_cutoff Numeric value representing the point at or above which
#' copy number alterations are considered gains. Only used if x represents CN
#' values.
#' @param plot_title Character string specifying the title to display on the
#' plot.
#' @param CN_Loss_colour Character string specifying the colour value for copy
#' number losses.
#' @param CN_Gain_colour Character string specifying the colour value for copy
#' number gains.
#' @param x_title_size Integer specifying the size of the x-axis title.
#' @param y_title_size Integer specifying the size of the y-axis title.
#' @param facet_lab_size Integer specifying the size of the faceted labels
#'  plotted.
#' @param plotLayer Valid ggplot2 layer to be added to the plot.
#' @param plotType Character string specifying the type of values to plot. 
#' One of "proportion" or "frequency"
#' @param genome Character string specifying a valid UCSC genome (see details).
#' @param out Character vector specifying the the object to output, one of
#' "data", "grob", or "plot", defaults to "plot" (see returns).
#' @details cnFreq requires the location of chromosome boundaries for a given
#' genome assembly in order to ensure the entire chromosome space is plotted.
#' As a convenience this information is available to cnSpec for
#' the following genomes "hg19", "hg38", "mm9", "mm10", "rn5" and can be
#' retrieved by supplying one of the afore mentioned assemblies via the `genome`
#' parameter. If a genome assembly is supplied to the `genome` parameter and is
#' unrecognized cnSpec will attempt to query the UCSC MySQL database for the
#' required information. If genomic segments are not identical across all samples
#' the algorithm will attempt to perform a disjoin operation
#' splitting existing segments such that there are no overlaps. The `plotLayer`
#' parameter can be used to add an additional layer to the ggplot2 graphic
#' (see vignette).
#' @return One of the following, a dataframe containing data to be
#' plotted, a grob object, or a plot.
#' @importFrom gtools mixedsort
#' @examples
#' # plot on internal GenVisR dataset
#' cnFreq(LucCNseg)
#' @export

cnFreq <- function(x, CN_low_cutoff=1.5, CN_high_cutoff=2.5, plot_title=NULL,
                   CN_Loss_colour='#002EB8', CN_Gain_colour='#A30000',
                   x_title_size=12, y_title_size=12, facet_lab_size=10,
                   plotLayer=NULL, plotType="proportion", genome="hg19", out="plot")
{
    # Perform quality check on input data
    x <- cnFreq_qual(x)
    samples <- unique(x$sample)

    # Calculate a columns of Observed CN gains/losses/and obs samples in the
    # cohort for each segment
    gainFreq <- function(x){length(x[x >= CN_high_cutoff])}
    gainFrequency <- aggregate(segmean ~ chromosome + start + end, data=x, gainFreq)$segmean
    
    lossFreq <- function(x){length(x[x <= CN_low_cutoff])}
    lossFrequency <- aggregate(segmean ~ chromosome + start + end, data=x, lossFreq)$segmean
    
    x <- aggregate(segmean ~ chromosome + start + end, data=x, length)
    colnames(x)[which(colnames(x) %in% "segmean")] <- "sampleFrequency"
    x$gainFrequency <- gainFrequency
    x$lossFrequency <- lossFrequency
    
    # Calculate the proportion
    x$gainProportion <- x$gainFrequency/length(samples)
    x$lossProportion <- x$lossFrequency/length(samples)
    
    # get the dummy data for plot boundaries
    preloaded <- c("hg38", "hg19", "mm10", "mm9", "rn5")
    if(any(genome == preloaded)){
        message("genome specified is preloaded, retrieving data...")
        UCSC_Chr_pos <- GenVisR::cytoGeno[GenVisR::cytoGeno$genome == genome,]
        UCSC_Chr_pos <- multi_chrBound(UCSC_Chr_pos)
    } else {
        # Obtain data for UCSC genome and extract relevant columns
        memo <- paste0("attempting to query UCSC mySQL database for chromosome",
                       " positions")
        message(memo)
        cyto_data <- suppressWarnings(multi_cytobandRet(genome))
        UCSC_Chr_pos <- multi_chrBound(cyto_data)        
    }     
    
    # check that the dummy data has a size
    if(nrow(UCSC_Chr_pos) < 1)
    {
        memo <- paste0("did not recognize genome ", genome,
                       ", plotting provided data and ignoring chromosome ",
                       "boundaries! Output could be decieving!")
        warning(memo)
    }
    
    dummy_data <- lapply(unique(x$sample),
                         function(sample, chr_pos) cbind(chr_pos, sample),
                         UCSC_Chr_pos)
    dummy_data <- do.call("rbind", dummy_data)
    chr_order <- gtools::mixedsort(unique(dummy_data$chromosome))
    dummy_data$chromosome <- factor(dummy_data$chromosome, levels=chr_order)
    
    # build the plot
    p1 <- cnFreq_buildMain(x, plotType, dummy_data = dummy_data, plot_title=plot_title,
                           CN_low_colour=CN_Loss_colour,
                           CN_high_colour=CN_Gain_colour,
                           x_lab_size=x_title_size,
                           y_lab_size=y_title_size,
                           facet_lab_size=facet_lab_size,
                           plotLayer=plotLayer)
    
    # Decide what to output
    output <- multi_selectOut(data=list("data"=x), plot=p1, out=out)
    return(output)
}
