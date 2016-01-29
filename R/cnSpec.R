#' Construct copy-number cohort plot
#'
#' Given a data frame construct a plot to display copy-number calls for a cohort
#' of samples.
#' @name cnSpec
#' @param x Object of class data frame with rows representing copy-number
#' segment calls. The data frame must contain columns with the following names
#' "chromosome", "start", "end", "segmean", "sample".
#' @param y Object of class data frame with rows representing chromosome
#' boundaries for a genome assembly. The data frame must contain columns with
#' the following names "chromosome", "start", "end" (optional: see details).
#' @param genome Character string specifying a valid UCSC genome (see details).
#' @param plot_title Character string specifying title to display on the plot.
#' @param CN_Loss_colour Character string specifying the colour value of copy
#' number losses.
#' @param CN_Gain_colour Character string specifying the colour value of copy
#' number gains.
#' @param x_title_size Integer specifying the size of the x-axis title.
#' @param y_title_size Integer specifying the size of the y-axis title.
#' @param facet_lab_size Integer specifying the size of the faceted labels
#' plotted.
#' @param plotLayer Valid ggplot2 layer to be added to the plot.
#' @param dataOut Boolean specifying whether to output the data to be passed to
#' ggplot instead of plotting.
#' @param CNscale Character string specifying if copy number calls supplied are
#' relative (i.e.copy neutral == 0) or absolute (i.e. copy neutral ==2). One of 
#' "relative" or "absolute"
#' @details cnSpec requires the location of chromosome boundaries for a given
#' genome assembly in order to ensure the entire chromosome space is plotted.
#' As a convenience this information is available to cnSpec for
#' the following genomes "hg19", "hg38", "mm9", "mm10", "rn5" and can be
#' retrieved by supplying one of the afore mentioned assemblies via the `genome`
#' parameter. If a genome assembly is supplied to the `genome` parameter and is
#' unrecognized cnSpec will attempt to query the UCSC MySQL database for the
#' required information. If chromosome boundary locations are unavailable for
#' a given assembly or if it is desireable to plot a specific region
#' encapsulating the copy number data these boundaries can be supplied to the
#' `y` paramter which has priority of the parameter `genome`.
#' 
#' The `plotLayer` parameter can be used to add an additional layer to the
#' ggplot2 graphic (see vignette).
#' @return ggplot object
#' @importFrom reshape2 dcast
#' @importFrom reshape2 melt
#' @examples
#' cnSpec(LucCNseg, genome="hg19")
#' @export

cnSpec <- function(x, y=NULL, genome='hg19', plot_title=NULL,
                   CN_Loss_colour='#002EB8', CN_Gain_colour='#A30000',
                   x_title_size=12, y_title_size=12, facet_lab_size=10,
                   plotLayer=NULL, dataOut=FALSE, CNscale="absolute")
{
    # Perform quality check on input data
    data <- cnSpec_qual(x, y, genome, CNscale=CNscale)
    x <- data[[1]]
    y <- data[[2]]

    # Get dummy data for genome 
    # (this is used to set chromosome boundaries in plot)
    
    # Check to see if y is specified if not check if genome is preloaded
    # else attempt to query UCSC, if unsuccessful report an error
    preloaded <- c("hg38", "hg19", "mm10", "mm9", "rn5")
    if(!is.null(y))
    {
        message("detected value in y, reformating...")
        # reformat the input
        temp <- y
        temp1 <- y
        temp2 <- y
        temp1$end <- temp$start
        temp2$start <- temp$end
        UCSC_Chr_pos <- rbind(temp1, temp2)
        
    } else if(is.null(y) && any(genome == preloaded)){
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

    # Check that dummy data has a size, if not report an error
    if(nrow(UCSC_Chr_pos) < 1)
    {
        memo <- paste0("query to UCSC unsuccessful, please supply data frame",
                       " to argument y")
        stop(memo)
    }

    # Dcast the input data into a recognizable format this will set as NA
    # segment calls in one sample but not the other
    CN_data <- reshape2::dcast(x, chromosome + start + end ~ sample,
                               value.var = "segmean")

    # Rbind fill the dummy and CN data and remove any chr appendices
    CN_data <- plyr::rbind.fill(CN_data, UCSC_Chr_pos)
    CN_data[is.na(CN_data)] <- NA
    CN_data$chromosome <- gsub('chr', '', CN_data$chromosome)

    # melt the data for ggplot2 call
    CN_data <- reshape2::melt(CN_data, id.vars=c('chromosome', 'start', 'end'))
    colnames(CN_data) <- c('chromosome', 'start', 'end', 'sample', 'cn')

    # Change the order of chromosomes and samples (natural sort order)
    chromosome_sorted <- as.vector(unique(CN_data$chromosome))
    chromosome_sorted <- gtools::mixedsort(chromosome_sorted)
    CN_data$chromosome <- factor(CN_data$chromosome, levels=chromosome_sorted)
    sample_sorted <- as.vector(unique(CN_data$sample))
    sample_sorted <- gtools::mixedsort(sample_sorted)
    CN_data$sample <- factor(CN_data$sample, levels=sample_sorted)
    
    # if requested output data instead of plotting
    if(isTRUE(dataOut))
    {
        return(list(cnData=CN_data))
    }
    
    # Construct the plot
    p1 <- cnSpec_buildMain(CN_data, plot_title=plot_title,
                           CN_low_colour=CN_Loss_colour,
                           CN_high_colour=CN_Gain_colour,
                           x_lab_size=x_title_size,
                           y_lab_size=y_title_size,
                           facet_lab_size=facet_lab_size,
                           layers=plotLayer, CNscale=CNscale)

    return(p1)
}
