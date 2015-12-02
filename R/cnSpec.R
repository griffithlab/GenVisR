#' Construct CN cohort plot
#'
#' given a data frame construct a plot to display CN information for a group of
#' samples
#' @name cnSpec
#' @param x object of class data frame containing columns "chromosome", "start",
#' "end", "segmean", "sample" consisting of CN segment calls
#' @param y object of class data frame containing columns "chromosome", "start",
#' "end" specifying chromosome boundary coordinates for all chromosomes in a
#' genome (optional)
#' @param genome character string specifying UCSC genome from which data is
#' derived, defaults to "hg19"
#' @param title character string specifying title of plot
#' @param CN_low_colour character string specifying low value of colour gradient
#' to plot
#' @param CN_high_colour character string specifying high value of colour
#' gradient to plot
#' @param x_lab_size integer specifying the size of the x labels on the plot
#' @param y_lab_size integer specifying the size of the y label on the plot
#' @param facet_lab_size integer specifying the size of the faceted labels
#' @param layers Additional ggplot2 layers to plot
#' @return ggplot object
#' @importFrom reshape2 dcast
#' @importFrom reshape2 melt
#' @examples
#' cnSpec(LucCNseg, genome="hg19")
#' @export

cnSpec <- function(x, y=NULL, genome='hg19', title=NULL,
                   CN_low_colour='#002EB8', CN_high_colour='#A30000',
                   x_lab_size=12, y_lab_size=12, facet_lab_size=10, layers=NULL)
{
    # Perform quality check on input data
    data <- cnSpec_qual(x, y, genome)
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

    # Dcast the input data into a recognizable format
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

    # Construct the plot
    p1 <- cnSpec_buildMain(CN_data, plot_title=title,
                           CN_low_colour=CN_low_colour,
                           CN_high_colour=CN_high_colour,
                           x_lab_size=x_lab_size,
                           y_lab_size=y_lab_size,
                           facet_lab_size=facet_lab_size,
                           layers=layers)

    return(p1)
}
