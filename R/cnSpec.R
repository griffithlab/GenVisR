#' Construct CN cohort plot
#' 
#' given a data frame construct a plot to display CN information for a group of samples
#' @name cnSpec
#' @param x object of class data frame containing columns "chromosome", "start", "end", "segmean", "sample" consisting of CN segment calls
#' @param y object of class data frame containing columns "chromosome", "start", "end" specifying chromosome boundary coordinates for all chromosomes in a genome (optional)
#' @param genome character string specifying UCSC genome from which data is derived, defaults to "hg19"
#' @param title character string specifying title of plot
#' @param CN_low_colour character string specifying low value of colour gradient to plot
#' @param CN_high_colour character string specifying high value of colour gradient to plot
#' @param x_lab_size integer specifying the size of the x labels on the plot
#' @param y_lab_size integer specifying the size of the y label on the plot
#' @param facet_lab_size integer specifying the size of the faceted labels
#' @param layers Additional ggplot2 layers to plot
#' @return ggplot object
#' @examples
#' cnSpec(LucCNseg, genome="hg19")
#' @export
#' @import plyr
#' @import reshape2
#' @import gtools

cnSpec <- function(x, y=NULL, genome='hg19', title=NULL, CN_low_colour='#002EB8', CN_high_colour='#A30000', x_lab_size=12, y_lab_size=12, facet_lab_size=10, layers=NULL)
{
  # Perform quality check on input data
  data <- cnSpec.qual(x, y, genome)
  x <- data[[1]]
  y <- data[[2]]
  
  # Get dummy data for genome
  # Check to see if y is specified if not check if genome is preloaded
  # else attempt to query UCSC if unsuccessful report an error
  preloaded <- c("hg38", "hg19", "mm10", "mm9", "rn5")
  if(!is.null(y))
  {
    message("detected value in y, reformating...")
    # reformat the input (i.e. put)
    temp <- y
    temp1 <- y
    temp2 <- y
    temp1$end <- temp$start
    temp2$start <- temp$end
    UCSC_Chr_pos <- rbind(temp1, temp2)
  } else if(is.null(y) && any(genome == preloaded)){
    message("genome specified is preloaded, retrieving data...")
    UCSC_Chr_pos <- GenVisR::cytoGeno[GenVisR::cytoGeno$genome == genome,]
    UCSC_Chr_pos <- CN_dummy_data(UCSC_Chr_pos)
  } else {
    # Obtain data for UCSC genome and extract relevant columns
    message("attempting to query UCSC sql database for chromosome positions")
    cyto_data <- suppressWarnings(get_cytobands(genome))
    UCSC_Chr_pos <- CN_dummy_data(cyto_data)
  }
  
  # Check that dummy data has a size, if not report an error
  if(nrow(UCSC_Chr_pos) < 1)
  {
    stop("Could not retrieve chromosome positions from UCSC, please specify y")
  }
  
  # Dcast the input data into a recognizable format
  CN_data <- dcast(x, chromosome + start + end ~ sample, value.var = "segmean")
  
  # Rbind fill the dummy and CN data and remove any chr appendices
  CN_data <- rbind.fill(CN_data, UCSC_Chr_pos)
  CN_data[is.na(CN_data)] <- NA
  CN_data$chromosome <- gsub('chr', '', CN_data$chromosome)
  
  # melt the data for ggplot2 call
  CN_data <- melt(CN_data, id.vars=c('chromosome', 'start', 'end'))
  colnames(CN_data) <- c('chromosome', 'start', 'end', 'sample', 'cn')
  
  # Change the order of chromosomes and samples (natural sort order)
  chromosome_sorted <- unique(mixedsort(CN_data$chromosome))
  CN_data$chromosome <- factor(CN_data$chromosome, levels=chromosome_sorted)
  sample_sorted <- unique(mixedsort(CN_data$sample))
  CN_data$sample <- factor(CN_data$sample, levels=sample_sorted)
  
  # Construct the plot
  p1 <- build.cnSpec(CN_data, plot_title=title, CN_low_colour=CN_low_colour, CN_high_colour=CN_high_colour, x_lab_size=x_lab_size, y_lab_size=y_lab_size, facet_lab_size=facet_lab_size, layers=layers)
  
  return(p1)
}