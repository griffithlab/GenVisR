#' produce a coverage plot
#' 
#' produce a coverage plot displaying gene and coverage information
#' @name genCov
#' @param x named list containing data frames with columns stop and cov
#' @param txdb A TxDb object for a genome
#' @param gr A Granges object specifying a region of interest
#' @param genome Object of class BSgenome specifying the genome
#' @param reduce Boolean specifying whether to collapse isoforms in the ROI
#' @param gene_colour character string specifying the colour of the gene to be plotted
#' @param gene_name character string specifying the name of the gene or ROI
#' @param bg_fill character string giving the colour to fill the label
#' @param text_fill character string giving the colour to fill the text
#' @param border character string specifying the colour to fill the border of the label
#' @param size integer specifying the size of the text within the label
#' @param width_ratio vector of length 2 giving the ratio of track labels to plot
#' @param colour character string specifying the color of the data in the plot
#' @param plot_type character string specifying one of line, area for data display
#' @param cores Integer specifying the number of cores to use for processing
#' @param base A vector of log bases to transform the data, corresponding to the elements of transform 
#' @param transform A vector of strings designating what objects to log transform
#' @return ggplot object
#' @export
#' @import GenomicRanges
#' @import plyr

genCov <- function(x, txdb, gr, genome, reduce=F, gene_colour=NULL, gene_name='Gene', bg_fill="black", 
                          text_fill="white", border="black", size=10, width_ratio=c(1, 10), colour="blue",
                          plot_type="line", base=c(10,2,2), transform=c('Intron','CDS','UTR')){
  
  # Obtain a plot for the gene overlapping the Granges object and covert to a named list
  gp_result <- gene_plot(txdb, gr, genome, reduce=reduce, gene_colour=gene_colour,
                    base=base, transform=transform)
  gene <- gp_result$plot
  gene_list <- list()
  gene_list[[gene_name]] <- gene
  
  # Remove entries in coverage data file that are not within the GRanges object specified
  test2 <- function(x, min, max)
  {
    x <- x[-which(x$end <= min | x$end >= max),]
    return(x)
  }
  coverage_data <- lapply(x, test2, min(ranges(gr)), max(ranges(gr)))
  
  # obtain xlimits for gene plot, this is overwritten if transforms
  xlimits <- c(start(gr), end(gr))
  
  # set flag to display x axis labels, overwritten if transforms
  display_x_axis <- TRUE
  
  # perform the intronic transform on the coverage data
  if(!is.null(transform) && length(transform) > 0)
  {
    # Obtain a copy of the master gene file
    master <- gp_result$master
    
    # Format coverage file so that there is a start column, then map coord into transformed intronic space
    test <- function(x)
    {
      x$start <- x$end
      return(x)
    }
    coverage_data <- lapply(coverage_data, test)
    message("Mapping coverage data onto transformed gene-space")
    coverage_data <- lapply(coverage_data, function(x, master) adply(x, 1, map_coord_space, master=master), master=master)

    # Replace original coordinates with transformed coordinates
    for(i in 1:length(coverage_data))
    {
      coverage_data[[i]] <- coverage_data[[i]][,c('trans_start', 'trans_end', 'cov')]
      colnames(coverage_data[[i]]) <- c('start', 'end', 'cov') 
    }
    
    # Obtain x limits for gene plot based on granges object
    start <- cbind(start(gr), start(gr))
    end <- cbind(end(gr), end(gr))
    temp <- as.data.frame(rbind(start, end))
    colnames(temp) <- c('start', 'end')
    temp <- adply(temp, 1, map_coord_space, master=master)
    xlimits <- c(min(temp$trans_start), max(temp$trans_end))
    
    # set flag to not display x axis labels
    display_x_axis <- FALSE
  }
  
  # obtain coverage plots for the data input as a list
  coverage_plot <- lapply(coverage_data, build_coverage, colour=colour, plot_type=plot_type, x_limits=xlimits, display_x_axis=display_x_axis)
  # Combine both gene and coverage plot lists
  merged_data <- c(gene_list, coverage_plot)
  
  # Plot the data on a track
  track_coverage_plot <- plot_track(merged_data, gene_name=gene_name, bg_fill=bg_fill, text_fill=text_fill,
                                    border=border, size=size, axis_align='width', width_ratio=width_ratio, nested_list=T)
  
  return(track_coverage_plot)
}