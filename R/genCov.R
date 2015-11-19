#' produce a coverage plot
#'
#' produce a coverage plot displaying gene and coverage information for a region of interest
#' @name genCov
#' @param x named list containing data frames with columns "end" and "cov" consisting of read pileups at bases of interest
#' @param txdb A TxDb object for a genome
#' @param gr A Granges object specifying a region of interest
#' @param genome Object of class BSgenome specifying the genome
#' @param reduce Boolean specifying whether to collapse isoforms in the ROI
#' @param gene.colour character string specifying the colour of the gene to be plotted
#' @param gene.layers additional ggplot2 layers for the gene plot
#' @param gene_name character string specifying the name of the gene or ROI
#' @param label.bg_fill character string giving the colour to fill the label
#' @param label.text_fill character string giving the colour to fill label text
#' @param label.border character string specifying the colour to fill the border of the label
#' @param label.size integer specifying the size of the text within the label
#' @param label.width_ratio integer vector of length 2 giving the ratio of track labels to plot
#' @param cov.colour character string specifying the colour of the data in the coverage plot
#' @param cov.plot_type character string specifying one of line, area or bar for coverage data display
#' @param cov.layers additional ggplot2 layers for the coverage plot
#' @param base A vector of log bases to transform the data, corresponding to the elements of transform
#' @param transform A vector of strings designating what objects to log transform accepts "Intron", "CDS", and "UTR"
#' @param gene.plot_transcript_name Boolean specifying whether to plot the transcript name
#' @param gene.transcript_name_size Integer specifying the size of the transcript name text
#' @param isoformSel Character vector giving the names (from the txdb object) of isoforms within the region of interest to display 
#' @return ggplot object
#' @examples
#' # Load transcript meta data
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#'
#' # Load BSgenome
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' genome <- BSgenome.Hsapiens.UCSC.hg19
#'
#' # Define a region of interest
#' gr <- GRanges(seqnames=c("chr10"),
#'ranges=IRanges(start=c(89622195), end=c(89729532)), strand=strand(c("+")))
#'
#' # Define preloaded coverage data as named list
#' data <- list("PTEN" = ptenCOV)
#'
#' # Call genCov
#' genCov(data, txdb, gr, genome, gene.transcript_name_size=3)
#' @export
#' @import GenomicRanges
#' @import plyr

genCov <- function(x, txdb, gr, genome, reduce=FALSE, gene.colour=NULL, gene_name='Gene', gene.layers=NULL, label.bg_fill="black",
                   label.text_fill="white", label.border="black", label.size=10, label.width_ratio=c(1, 10), cov.colour="blue",
                   cov.plot_type="line", cov.layers=NULL, base=c(10,2,2), transform=c('Intron','CDS','UTR'),
                   gene.plot_transcript_name=TRUE, gene.transcript_name_size=4, isoformSel=NULL){

  # Perform data quality checks
    data <- genCov.qual(x, txdb, gr, genome)
    x <- data[[1]]
    txdb <- data[[2]]
    gr <- data[[3]]
    genome <- data[[4]]

  # Obtain a plot for the gene overlapping the Granges object and covert to a named list
    gp_result <- geneViz(txdb, gr, genome, reduce=reduce, gene_colour=gene.colour,
                         base=base, transform=transform, isoformSel=isoformSel, transcript_name_size=gene.transcript_name_size,
                         plot_transcript_name=gene.plot_transcript_name, layers=gene.layers)
    gene <- gp_result$plot
    gene_list <- list()
    gene_list[[gene_name]] <- gene

  # Remove entries in coverage data file that are not within the GRanges object specified
    test2 <- function(x, min, max)
    {
        idx <- -which(x$end <= min | x$end >= max)
        if(length(idx)){
          x <- x[idx,]
        }
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
        coverage_data <- lapply(coverage_data, function(x, master) adply(x, 1, map_coord_space, master=master, .progress='text'), master=master)

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
    coverage_plot <- lapply(coverage_data, build.genCov.coverage, colour=cov.colour, plot_type=cov.plot_type, x_limits=xlimits, display_x_axis=display_x_axis, layers=cov.layers)
  # Combine both gene and coverage plot lists
    merged_data <- c(gene_list, coverage_plot)

  # Plot the data on a track
    track_coverage_plot <- trackViz(merged_data, gene_name=gene_name, bgFill=label.bg_fill, textFill=label.text_fill,
                                    border=label.border, size=label.size, axis_align='width', widthRatio=label.width_ratio, list=TRUE)

    return(track_coverage_plot)
}
