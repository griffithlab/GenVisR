#' Construct a region of interest coverage plot
#'
#' Given a list of data frames construct a sequencing coverage view over a
#' region of interest.
#' @name genCov
#' @param x Named list with list elements containing data frames
#' representing samples. Data frame rows should represent read pileups observed
#' in sequencing data. Data frame column names must include "end" and "cov"
#' corresponding to the base end position and coverage of a pileup respectively.
#' Data within data frames must be on the same chromosome as the region of
#' interest, see details!
#' @param txdb Object of class TxDb giving transcription meta data for a genome
#' assembly. See Bioconductor annotation packages.
#' @param gr Object of class GRanges specifying the region of interest and 
#' corresponding to a single gene. See Bioconductor package GRanges.
#' @param genome Object of class BSgenome specifying the genome sequence of
#' interest. See Bioconductor annotation packages.
#' @param reduce Boolean specifying whether to collapse gene isoforms within the
#' region of interest into one representative transcript. Experimental use with
#' caution!
#' @param gene_colour Character string specifying the colour of the gene to be
#' plotted in the gene track.
#' @param gene_plotLayer Valid ggplot2 layer to be added to the gene sub-plot.
#' @param gene_name Character string specifying the name of the gene or region
#' of interest.
#' @param label_bgFill Character string specifying the desired background colour
#' of the track labels.
#' @param label_txtFill Character string specifying the desired text colour of
#' the track labels.
#' @param label_borderFill Character string specifying the desired border colour
#' of the track labels.
#' @param label_txtSize Integer specifying the size of the text within the track
#' labels.
#' @param lab2plot_ratio Numeric vector of length 2 specifying the ratio of
#' track labels to plot space.
#' @param cov_colour Character string specifying the colour of the data in the
#' coverage plots.
#' @param cov_plotType Character string specifying one of "line",
#' "bar" or "point". Changes the ggplot2 geom which constructs the data display.
#' @param cov_plotLayer Valid ggplot2 layer to be added to the coverage
#' sub-plots.
#' @param base Numeric vector of log bases to transform the data
#' corresponding to the elements supplied to the variable transform See details.
#' @param transform Character vector specifying what objects to log transform,
#' accepts "Intron", "CDS", and "UTR" See details.
#' @param gene_labelTranscript Boolean specifying whether to plot the transcript
#' names in the gene plot.
#' @param gene_labelTranscriptSize Integer specifying the size of the transcript
#' name text in the gene plot.
#' @param gene_isoformSel Character vector specifying the names
#' (from the txdb object) of isoforms within the region of interest to display.
#' @param out Character vector specifying the object to output, one of
#' "data", "grob", or "plot", defaults to "plot" (see returns).
#' @param subsample Boolean value specifying whether to reduce the provided
#' coverage data to a subset of approximately 1000 points. Used to generate
#' sparse plots that use less disk space and are faster to render.
#' @details genCov is a function designed construct a series of tracks based on 
#' a TxDb object giving transcript features, and coverage data supplied to
#' parameter `x`. The function will look at a region of interest specified by
#' the argument supplied to gr and plot transcript features and the
#' corresponding coverage information. The argument supplied to `genome` enables
#' gc content within genomic features to be calculated and displayed. The
#' argument supplied to x must contain data on the same chromosome as the region
#' of interest specified in the parameter `gr`!
#' 
#' Typically, introns of a transcript are much larger than exons, while exons are
#' sometimes of greater interest. To address this, genCov will by default
#' scale the x-axis to expand track information according to region type: coding
#' sequence (CDS), untranslated region (UTR), or intron / intergenic (Intron).
#' The amount by which each region is scaled is controlled by the `base` and
#' `transform` arguments. `transform` specifies which regions to scale, and `base`
#' corresponds to the log base transform to apply to those regions. To keep one or
#' more region types from being scaled, omit the corresponding entries from the 
#' `base` and `transform` vectors.
#' 
#' @return One of the following, a list of dataframes containing data to be
#' plotted, a grob object, or a plot.
#' @importFrom GenomicRanges start
#' @importFrom GenomicRanges end
#' @importFrom plyr adply
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
#' # Create Data for input
#' start <- c(89622194:89729524)
#' end <- c(89622195:89729525)
#' chr <- 10
#' cov <- c(rnorm(100000, mean=40), rnorm(7331, mean=10))
#' cov_input_A <- as.data.frame(cbind(chr, start, end, cov))
#'
#' start <- c(89622194:89729524)
#' end <- c(89622195:89729525)
#' chr <- 10
#' cov <- c(rnorm(50000, mean=40), rnorm(7331, mean=10), rnorm(50000, mean=40))
#' cov_input_A <- as.data.frame(cbind(chr, start, end, cov))
#'
#' # Define the data as a list
#' data <- list("Sample A"=cov_input_A)
#' 
#' # Call genCov
#' genCov(data, txdb, gr, genome, gene_labelTranscriptSize=3)
#' @export

genCov <- function(x, txdb, gr, genome, reduce=FALSE, gene_colour=NULL,
                   gene_name='Gene', gene_plotLayer=NULL, label_bgFill="black",
                   label_txtFill="white", label_borderFill="black",
                   label_txtSize=10, lab2plot_ratio=c(1, 10),
                   cov_colour="blue", cov_plotType="point", cov_plotLayer=NULL,
                   base=c(10, 2, 2),
                   transform=c("Intron", "CDS", "UTR"),
                   gene_labelTranscript=TRUE,
                   gene_labelTranscriptSize=4, gene_isoformSel=NULL, out="plot",
                   subsample=FALSE)
{
    # Perform data quality checks
    data <- genCov_qual(x=x, txdb=txdb, gr=gr, genome=genome)
    x <- data$x
    txdb <- data$txdb
    gr <- data$gr
    genome <- data$genome

    # Obtain a plot for the gene overlapping the Granges object and information
    # Used to make the plot
    gp_result <- geneViz(txdb, gr, genome, reduce=reduce,
                         gene_colour=gene_colour, base=base,
                         transform=transform, isoformSel=gene_isoformSel,
                         labelTranscriptSize=gene_labelTranscriptSize,
                         labelTranscript=gene_labelTranscript,
                         plotLayer=gene_plotLayer)
    gene <- gp_result$plot
    gene_list <- list()
    gene_list[[gene_name]] <- gene

    # Remove entries in coverage data file that are not within the GRanges
    # object specified
    test2 <- function(x, min, max)
    {
        idx <- -which(x$end <= min | x$end >= max)
        if(length(idx)){
          x <- x[idx,]
        }
        return(x)
    }
    
    coverage_data <- lapply(x, test2, min(GenomicRanges::ranges(gr)),
                            max(GenomicRanges::ranges(gr)))

    # obtain xlimits for gene plot, this is overwritten if transforms
    xlimits <- c(GenomicRanges::start(gr), GenomicRanges::end(gr))

    # set flag to display x axis labels, overwritten if transforms
    display_x_axis <- TRUE

    # perform the intronic transform on the coverage data
    if(!is.null(transform) && length(transform) > 0)
    {
        # Obtain a copy of the master gene file
        master <- gp_result$master

        # Format coverage file so that there is a start column, then map coord
        # into transformed intronic space
        test <- function(x)
        {
            x$start <- x$end
            return(x)
        }
        
        coverage_data <- lapply(coverage_data, test)
        message("Mapping coverage data onto transformed gene-space")
        coverage_data <- lapply(coverage_data,
                                geneViz_mapCovCoordSpace,
                                master=master)

        # Replace original coordinates with transformed coordinates
        for(i in 1:length(coverage_data))
        {
            coverage_data[[i]] <- coverage_data[[i]][,c('trans_start',
                                                        'trans_end', 'cov')]
            colnames(coverage_data[[i]]) <- c('start', 'end', 'cov')
        }

        # Obtain x limits for gene plot based on granges object
        start <- cbind(GenomicRanges::start(gr), GenomicRanges::start(gr))
        end <- cbind(GenomicRanges::end(gr), GenomicRanges::end(gr))
        temp <- as.data.frame(rbind(start, end))
        colnames(temp) <- c('start', 'end')
        temp <- plyr::adply(temp, 1, geneViz_mapCoordSpace, master=master)
        xlimits <- c(min(temp$trans_start), max(temp$trans_end))

        # set flag to not display x axis labels
        display_x_axis <- FALSE
    }
    
    
    if(subsample){
        transform.coordinates <- function(master_row, x, tx.graph.width){
            lt <- x$start <= master_row$trans_end
            gt <- x$end >= master_row$trans_start
            idx = which(lt & gt)
            x.sub <- x[idx,]
            width <- dim(x.sub)[1]
            adjusted.width <- round(1000 * (master_row$width / tx.graph.width))
            if( width > adjusted.width){
                s <- seq.int(from = 1, to = width, length.out = adjusted.width)
                x.sub <- x.sub[s,]
            }
            return(x.sub)
        }
        perform_subsample <- function(x, master){
            tx.graph.width <- sum(master$width)
            y <- plyr::adply(master, 1, transform.coordinates, x=x, tx.graph.width=tx.graph.width)
            y
        }
        
        coverage_data <- lapply(coverage_data, perform_subsample, master=master)
    }

    # obtain coverage plots for the data input as a list
    coverage_plot <- lapply(coverage_data, genCov_buildCov,
                            colour=cov_colour, plot_type=cov_plotType,
                            x_limits=xlimits, display_x_axis=display_x_axis,
                            layers=cov_plotLayer)
    
    # Combine both gene and coverage plot lists
    merged_data <- c(gene_list, coverage_plot)

    # Plot the data on a track
    track_coverage_plot <- genCov_trackViz(merged_data, gene_name=gene_name,
                                           bgFill=label_bgFill,
                                           textFill=label_txtFill,
                                           border=label_borderFill,
                                           size=label_txtSize,
                                           axis_align='width',
                                           widthRatio=lab2plot_ratio,
                                           list=TRUE)
    # Decide what to output
    output <- multi_selectOut(data=merged_data, plot=track_coverage_plot,
                              draw=TRUE, out=out)
    return(output)
}
