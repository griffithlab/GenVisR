#' Plot LOH data
#'
#' Construct a graphic visualizing Loss of Heterozygosity in a cohort
#' @name lohSpec
#' @param x object of class data frame with rows representing germline calls.
#' The data frame must contain columns with the following names "chromosome",
#' "position", "n_vaf", "t_vaf", "sample". required if path is set to NULL (see
#' details). vaf should range from 0-1.
#' @param y Object of class data frame with rows representing chromosome
#' boundaries for a genome assembly. The data frame must contain columns with
#' the following names "chromosome", "start", "end" (optional: see details).
#' @param genome Character string specifying a valid UCSC genome (see details).
#' @param gender Character vector of length equal to the number of samples,
#' consisting of elements from the set {"M", "F"}. Used to suppress the plotting
#' of allosomes where appropriate.
#' @param path Character string specifying the path to a directory containing
#' germline calls for each sample. Germline calls are expected to be stored as
#' tab-seperated files which contain the following column names "chromosome", 
#' "position", "n_vaf", "t_vaf", and "sample". required if x is set to null
#' (see details).
#' @param fileExt Character string specifying the file extensions of files
#' within the path specified. Required if argument is supplied to path
#' (see details).
#' @param step Integer value specifying the step size (i.e. the number of base
#' pairs to move the window). required when method is set to slide
#' (see details).
#' @param window_size Integer value specifying the size of the window in base
#' pairs in which to calculate the mean Loss of Heterozygosity (see details).
#' @param normal Numeric value within the range 0-1 specifying the expected
#' normal variant allele frequency to be used in Loss of Heterozygosity 
#' calculations. defaults to .50\%
#' @param colourScheme Character vector specifying the colour scale to use from
#' the viridis package. One of "viridis", "magma", "plasma", or "inferno".
#' @param plotLayer Valid ggpot2 layer to be added to the plot.
#' @param method character string specifying the approach to be used for 
#' displaying Loss of Heterozygosity, one of "tile" or "slide" (see details).
#' @param out Character vector specifying the the object to output, one of
#' "data", "grob", or "plot", defaults to "plot" (see returns).
#' @return One of the following, a list of dataframes containing data to be
#' plotted, a grob object, or a plot.
#' @details lohSpec is intended to plot the loss of heterozygosity (LOH) within
#' a sample. As such lohSpec expects input data to contain only LOH calls. Input
#' can be supplied as a single data frame given to the argument x with rows
#' containing germline calls and variables giving the chromosome, position, 
#' normal variant allele frequency, tumor variant allele frequency, and the
#' sample. In lieu of this format a series of .tsv files can be supplied via the 
#' path and fileExt arguments. If this method is choosen samples will be infered
#' from the file names. In both cases columns containing the variant allele
#' frequency for normal and tumor samples should range from 0-1.
#' Two methods exist to calculate and display LOH events. If the method is set
#' to "tile" mean LOH is calculated based on the window_size argument with
#' windows being placed next to each other. If the method is set to slide the
#' widnow will slide and calculate the LOH based on the step parameter.
#' In order to ensure the entire chromosome is plotted lohSpec requries the
#' location of chromosome boundaries for a given genome assembly. As a
#' convenience this information is available for the following genomes "hg19",
#' "hg38", "mm9", "mm10", "rn5" and can be tetrieved by supplying one of the
#' afore mentioned assemblies via the 'genome'paramter. If an argument is
#' supplied to the 'genome' parameter and is unrecognized a query to the UCSC
#' MySQL database will be attempted to obtain the required information. If
#' chromosome boundary locations are unavailable for a given assembly this
#' information can be supplied to the 'y' parameter which has priority over the
#' 'genome' parameter. 
#' @importFrom gtools mixedsort
#' @examples 
#' # plot loh within the example dataset
#' lohSpec(x=HCC1395_Germline)
#' @export

lohSpec <- function(x=NULL, path=NULL, fileExt=NULL, y=NULL, genome='hg19',
                    gender=NULL, step=1000000, window_size=2500000, 
                    normal=.50, colourScheme="inferno", plotLayer=NULL,
                    method="slide", out="plot")
{
    # Grab data if necessary
    if(!is.null(path))
    {
        if(is.null(fileExt))
        {
            memo <- paste0("argument required to variable fileExt if argument ",
                           "is supplied to path")
            stop(memo)
        }
        x <- lohSpec_fileGlob(path=path, fileExt=fileExt, step=step, 
                              window_size=window_size, gender=gender)        
    }
    if (is.null(path)) {
        if (is.null(gender) == FALSE) {
            x <- x[x$chromosome !="Y",]
        }
        if(is.null(gender) == TRUE) {
            x <- x[(x$chromosome != "X" &
                             x$chromosome != "Y"),]
        }
    }

    # Data Quality Check
    input <- lohSpec_qual(x, y, genome)
    x <- input[[1]]
    y <- input[[2]]

    # Obtain dummy data for genome
    preloaded <- c('hg38', 'hg19', 'mm10', 'mm9', 'rn5')
    if(!is.null(y))
    {
        message("detected input to y, using supplied positions for chromosome
                boundaries")
        chr_pos <- y
    } else if(is.null(y) && any(genome == preloaded)) {
        message("genome specified is preloaded, retrieving data...")
        chr_pos <- GenVisR::cytoGeno[GenVisR::cytoGeno$genome == genome,]
        chr_pos <- multi_chrBound(chr_pos)
    } else {
        message("attempting to query UCSC sql database for chromosome
                positions")
        cyto_data <- suppressWarnings(multi_cytobandRet(genome))
        chr_pos <- multi_chrBound(cyto_data)
    }

    # Quality check for dummy data
    if(nrow(chr_pos) < 1)
    {
        memo <- paste0("Could not retrieve chromosome boundaries from",
                       " UCSC, please specify this information via ",
                       "the y paramter")
        stop(memo)
    }
    
    # Produce dataset with loh mean absolute differences 
    if (toupper(method) == 'SLIDE') {
        # Calculate loh via sliding window
        loh <- lohSpec_slidingWindow(loh_data=x, step, window_size, normal)
    }
    else if(toupper(method) == 'TILE') {
        # Calculate loh via tiled window
        ## Insert code
        loh <- lohSpec_tileWindow(loh_data=x, window_size, normal)
    }
    else {
        memo <- paste0("Did not recognize input to parameter method.", 
                       "Please specify one of \"Tile\" or \"Slide\".")
        stop(memo)
    }
    
    # set order of x axis labels in plot
    chromosomes <- gtools::mixedsort(as.character(unique(loh$chromosome)))
    
    # remove X and/or Y chromosomes
    if (is.null(gender) == FALSE) {
        chromosomes <- chromosomes[chromosomes != "Y"]
        chr_pos <- chr_pos[(chr_pos$chromosome != "chrY"),]
        loh <- loh[loh$chromosome != "Y",]
    }
    if (is.null(gender) == TRUE) {
        chromosomes <- chromosomes[chromosomes != "X" & chromosomes != "Y"]
        chr_pos <- chr_pos[(chr_pos$chromosome != "chrX" & 
                                                chr_pos$chromosome != "chrY"),]
        loh <- loh[(loh$chromosome != "Y" & loh$chromosome != "X"),]
    }
    loh$chromosome <- factor(loh$chromosome, levels=chromosomes)
    chr_pos_levels <- gtools::mixedsort(as.character(unique(chr_pos$chromosome)))
    chr_pos$chromosome <- 
        factor(chr_pos$chromosome, levels=chr_pos_levels)
    
    # set order of y axis labels in plot
    samples <- gtools::mixedsort(as.character(unique(loh$sample)))
    loh$sample <- factor(loh$sample, levels=samples)
    
    #build  the plot
    loh_plot <- lohSpec_buildMain(loh, dummyData=chr_pos,
                                  colourScheme=colourScheme,
                                  plotLayer=plotLayer)
    
    # Decide what to output
    output <- multi_selectOut(data=loh, plot=loh_plot, draw=FALSE, out=out)
    return(output)
}
