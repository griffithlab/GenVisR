#' plot LOH data
#'
#' Produce a graphic visualizing Loss of Heterozygosity in a cohort
#' @name lohView
#' @param x object of class data frame with columns 'chromosome', 'position',
#' 'n_vaf', 't_vaf', 'sample'
#' @param y object of class data frame with columns 'chromosome', 'start',
#' 'end' specifying chromosomal boundaries for a genome assembly
#' (required if genome is not specified)
#' @param genome character string specifying the genome assembly from which
#' input data is based
#' @param gender character string specifying whether loh on the X chromosome
#' should be calculated for and plotted
#' @param path character string specifying which directory contains 
#' the sample information stored as datasets with columns "chromosome", 
#' "position", "n_freq", "t_freq", and "sample" (required if x is not specified)
#' @param fileExt character string specifying the file extensions of files
#' within the path specified (required if path argument is specified)
#' @param step integer with the length of divisions (bp) in chromosomes
#' @param window_size integer with the size of the sliding window (bp) to be 
#' applied
#' @param normal integer specifying the normal VAF frequency used in LOH 
#' calculations#' @return ggplot object
#' @param gradient_midpoint object of class numeric specifying the midpoint 
#' for legend's gradient scale
#' @param gradient_low object of class character for hex color code for 
#' gradient's lower values
#' @param gradient_mid object of class character for hex color code for 
#' gradient's middle values
#' @param gradient_high object of class character for hex color code for 
#' gradient's upper values
#' @param plotLayer ggplot theme object specifying parameters for non data
#' elements
#' @return grid object
#' @importFrom gtools mixedsort
#' @examples
#' ## Insert an example once the first part of this is understood
#' @export

lohView <- function(x=NULL, path=NULL, fileExt=NULL, y=NULL, genome='hg19',
                    gender=FALSE, step=500000, window_size=1000000, normal=50,
                    gradient_midpoint=20, gradient_low="#ffffff",
                    gradient_mid="#b2b2ff", gradient_high="#000000",
                    plotLayer=NULL)
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
        x <- lohView_fileGlob(path=path, fileExt=fileExt, step=step, 
                              gender=gender)        
    } 
    
    # Data Quality Check
    input <- lohView_qual(x, y, genome)
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
    
    # Calculate loh via sliding window
    loh <- lohView_slidingWindow(loh_data=x, step, window_size, normal)
    
    # set order of x axis variables in plot
    chromosomes <- gtools::mixedsort(as.character(unique(loh$chromosome)))
    # remove X and/or Y variables
    if (gender == TRUE) {
        chromosomes <- chromosomes[chromosomes != "Y"]
    }
    if (gender == FALSE) {
        chromosomes <- chromosomes[chromosomes != "X" & chromosomes != "Y"]
    }
    loh$chromosome <- factor(loh$chromosome, levels=chromosomes)
    chr_pos$chromosome <- factor(chr_pos$chromosome, levels=chromosomes)
    
    # set order of y axis variables in plot
    samples <- gtools::mixedsort(as.character(unique(loh$sample)))
    loh$sample <- factor(loh$sample, levels=samples)
    
    #build  the plot
    loh_plot <- lohView_buildMain(loh, dummyData=chr_pos,
                                  gradient_midpoint=gradient_midpoint,
                                  gradient_low=gradient_low,
                                  gradient_mid=gradient_mid,
                                  gradient_high=gradient_high,
                                  theme_layer=plotLayer)
    
    #return the final plot
    return(loh_plot)
}
