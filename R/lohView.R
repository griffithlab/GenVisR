#' Plot LOH data
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
#' @param gender vector of length equal to the number of samples, consisting of 
#' elements from the set {"M", "F"}
#' @param path character string specifying which directory contains 
#' the sample information stored as datasets with columns "chromosome", 
#' "position", "n_vaf", "t_vaf", and "sample" (required if x is not specified)
#' @param fileExt character string specifying the file extensions of files
#' within the path specified (required if path argument is specified)
#' @param step integer with the length of divisions (bp) in chromosomes (only 
#' required when tile=FALSE
#' @param window_size integer with the size of the sliding window (bp) to be 
#' applied
#' @param normal integer specifying the normal VAF frequency used in LOH 
#' calculations 
#' @param gradient_midpoint object of class numeric specifying the midpoint 
#' for legend's gradient scale
#' @param gradient_low object of class character for hex color code for 
#' gradient's lower values
#' @param gradient_mid object of class character for hex color code for 
#' gradient's middle values
#' @param gradient_high object of class character for hex color code for 
#' gradient's upper values
#' @param theme_layer ggplot theme object specifying parameters for non data
#' elements
#' @param method Character string specifying the approach to be used for
#' displaying Loss of Heterozygosity, one of "tile" or "slide".
#' @return grid object
#' @importFrom gtools mixedsort
#' @examples
#' ## Insert an example once the first part of this is understood

lohView <- function(x=NULL, path=NULL, fileExt=NULL, y=NULL, genome='hg19',
                    gender=NULL, step=500000, window_size=1000000, 
                    normal=50, gradient_midpoint=20, gradient_low="#ffffff",
                    gradient_mid="#b2b2ff", gradient_high="#000000",
                    theme_layer=NULL, method="tile")
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
    
    # Produce dataset with loh mean absolute differences 
    if (toupper(method) == 'SLIDE') {
        # Calculate loh via sliding window
        loh <- lohView_slidingWindow(loh_data=x, step, window_size, normal)
    }
    else if(toupper(method) == 'TILE') {
        # Calculate loh via tiled window
        ## Insert code
        loh <- lohView_tileWindow(loh_data=x, window_size, normal)
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
    }
    if (is.null(gender) == TRUE) {
        chromosomes <- chromosomes[chromosomes != "X" & chromosomes != "Y"]
        chr_pos <- chr_pos[(chr_pos$chromosome != "chrX" & 
                                                chr_pos$chromosome != "chrY"),]
    }
    loh$chromosome <- factor(loh$chromosome, levels=chromosomes)
    chr_pos$chromosome <- factor(chr_pos$chromosome, levels=chromosomes)
    
    # set order of y axis labels in plot
    samples <- gtools::mixedsort(as.character(unique(loh$sample)))
    loh$sample <- factor(loh$sample, levels=samples)
    
    #build  the plot
    loh_plot <- lohView_buildMain(loh, dummyData=chr_pos,
                                  gradient_midpoint=gradient_midpoint,
                                  gradient_low=gradient_low,
                                  gradient_mid=gradient_mid,
                                  gradient_high=gradient_high,
                                  theme_layer=theme_layer)
    
    #return the final plot
    return(loh_plot)
}
