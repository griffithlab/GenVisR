#' plot LOH data
#' 
#' Produce a graphic visualizing Loss of Heterozygosity in a cohort
#' @name lohView
#' @param x object of class data frame with columns 'chromosome', 'position', 
#' 'n_vaf', 't_vaf', 'sample'
#' @param y object of class data frame with columns 'chromosome', 'start', 
#' 'end' specifying chromosomal boundaries for a genome assembly (optional)
#' @param genome character string specifying the genome assembly from which 
#' input data is based
#' @param path character string specifying which directory contains 
#' the sample information stored as datasets with columns "chromosome", 
#' "position", "n_freq", "t_freq", and "sample"
#' @param step integer with the length of divisions (bp) in chromosomes
#' @param window_size integer with the size of the sliding window (bp) to be 
#' applied
#' @param normal integer specifying the normal VAF frequency used in LOH 
#' calculations#' @return ggplot object
#' @param output_path character string specifying the directory the loh heatmap 
#' plot should be exported to
#' @param file_name character string specifying the directory and the name of 
#' the loh heatmap plot produced (pdf extension)
#' @param width integer specfying the width of the pdf produced
#' @param height integer specfying the height of the pdf produced
#' @param 'data' object of class dataframe with loh difference 
#' calculations and column names "window_start", "window_stop", "chromosome", 
#' "sample", and "loh_diff"
#' @param 'gradient_midpoint' object of class numeric specifying the midpoint 
#' for legend's gradient scale
#' @param 'gradient_low' object of class character for hex color code for 
#' gradient's lower values
#' @param 'gradient_mid' object of class character for hex color code for 
#' gradient's middle values
#' @param 'graident_high' object of class character for hex color code for 
#' gradient's upper values
#' @param 'y_facet_lab.size' object of class numeric specifying the facet label 
#' size of y axis
#' @param 'y_facet_lab.angle' object of class numeric specifying the facet label
#' angle of y axis
#' @param 'x_facet_lab.size' object of class numeric specifying the facet label 
#' size of x axis
#' @param 'x_facet_lab.angle' object of class numeric specifying the facet label
#'  angle of x axis
#' @param 'background_color' object of class character specifying the background
#'  color of the facet labels
#' @param 'axis.title.size.x' object of class numeric specifying the x axis 
#' label size 
#' @param 'axis.title.size.y' object of class numeric specifying the y axis 
#' label size 
#' @param 'xlabel' object of class character specifying the title for the x 
#' axis label
#' @param 'ylabel' object of class character specifying the title for the y 
#' axis label
#' @param 'legend.text.size' object of class numeric specifying the text size 
#' for the legend
#' @param 'legend.title.size' object of class numeric specifying the legend 
#' title's size
#' @param 'legend.key.size' object of class numeric specifying the legend 
#' key's size
#' @param 'legend.key.size.unit' object of class cahracter specifying the units 
#' for the legend's parameters
#' @import dplyr, gtools, grid
#' @examples
#' ## Insert an example once the first part of this is understood
#' @export

path <- "/Users/jasonkunisaki/Desktop/Genome Institute/HCC Tasks/HCC_LOH/HCC_LOH_sample/"
file_name <- "/Users/jasonkunisaki/Desktop/final_HCC.pdf"
lohView <- function(x, y=NULL, genome='hg19', path, step=500000, 
                    window_size=1000000, normal=50, file_name, 
                    width=70, height=35, gradient_midpoint=20, 
                    gradient_low="#ffffff", gradient_mid="#b2b2ff", 
                    gradient_high="#000000", 
                    y_facet_lab.size=50, y_facet_lab.angle=0, 
                    x_facet_lab.size=40, x_facet_lab.angle=0, 
                    background_color="Light Grey", axis.title.size.x=70, 
                    axis.title.size.y=70, xlabel="Chromosome", 
                    ylabel="Sample",
                    legend.text.size=45, legend.title.size=50, 
                    legend.key.size.x=4.5, legend.key.size.unit="cm") {
    # Data Quality Check
     input <- lohView_qual(x)
     x <- input[[1]]
     
     # Obtain dummy data for genome
     preloaded <- c('hg38', 'hg19', 'mm10', 'mm9', 'rn5')
     if(!is.null(y))
     {
         message("detected input to y, using supplied positions for chromosome
                 boundaries")
         chr_pos <- y
     } else if(is.null(y) && any(genome == preloaded)) {
         message("genome specified is preloaded, retrieving data...")
         chr_pos <- cytoGeno[cytoGeno$genome == genome,]
         chr_pos <- CN_dummy_data(chr_pos)
         message("developer note: this means you Zach, rename function above")
     } else {
         message("attempting to query UCSC sql database for chromosome
                 positions")
         cyto_data <- suppresswarnings(get_cytobands(genome))
         chr_pos <- CN_dummy_data(cyto_data)
         message("developer note: this means you Zach, rename function above")
     }
     
     # Quality check for dummy data
     if(nrow(chr_pos) < 1)
     {
         stop("Could not retrieve chromosome boundaries from UCSC, please
              specify this information via y")
     }
     loh <- lohView_slidingWindow(path, step, window_size, normal)
     loh$chromosome <- factor(loh$chromosome, 
                                   levels=c(1:22, "X", "Y"))
     samples <- gtools::mixedsort(as.character(unique(loh$sample)))
     loh$sample <- factor(loh$sample, levels=samples)
     loh_plot <- lohView_ggplot2(data=loh, gradient_midpoint, gradient_low, 
                                 gradient_mid, 
                                 gradient_high, y_facet_lab.size, 
                                 y_facet_lab.angle, 
                                 x_facet_lab.size, x_facet_lab.angle, 
                                 background_color, axis.title.size.x, 
                                 axis.title.size.y, xlabel, ylabel, 
                                 legend.text.size, legend.title.size, 
                                 legend.key.size.x, legend.key.size.unit)
     pdf(file_name, width, height)
     print(loh_plot)
     dev.off()
}

