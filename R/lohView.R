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
#' @return ggplot object
#' @export

lohView <- function(x, y=NULL, genome='hg19')
{
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
     

}