#' Construct an ideogram
#'
#' Given a data frame with cytogenetic information, construct an ideogram.
#' @name ideoView
#' @param x Object of class data frame with rows representing cytogenetic bands.
#' The data frame must contain the following column names "chrom", "chromStart",
#' "chromEnd", "name", "gieStain"
#' @param chromosome Character string specifying which chromosome from the
#' "chrom" column in the argument supplied to parameter x to plot.
#' @param txtAngle Integer specifying the angle of text labeling cytogenetic
#' bands.
#' @param txtSize Integer specifying the size of text labeling cytogenetic
#' bands.
#' @param plotLayer additional ggplot2 layers for the ideogram
#' @details ideoView is a function designed to plot cytogenetic band
#' inforamtion. Modifications to the graphic object can be made via the
#' `plotLayer` parameter, see vignette for details.
#' @examples
#' # Obtain cytogenetic information for the genome of interest from attached
#' # data set cytoGeno
#' data <- cytoGeno[cytoGeno$genome == 'hg38',]
#'
#' # Call ideoView for chromosome 1
#' ideoView(data, chromosome='chr1', txtSize=4)
#' @return ggplot object
#' @export

ideoView <- function(x, chromosome='chr1', txtAngle=45, txtSize=5,
                     plotLayer=NULL)
{ 
    # Perform quality check
    cytobands <- ideoView_qual(x)

    # format the data obtained from ucsc into something ggplot can understand
    cytobands <- ideoView_formatCytobands(cytobands, chromosome=chromosome)

    # plot the chromosome
    chr_plot <- ideoView_buildMain(cytobands, chromosome=chromosome,
                                   chr_txt_angle=txtAngle,
                                   chr_txt_size=txtSize, layers=plotLayer)

    return(chr_plot)
}
