#' reformat cytogenetic band data frame
#'
#' given a data frame of cytogenetic bands, format it for ggplot call
#' @name ideoView_formatCytobands
#' @param data_frame a data frame retrieved from UCSC giving cytogenetic
#' information
#' @param chromosome character string specifying chromosome of interest from
#' UCSC data frame
#' @return object of class data frame

ideoView_formatCytobands <- function(data_frame, chromosome)
{
    # set height of chromosome for plotting
    data_frame$height_min <- -.5
    data_frame$height_max <- .5

    # subset bands retrieve pulling just the chromosome requested, also creating
    # an alternating vector to map text to
    data_frame <- data_frame[data_frame$chrom == chromosome,]
    alternate <- c("top", "bottom")
    alternate <- rep_len(alternate, nrow(data_frame))
    data_frame$alternate <- alternate

    # Get median of band for geom_text
    band_median <- apply(data_frame[,c("chromStart", "chromEnd")],
                         1,
                         stats::median)
    
    data_frame$band_center <- band_median

    # set up y positions for geom_text call
    text_y <- rep_len(c(.7,-.7), nrow(data_frame))
    data_frame$text_y <- text_y

    # add p/q arm distinction
    data_frame$arm <- NA
    Parm_boolean <- grepl("^p", data_frame[,c('name')])
    data_frame[Parm_boolean, c('arm')] <- "p"
    Qarm_boolean <- grepl("^q", data_frame[,c('name')])
    data_frame[Qarm_boolean, c('arm')] <- "q"

    return(data_frame)
}
