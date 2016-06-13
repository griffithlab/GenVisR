#' construct loh plot
#'
#' given a loh data frame plot points in ggplot
#' @name lohView_buildMain
#' @param x a data frame with columns chromosome, position, n_vaf, t_vaf, sample
#' @param y a data frame with columns chromosome, coordinate for plotting
#' chromosome boundaries
#' @param chr a character string specifying chromosome
#' @param layers additional ggplot2 layers to add
#' @return ggplot2 object
#' @import ggplot2

lohView_buildMain <- function(x, y, chr, layers=NULL)
{
    # Define various parameters of the plot
    dummy_data <- geom_point(data=y, mapping=aes_string(x='coordinate', y=0),
                             alpha=0)
    
    theme <- theme(axis.text.x=element_text(angle=30, hjust=1))

    # y label
    ylabel <- ylab('Variant Allele Fraction')
    
    # provide x label
    xlabel <- xlab('Coordinate')
    
    # allow an extra layer in the plot
    if(!is.null(layers))
    {
        layers <- layers
    } else {
        layers <- geom_blank()
    }
    
    lohpoints <- geom_point(data=x, mapping=aes_string(x='position',
                                                       y='vaf',
                                                       colour='Tissue'))
    
    # build the plot
    tmp <- data.frame(x=0, y=0)
    p1 <- ggplot(data=tmp, aes(x=0)) + lohpoints + ylabel + xlabel +
            theme_bw() + theme + dummy_data + layers
    
    # Plot all chromosomes at once if specified
    if(chr == 'all')
    {
        facet <- facet_wrap(~chromosome, scales='free')
        p1 <- p1 + facet
    }
    
    return(p1)
}