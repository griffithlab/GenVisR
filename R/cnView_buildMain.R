#' construct CN plot
#'
#' given a CN data frame plot points in ggplot
#' @name cnView_buildMain
#' @param x a data frame with columns chromosome, coordinate, cn, p_value
#' @param y a data frame with columns chromosome, coordinate for plotting
#' boundaries
#' @param z a data frame with columns chromsome, start, end, segmean specifying
#' segments called from copy number (optional)
#' @param chr a character string specifying chromosome
#' @param CNscale Character string specifying if copy number calls supplied are
#' relative (i.e.copy neutral == 0) or absolute (i.e. copy neutral ==2). One of 
#' "relative" or "absolute"
#' @param layers additional ggplot2 layers to add
#' @return ggplot2 object
#' @import ggplot2

cnView_buildMain <- function(x, y, z=NULL, chr, CNscale=FALSE, layers=NULL)
{
    # Define various parameters of the plot
    dummy_data <- geom_point(data=y, mapping=aes_string(x='coordinate', y=2),
                             alpha=0)

    theme <- theme(axis.text.x=element_text(angle=30, hjust=1))
    if(CNscale == "relative")
    {
        # cn fill colors
        shade_cn <- scale_color_gradient2(midpoint=0,
                                          low='#009ACD',
                                          mid='#646082',
                                          high='#C82536',
                                          space='Lab')
        # y label
        ylabel <- ylab('Copy Number Difference')
    } else if(CNscale == "absolute") {
        # cn fill colors
        shade_cn <- scale_color_gradient2(midpoint=2,
                                          low='#009ACD',
                                          mid='#646082',
                                          high='#C82536',
                                          space='Lab')
        # y label
        ylabel <- ylab('Absolute Copy Number')
    } else {
        memo <- paste0("Did not recognize input to CNscale... defaulting to", 
                       "absolute scale, please specify \"relative\"",
                       "if copy neutral calls == 0")
        warning(memo)
    }
    
    # provide x label
    xlabel <- xlab('Coordinate')
    
    # allow an extra layer in the plot
    if(!is.null(layers))
    {
        layers <- layers
    } else {
        layers <- geom_blank()
    }

    # if x contains a p_value column set an alpha for it and plot points
    if(any('p_value' %in% colnames(x)))
    {
        x$transparency <- 1-x$p_value
        cnpoints <- geom_point(data=x, mapping=aes_string(x='coordinate',
                                                          y='cn',
                                                          colour='cn',
                                                          alpha='transparency'))
        transparency <- scale_alpha(guide='none')
    } else {
        cnpoints <- geom_point(data=x, mapping=aes_string(x='coordinate',
                                                          y='cn',
                                                          colour='cn'))
        transparency <- geom_blank()
    }

    # Define segments for main plot
    if(!is.null(z))
    {
        cnseg <- geom_segment(data=z, mapping=aes_string(x='start',
                                                         xend='end',
                                                         y='segmean',
                                                         yend='segmean'),
                              colour='green', size=2)
    } else {
        cnseg <- geom_blank()
    }

    # build the plot
    p1 <- ggplot() + cnpoints + shade_cn + ylabel + xlabel + theme_bw() +
        theme + cnseg + dummy_data + transparency + layers

    if(chr == 'all')
    {
        facet <- facet_wrap(~chromosome, scales='free')
        p1 <- p1 + facet
    }

    return(p1)
}
