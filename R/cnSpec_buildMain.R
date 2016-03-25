#' Construct CN cohort plot
#'
#' given a data frame construct a plot to display CN information for a group
#' of samples
#' @name cnSpec_buildMain
#' @param data_frame object of class data frame containing columns chromosome,
#' start, end, cn, sample
#' @param dummy_data Object of class data frame containing columns chromosome,
#' start, end, cn, sample. Used for defining chromosome boundaries
#' @param plot_title character string for title of plot
#' @param CN_low_colour character string specifying low value of colour gradient
#' @param CN_high_colour character string specifying high value of colour
#' gradient
#' @param x_lab_size integer specifying the size of the X label
#' @param y_lab_size integer specifying the size of the Y label
#' @param facet_lab_size integer specifying the size of the faceted labels
#' @param layers Additional layers to be plotted, can be a theme but must be a
#' ggplot layer
#' @param CNscale Character string specifying if copy number calls supplied are
#' relative (i.e.copy neutral == 0) or absolute (i.e. copy neutral ==2). One of 
#' "relative" or "absolute"
#' @return ggplot object
#' @import ggplot2
#' @importFrom scales squish
#' @importFrom scales rescale

cnSpec_buildMain <- function(data_frame, dummy_data, plot_title=NULL,
                             CN_low_colour='#002EB8', CN_high_colour='#A30000',
                             x_lab_size=12, y_lab_size=12, facet_lab_size=10,
                             layers=NULL, CNscale="absolute")
{
    CN_data <- data_frame
    dummy_data <- dummy_data

    # Define Theme of plot
    theme <- theme(strip.text.y=element_text(angle=0, size=facet_lab_size),
                   strip.text.x=element_text(size=facet_lab_size),
                   axis.text.y=element_blank(), axis.ticks.y=element_blank(),
                   axis.text.x=element_blank(), axis.ticks.x=element_blank(),
                   legend.position='right',
                   axis.title.x=element_text(size=x_lab_size, face='bold'),
                   axis.title.y=element_text(size=y_lab_size, face='bold'))

    # Define parameters of plot
    facet <- facet_grid(sample ~ chromosome, scales='free', space='free')
    if(CNscale == "absolute")
    {
        fill_gradient <- scale_fill_gradientn("Copy Number",
                                              colours=c(CN_low_colour,
                                                        'white',
                                                        CN_high_colour),
                                              values=scales::rescale(c(0, 2, 4)),
                                              limits=c(0, 4),
                                              oob=scales::squish)
    } else if(CNscale == "relative") {
        fill_gradient <- scale_fill_gradientn("Copy Number",
                                              colours=c(CN_low_colour,
                                                        "white",
                                                        CN_high_colour),
                                              values=scales::rescale(c(-2, 0, 4)),
                                              limits=c(-2, 4),
                                              oob=scales::squish)
    }

    ylabel <- ylab('Sample')
    xlabel <- xlab('Chromosome')
    title <- ggtitle(plot_title)
    
    # allow the addition of an extra layer
    if(!is.null(layers))
    {
        layers <- layers
    } else {
        layers <- geom_blank()
    }

    # Define main plot using boundaries in dummy data and then plot CN data
    p1 <- ggplot(data=dummy_data,
                 mapping=aes_string(xmin='start',
                                    xmax='end',
                                    ymin=0,
                                    ymax=1)) + geom_rect(alpha=0) +
        scale_x_continuous(expand=c(0,0)) +
        scale_y_continuous(expand=c(0,0))
    
    p1 <- p1 + geom_rect(data=CN_data,
                         mapping=aes_string(xmin='start',
                                            xmax='end',
                                            ymin=0,
                                            ymax=1,
                                            fill='cn'))

    # build the plot
    p1 <- p1 + fill_gradient + ylabel + xlabel + facet + theme_bw() +
        theme + layers

    # if title is supplied plot it
    if(!is.null(plot_title))
    {
        p1 <- p1 + title
    }

    return(p1)
}
