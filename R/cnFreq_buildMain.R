#' Construct CN frequency plot
#'
#' given a data frame construct a plot to display proportions of losses and
#' gains across the genome
#' @name cnFreq_buildMain
#' @param x object of class data frame containing columns chromosome,
#'  start, end, gain, and loss
#' @param dummy_data Object of class data frame containing columns chromosome,
#' start, end, cn, sample. Used for defining chromosome boundaries
#' @param plotType character string to determine whether to plot proportions or
#'  frequencies
#' @param plot_title character string for title of plot
#' @param CN_low_colour character string specifying low value of colour gradient
#' @param CN_high_colour character string specifying high value of colour
#'  gradient
#' @param x_lab_size integer specifying the size of the X label
#' @param y_lab_size integer specifying the size of the Y label
#' @param facet_lab_size integer specifying the size of the faceted labels
#' @param plotLayer Additional layer to be plotted, can be a theme but must be a
#'  ggplot layer
#' @return ggplot object
#' @import ggplot2

cnFreq_buildMain <- function(x, plotType, dummy_data, plot_title=NULL,
                             CN_low_colour='#002EB8', CN_high_colour='#A30000',
                             x_lab_size=12, y_lab_size=12, facet_lab_size=10,
                             plotLayer=NULL)
{
    # Transform losses to be negative values for plotting purposes
    x$lossFrequency <- -1*x$lossFrequency
    x$lossProportion <- -1*x$lossProportion

    # Define parameters of plot
    theme <- theme(strip.text.x=element_text(size=facet_lab_size),
                   axis.text.x=element_blank(),
                   axis.ticks.x=element_blank(),
                   legend.position='right',
                   axis.title.x=element_text(size=x_lab_size, face='bold'),
                   axis.title.y=element_text(size=y_lab_size, face='bold'))
    facet <- facet_grid(. ~ chromosome, scales='free', space='free')
    xlabel <- xlab('Chromosomes')
    
    # Choose whether to plot aesthetics for proportion or frequency
    if(grepl("^PROP", plotType, ignore.case=TRUE)){
        ylabel <- ylab("Proportion of Copy Number Gains/Losses")
        ymax <- 1
        x$gain <- x$gainProportion
        x$loss <- x$lossProportion 
    } else if(grepl("^FREQ", plotType, ignore.case=TRUE)){
        ylabel <- ylab("Frequency of Copy Number Gains/Losses")
        ymax <- max(as.numeric(as.character(x$sampleFrequency)), na.rm=TRUE)
        x$gain <- x$gainFrequency
        x$loss <- x$lossFrequency 
    } else {
        memo <- paste0("did not recognize plotType ", plotType,
                       ", please specify one of \"proportion\" or \"frequency\"")
        stop(memo)
    }
    
    # Define the initial plot
    p1 <- ggplot(data=dummy_data,
                 mapping=aes_string(xmin='start',
                                    xmax='end',
                                    ymin=-1*ymax,
                                    ymax=ymax)) + geom_rect(alpha=0) +
        scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0))
    
    # add copy number data 
    p1 <- p1 + geom_rect(data=x, mapping=aes_string(xmin='start',
                                                    xmax='end',
                                                    ymin='loss',
                                                    ymax=0), fill=CN_low_colour)
    p1 <- p1 + geom_rect(data=x, mapping=aes_string(xmin='start',
                                                    xmax='end', 
                                                    ymin=0,
                                                    ymax='gain'), fill=CN_high_colour)
    
    p1 <- p1 + geom_hline(aes(yintercept=0), linetype="dotted")

    # build the plot
    p1 <- p1 + ylabel + xlabel + facet + theme_bw() + theme

    # if there are other layers, add them
    if(!is.null(plotLayer))
    {
        p1 <- p1 + plotLayer
    }
    
    # if title is supplied plot it
    if(!is.null(plot_title))
    {
        p1 <- p1 + ggtitle(plot_title)
    }
    
    return(p1)
}
