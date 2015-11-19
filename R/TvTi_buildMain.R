#' build transitions/transversions
#'
#' Given a data frame with columns 'trans_tranv', 'sample', 'Freq', and 'Prop',
#' build a transition/transversion plot
#' @name TvTi_buildMain
#' @param x Object of class data frame containing columns 'trans_tranv',
#' 'sample', 'Freq', and 'Prop'
#' @param y Object of class data frame containing columns 'Prop', 'trans_tranv'
#' for display of expected results
#' @param type Object of class character specifying whether to plot the
#' Proportion or Frequency, one of "Prop"
#' @param label_x_axis boolean specifying wheter to label x axis
#' @param x_axis_text_angle Integer specifying the angle to labels on x_axis
#' @param palette Character vector of length 6 specifying colors for
#' trans/tranv type
#' @param plot_expected Boolean specifying if this is the main TvTi plot or a
#' sub plot for expected values
#' @param tvti.layers Additional ggplot2 layers for the main plot
#' @param expec.layers Additional ggplot2 layers for the expected values plot
#' @param title_x_axis boolean specifying whether to display an x axis title
#' @return GGplot Object
#' @import ggplot2

TvTi_buildMain <- function(x, y=NULL, type='Proportion', label_x_axis=TRUE,
                           x_axis_text_angle=45,
                           palette=c('#D53E4F', '#FC8D59', '#FEE08B', '#E6F598',
                                     '#99D594', '#3288BD'),
                           plot_expected=FALSE, tvti.layers=NULL,
                           expec.layers=NULL, title_x_axis=TRUE)
{

    if(!is.null(y))
    {
        # cumulativley sum the expected values and plot
        y$cumsum <- cumsum(y$Prop)
        expected <- geom_hline(data=y, mapping=aes_string(yintercept='cumsum'),
                               linetype="longdash", size=.5)
    } else {
        expected <- geom_blank()
    }

    # Define various parameters of plot
    if(plot_expected == TRUE)
    {
        bar <- geom_bar(data=x, mapping=aes_string(x=shQuote('Expected'),
                                                   y='Prop',
                                                   fill='trans_tranv'),
                        stat='identity', width=1)
    } else if(toupper(type) == 'PROPORTION') {
        bar <- geom_bar(data=x,
                        mapping=aes_string(x='sample', y='Prop',
                                           fill='trans_tranv'),
                        stat='identity', width=1)
    } else if(toupper(type) == 'FREQUENCY') {
        bar <- geom_bar(data=x,
                        mapping=aes_string(x='sample', y='Freq',
                                           fill='trans_tranv'),
                        stat='identity', width=1)
    }

    ylabel <- ylab(type)
    
    if(title_x_axis == TRUE)
    {
        xlabel <- xlab(paste0("Sample: n=", length(unique(x$sample))))
    } else {
        xlabel <- xlab('')
    }
    
    fill_palette <- scale_fill_manual(name='Transistion/Transversion',
                                      values=palette)
    
    # additional layers to plot?
    if(!is.null(tvti.layers))
    {
        layers <- tvti.layers
    } else {
        layers <- geom_blank()
    }

    # Define theme of plot
    if(plot_expected == TRUE)
    {
        theme <- theme(axis.title.y=element_blank(),
                       axis.text.y=element_blank(),
                       axis.ticks.y=element_blank(),
                       legend.position='none',
                       axis.title.x=element_blank(),
                       axis.text.x=element_text(angle=x_axis_text_angle,
                                                hjust=1, vjust=1))
        if(!is.null(expec.layers))
        {
            layers <- expec.layers
        } else {
            layers <- geom_blank()
        }
        
    } else if(label_x_axis == TRUE) {
        theme <- theme(axis.text.x=element_text(angle=x_axis_text_angle,
                                                hjust=1, vjust=1))
    } else {
        theme <- theme(axis.text.x=element_blank(),
                       axis.ticks.x=element_blank())
    }

    # Define plot
    p1 <- ggplot() + bar + xlabel + ylabel + theme_bw() + theme + fill_palette +
    expected + guides(fill=guide_legend(reverse=TRUE)) + layers

    return(p1)
}
