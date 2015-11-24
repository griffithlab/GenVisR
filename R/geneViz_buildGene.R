#' build gene plot
#'
#' given a data frame with gene feature information build the ggplot2 object
#' @name geneViz_buildGene
#' @param data_frame an object of class data frame specifying gene feature
#' information
#' @param display_x_axis Boolean specifying whether to display X axis coordinate
#' values
#' @param x_limits vector specifying x-axis limits of plot
#' @param gene_colour character specifying colour of gene to be plotted
#' @param transcript_name Boolean specifying whether to plot USCS transcript
#' names
#' @param transcript_name_size Integer specifying the size of the transcript
#' name text
#' @param layers additional ggplot2 layers to plot
#' @return ggplot object
#' @import ggplot2

geneViz_buildGene <- function(data_frame, display_x_axis=TRUE, x_limits=NULL,
                              gene_colour=NULL, transcript_name=FALSE,
                              transcript_name_size=4, layers=NULL)
{
    # Define various parameters of plot
    if(is.null(gene_colour))
    {
        gene_features <- geom_rect(data=data_frame,
                                   mapping=aes_string(xmin='start', xmax='end',
                                                      ymin='Upper',
                                                      ymax='Lower', fill='GC'))
    } else {
        gene_features <- geom_rect(data=data_frame,
                                   mapping=aes_string(xmin='start', xmax='end',
                                                      ymin='Upper',
                                                      ymax='Lower'),
                                   fill=gene_colour)
    }

    if(transcript_name == TRUE)
    {
        transcript_data_x <- aggregate(start ~ txname, data=data_frame, min)
        transcript_data_y <- aggregate(Mid ~ txname, data=data_frame, max)
        transcript_data <- merge(transcript_data_x, transcript_data_y,
                                 by="txname")
        transcript_data$label_pos <- transcript_data$start - 1.5
        
        transcript_name <- geom_text(data=transcript_data,
                                     mapping=aes_string(x='label_pos', y='Mid',
                                                        label='txname'),
                                     angle=90, vjust=0,
                                     size=transcript_name_size)
    } else {
        transcript_name <- geom_blank()
    }

    gene_track <- geom_segment(data=data_frame,
                               mapping=aes_string(x='segStart', xend='segEnd',
                                                  y='Mid', yend='Mid'))
    
    if(is.null(x_limits))
    {
        xlimits <- xlim(c(min(data_frame$start),
                          max(data_frame$end)))
    } else {
        xlimits <- xlim(x_limits)
    }

    # Define the theme of the plot
    if(display_x_axis == TRUE)
    {
        theme <- theme(axis.title.x=element_blank(),
                       axis.title.y=element_blank(),
                       axis.text.y=element_blank(),
                       axis.ticks.y=element_blank(), legend.position='top')
    } else {
        theme <- theme(axis.title.x=element_blank(),
                       axis.text.x=element_blank(),
                       axis.ticks.x=element_blank(),
                       axis.title.y=element_blank(),
                       axis.text.y=element_blank(),
                       axis.ticks.y=element_blank(), legend.position='top')
    }

    if(is.null(layers))
    {
        layers <- geom_blank()
    } else {
        layers <- layers
    }

    # Define the main plot
    gene_plot <- ggplot() + gene_track + gene_features + theme_bw() + theme +
        xlimits + transcript_name + layers

    return(gene_plot)
}
