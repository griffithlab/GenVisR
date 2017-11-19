#' Construct Lolliplot
#'
#' Construct Lolliplot given gene and mutation data
#' @name lolliplot_buildMain
#' @param gene_data object of class dataframe giving protien domain and gene
#' information
#' @param length integer specifying the length of the protien in amino acids
#' @param mutation_observed object of class data frame specifying mutations
#' observed in input file
#' @param mutation_observed2 optional object of class data frame specifying
#' additional mutations for bottom track
#' @param fill_value character string specifying the column on which to colour
#' mutation points
#' @param label_column character string specifying the column containing the
#' labels to attach to mutation points
#' @param plot_text_angle numeric value specifying the angle of text to be
#' plotted
#' @param plot_text_size numeric value specifying the size of text to be plotted
#' @param point_size numeric value specigying the size of mutation points
#' @param gene_colour color to shade plotted gene
#' @param sequence_data object of class dataframe giving AA sequence, sidechain,
#' and coord required if plot_sidechain is true
#' @param plot_sidechain boolean specifying whether to plot the AA sidechain
#' instead of domain information
#' @param layers additional ggplot2 layers to plot
#' @param paletteA Character vector specifying colours for gene features
#' @param paletteB Character vector specifying colours for lolli features
#' @return a ggplot2 object
#' @import ggplot2

lolliplot_buildMain <- function(gene_data, length, mutation_observed,
                                mutation_observed2, fill_value, label_column,
                                plot_text_angle, plot_text_size, point_size,
                                gene_colour, sequence_data,
                                plot_sidechain=FALSE,layers=NULL,
                                paletteA=NULL, paletteB=NULL)
{
    # build the various features of the lolliplot

    # Build gene base either using domain information
    # or AA sidechain information
    if(plot_sidechain == TRUE)
    {
        sequence_data$coord_start <-
        as.numeric(as.character(sequence_data$coord)) - 1
        sequence_data$coord_end <- as.numeric(as.character(sequence_data$coord))
        gene_plot <- geom_rect(data=sequence_data,
                               mapping=aes_string(xmin='coord_start',
                                                  xmax='coord_end', ymin=-.1,
                                                  ymax=.1, fill='sidechain'))
        domain_plot <- geom_blank()
    } else {
        gene_plot <- geom_rect(data=gene_data[1,],
                               mapping=aes_string(xmin='pos_from',
                                                  xmax='pos_to',
                                                  ymin='height_min',
                                                  ymax='height_max'),
                               fill='#999999', colour='#000000')
        # Take into account there might not be any domains
        if(nrow(gene_data) == 1)
        {
            domain_plot <- geom_blank()
        } else {
            domain_plot <- geom_rect(data=gene_data[-1,],
                                     mapping=aes_string(xmin='pos_from',
                                                        xmax='pos_to',
                                                        ymin='height_min',
                                                        ymax='height_max',
                                                        fill='Domain'),
                                     alpha=1, colour='black')
        }
    }

    # Build the Observed track
    observed_plot <- geom_point(data=mutation_observed,
                                mapping=aes_string(x='coord_x_dodge',
                                                   y='coord_y_dodge',
                                                   colour=fill_value),
                                size=point_size)
    observed_line <- geom_segment(data=mutation_observed,
                                  mapping=aes_string(x='mutation_coord',
                                                     y=.1, xend='coord_x_dodge',
                                                     yend=.3))
    observed_line_2 <- geom_segment(data=mutation_observed,
                                    mapping=aes_string(x='coord_x_dodge', y=.3,
                                                       xend='coord_x_dodge',
                                                       yend='coord_y_dodge'))

    # Miscelaneous features
    title <- ggtitle(gene_data[1,1])
    x_label <- xlab('Amino Acid Position')

    # add a theme and guide to the plot
    theme <- theme(legend.position='bottom',
                   legend.direction='vertical',
                   legend.box='horizontal',
                   axis.text.y=element_blank(),
                   axis.ticks.y=element_blank(),
                   axis.title.y=element_blank())
    guide <- guides(colour=guide_legend(ncol=2), fill=guide_legend(ncol=2))
    
    # set colours manually if these are specified
    if(!is.null(paletteA))
    {
        gene_features_fill <- scale_fill_manual(values=paletteA)
    } else {
        gene_features_fill <- geom_blank()
    }
    
    if(!is.null(paletteB))
    {
        lolli_features_fill <- scale_colour_manual(values=paletteB)
    } else {
        lolli_features_fill <- geom_blank()
    }

    # construct the plot with or without 2nd observed track
    if(is.null(mutation_observed2))
    {
        y_limits <- ylim(c(-.1, max(mutation_observed$coord_y_dodge) + .1))
        tmp <- data.frame(xmin=0, xmax=0, ymin=0, ymax=0, x=0, y=0, xend=0, yend=0)
        p1 <- ggplot(data=tmp, aes(xmin=0, xmax=0, ymin=0, ymax=0, x=0, y=0, xend=0, yend=0)) +
            gene_plot + domain_plot + observed_line_2 + observed_line +
            observed_plot + x_label + title + y_limits + theme_bw() + theme +
            guide + layers + gene_features_fill + lolli_features_fill
    } else {
        y_limits <- ylim(c(min(mutation_observed2$coord_y_dodge) - .1,
                           max(mutation_observed$coord_y_dodge) + .1))
        if(any(colnames(mutation_observed2) %in% fill_value))
        {
            observed2_plot <- geom_point(data=mutation_observed2,
                                         mapping=aes_string(x='coord_x_dodge',
                                                            y='coord_y_dodge',
                                                            colour=fill_value),
                                         size=point_size)
        } else {
            observed2_plot <- geom_point(data=mutation_observed2,
                                         mapping=aes_string(x='coord_x_dodge',
                                                            y='coord_y_dodge'),
                                         size=point_size)
        }
        observed2_line <- geom_segment(data=mutation_observed2,
                                       mapping=aes_string(x='mutation_coord',
                                                          y=-.1,
                                                          xend='coord_x_dodge',
                                                          yend=-.3))
        observed2_line_2 <- geom_segment(data=mutation_observed2,
                                         mapping=aes_string(x='coord_x_dodge',
                                                            y=-.3,
                                                            xend='coord_x_dodge',
                                                            yend='coord_y_dodge'))
        
        tmp <- data.frame(xmin=0, xmax=0, ymin=0, ymax=0, x=0, y=0, xend=0, yend=0)
        p1 <- ggplot(tmp, aes(xmin=0, xmax=0, ymin=0, ymax=0, x=0, y=0, xend=0, yend=0)) +
            gene_plot + domain_plot + observed_line + observed_line_2 +
            observed_plot + observed2_line + observed2_line_2 + observed2_plot +
            x_label + title + y_limits + theme_bw() + theme + guide + layers +
            gene_features_fill + lolli_features_fill
    }

    # If a label column is specified plot labels
    if(any(colnames(mutation_observed) %in% "labels"))
    {
        mutation_observed$y_label_offset <-
        mutation_observed$coord_y_dodge + .01
        p1 <- p1 + geom_text(data=mutation_observed,
                             mapping=aes_string(x='coord_x_dodge',
                                                y='y_label_offset',
                                                label='labels'),
                             angle=plot_text_angle, size=plot_text_size,
                             vjust=1, hjust=0)
    }
    if(any(colnames(mutation_observed2) %in% "labels"))
    {
        mutation_observed2$y_label_offset <-
        mutation_observed2$coord_y_dodge - .01
        p1 <- p1 + geom_text(data=mutation_observed2,
                             mapping=aes_string(x='coord_x_dodge',
                                                y='y_label_offset',
                                                label='labels'),
                             angle=plot_text_angle, size=plot_text_size,
                             vjust=0, hjust=1)
    }

    return(p1)
}
