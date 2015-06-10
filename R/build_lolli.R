#' Construct Lolliplot
#' 
#' Construct Lolliplot given gene and mutation data
#' @name build_lolli
#' @param gene_data object of class dataframe giving protien domain and gene information
#' @param length integer specifying the length of the protien in amino acids
#' @param mutation_observed object of class data frame specifying mutations observed in input file
#' @param mutation_cosmic optional object of class data frame specifying mutations observed in cosmic
#' @param fill_value character string specifying the column on which to colour mutation points
#' @param label_column character string specifying the column containing the labels to attach to mutation points
#' @param plot_text_angle numeric value specifying the angle of text to be plotted
#' @param plot_text_size numeric value specifying the size of text to be plotted
#' @param point_size numeric value specigying the size of mutation points
#' @param gene_colour color to shade plotted gene
#' @param sequence_data object of class dataframe giving AA sequence, sidechain, and coord required if plot_sidechain is true
#' @param plot_sidechain boolean specifying whether to plot the AA sidechain instead of domain information
#' @return a ggplot2 object
#' @import ggplot2

build_lolli <- function(gene_data, length, mutation_observed, mutation_cosmic, fill_value, label_column, plot_text_angle, plot_text_size, point_size, gene_colour, sequence_data, plot_sidechain=FALSE)
{
  ###################################################################################
  ####################### Function to make lolliplot type plot ######################
  ###################################################################################

  # build the various features of the plot
  
  # Build gene base either using domain information or AA sidechain information
  if(plot_sidechain == TRUE)
  {
    gene_plot <- geom_rect(data=sequence_data, mapping=aes(xmin=as.numeric(as.character(coord))-1, xmax=as.numeric(as.character(coord)), ymin=-1, ymax=1, fill=sidechain))
    domain_plot <- NULL
  } else {
    gene_plot <- geom_rect(data=gene_data[1,], mapping=aes(xmin=pos_from, xmax=pos_to, ymin=height_min, ymax=height_max), fill='#999999', colour='#000000')
    domain_plot <- geom_rect(data=gene_data[-1,], mapping=aes(xmin=pos_from, xmax=pos_to, ymin=height_min, ymax=height_max, fill=Domain), alpha=0.75, colour='black')
  }

  
  # Build the Observed track
  observed_plot <- geom_point(data=mutation_observed, mapping=aes_string(x="coord_x_dodge", y="coord_y_dodge", colour=fill_value), size=point_size)
  observed_line <- geom_segment(data=mutation_observed, mapping=aes(x=mutation_coord, y=1, xend=coord_x_dodge, yend=1.5))
  observed_line_2 <- geom_segment(data=mutation_observed, mapping=aes(x=coord_x_dodge, y=1.5, xend=coord_x_dodge, yend=coord_y_dodge))
  
  # Miscelaneous features
  title <- ggtitle(gene_data[1,1])
  x_label <- xlab('Amino Acid Position')
  
  # add a theme and guide to the plot
  theme <- theme(legend.position='bottom', legend.direction='vertical', legend.box='horizontal', axis.text.y=element_blank(), axis.ticks.y=element_blank())
  guide <- guides(colour=guide_legend(ncol=2))
  
  
  # construct the plot with or without cosmic track
  if(is.null(mutation_cosmic))
  {	
    y_limits <- ylim(c(-1, max(mutation_observed$coord_y_dodge) + 5))
    y_label <- ylab('Observed')
    p1 <- ggplot() + gene_plot + domain_plot + observed_line_2 + observed_line + observed_plot + x_label + y_label + title + y_limits + theme + guide
    
  } else {
    
    y_limits <- ylim(c(min(mutation_cosmic$coord_y_dodge) - 1, max(mutation_observed$coord_y_dodge) + 1))
    y_label <- ylab('Cosmic v Observed')
    cosmic_plot <- geom_point(data=mutation_cosmic, mapping=aes(x=coord_x_dodge, y=coord_y_dodge), size=point_size)
    cosmic_line <- geom_segment(data=mutation_cosmic, mapping=aes(x=mutation_coord, y=-1, xend=coord_x_dodge, yend=-1.5))
    cosmic_line_2 <- geom_segment(data=mutation_cosmic, mapping=aes(x=coord_x_dodge, y=-1.5, xend=coord_x_dodge, yend=coord_y_dodge))
      
    p1 <- ggplot() + gene_plot + domain_plot + observed_line + observed_line_2 + observed_plot + cosmic_line + cosmic_line_2 + cosmic_plot + x_label + y_label + title + y_limits + theme + guide
  }
  
  # If a label column is specified plot labels
  
  if(!is.null(label_column)) 
  {
    p1 <- p1 + geom_text(data=mutation_observed, mapping=aes(x=coord_x_dodge, y=coord_y_dodge + .06, label=labels), angle=plot_text_angle, size=plot_text_size, vjust=1, hjust=0)
  }
  
  return(p1)
}