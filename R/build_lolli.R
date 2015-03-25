#' Construct Lolliplot
#' 
#' Construct Lolliplot given gene and mutation data
#' @name build_lolli
#' @param gene_data object of class dataframe giving protien domain and gene information
#' @param length integer specifying the length of the protien in amino acids
#' @param mutation_observed object of class data frame specifying mutations observed in input file
#' @param mutation_cosmic optional object of class data frame specifying mutations observed in cosmic
#' @return a ggplot2 object

build_lolli <- function(gene_data, length, mutation_observed, mutation_cosmic, fill_value, label_column, plot_text_angle, plot_text_size, point_size, gene_colour)
{
  ###################################################################################
  ####################### Function to make lolliplot type plot ######################
  ###################################################################################
  library('ggplot2')
  
  # build the various features of the plot
  gene_plot <- geom_rect(data=gene_data[1,], mapping=aes(xmin=pos_from, xmax=pos_to, ymin=height_min, ymax=height_max), fill='#999999', colour='#000000')
  domain_plot <- geom_rect(data=gene_data[-1,], mapping=aes(xmin=pos_from, xmax=pos_to, ymin=height_min, ymax=height_max, fill=Domain), alpha=0.75, colour='black')
  observed_plot <- geom_point(data=mutation_observed, mapping=aes_string(x="coord_x_dodge", y="coord_y_dodge", colour=fill_value), size=point_size)
  observed_rug <- geom_rug(data=mutation_observed, mapping=aes(x=mutation_coord), sides="t", alpha=.25)
  observed_line <- geom_segment(data=mutation_observed, mapping=aes(x=mutation_coord, y=.5, xend=coord_x_dodge, yend=coord_y_dodge))
  x_axis_scale <- scale_x_continuous(limits=c(1, length))
  title <- ggtitle(gene_data[1,1])
  x_label <- xlab('Amino Acid Position')
  
  # add a theme to the plot
  theme <- theme(legend.position='bottom', legend.direction='vertical', legend.box='horizontal', axis.text.y=element_blank(), axis.ticks.y=element_blank())
  
  # construct the plot with or without cosmic track
  if(is.null(mutation_cosmic))
  {	
    y_limits <- ylim(c(-.5, max(mutation_observed$coord_y_dodge) + .5))
    y_label <- ylab('Observed')
    p1 <- ggplot() + gene_plot + domain_plot + observed_line + observed_plot + observed_rug + x_axis_scale + x_label + y_label + title + y_limits + theme
    
  } else {
    
    y_limits <- ylim(c(min(mutation_cosmic$coord_y_dodge) - .5, max(mutation_observed$coord_y_dodge) + .5))
    y_label <- ylab('Cosmic v Observed')
    cosmic_plot <- geom_point(data=mutation_cosmic, mapping=aes(x=coord_x_dodge, y=coord_y_dodge), size=point_size)
    cosmic_line <- geom_segment(data=mutation_cosmic, mapping=aes(x=mutation_coord, y=-.5, xend=coord_x_dodge, yend=coord_y_dodge))
    cosmic_rug <- geom_rug(data=mutation_cosmic, mapping=aes(x=mutation_coord), sides="b", alpha=.25)
    
    p1 <- ggplot() + gene_plot + domain_plot + observed_line + observed_plot + observed_rug + cosmic_line + cosmic_plot + cosmic_rug + x_axis_scale + x_label + y_label + title + y_limits + theme
  }
  
  if(!is.null(label_column)) 
  {
    p1 <- p1 + geom_text(data=mutation_observed, mapping=aes(x=coord_x_dodge, y=coord_y_dodge, label=labels), angle=plot_text_angle, size=plot_text_size, vjust=0, hjust=0)
  }
  
  return(p1)
}