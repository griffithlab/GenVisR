#' build chromosome 
#' 
#' given a data frame with cytogenetic band locations plot chromosome in ggplot
#' @name build_chromosome
#' @param data_frame a data frame with columns Chr, Coord, Tumor, Normal, Diff, p_value
#' @param chromosome character string specifying UCSC chromosome to plot one of chr... or all
#' @param chr_txt_angle integer specifying angle of text when plotting band text
#' @param chr_txt_size integer specifying size of text when plotting band text
#' @return ggplot object

build_chromosome <- function(data_frame, chromosome, chr_txt_angle=chr_txt_angle, chr_txt_size=chr_txt_size)
{
  require(ggplot2)
  
  # define theme layer for ggplot
  theme <- theme(axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), legend.position='right')
  # modify ggplot legend defaults
  legend <- scale_fill_brewer()
  
  # define additional layers for ggplot including p arm/q arm text labels, segments connecting label to chromosome
  P_arm_text <- geom_text(data=subset(data_frame, data_frame$alternate == 'top'), mapping=aes(x=band_center, y=text_y, label=name), angle=chr_txt_angle, hjust=0, size=chr_txt_size)
  Q_arm_text <- geom_text(data=subset(data_frame, data_frame$alternate == 'bottom'), mapping=aes(x=band_center, y=text_y, label=name), angle=chr_txt_angle, hjust=1, size=chr_txt_size)
  text_line_seg_p <- geom_segment(data=subset(data_frame, data_frame$alternate == 'top'), mapping=aes(x=band_center, y=text_y, xend=band_center, yend=height_max))
  text_line_seg_q <- geom_segment(data=subset(data_frame, data_frame$alternate == 'bottom'), mapping=aes(x=band_center, y=text_y, xend=band_center, yend=height_min))
  ylabel <- ylab(chromosome)
  
  # define the chromosome main body plot
  chr <- ggplot(data_frame, aes(xmin=chromStart, xmax=chromEnd, ymin=height_min, ymax=height_max)) + geom_rect(aes(fill=gieStain)) + ylim(-1.2, 1.2)
  
  #plot the resulting layers
  chr <- chr + P_arm_text + Q_arm_text + text_line_seg_p + text_line_seg_q + theme + ylabel + legend
  
  return(chr)
}