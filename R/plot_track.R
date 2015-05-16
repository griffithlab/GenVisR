#' plot with named tracka
#' 
#' given a named list of plots, plot them on tracks
#' @name plot_track
#' @param ... named list of plots 
#' @param bg_fill character string giving the colour to fill the label
#' @param text_fill character string giving the colour to fill the text
#' @param border character string specifying the colour to fill the border of the label
#' @param size integer specifying the size of the text within the label
#' @param axis_align character string specifying axis to align plotting space on, one of 'both', 'height', 'width', 'none'
#' @param width_ratio vector of length 2 giving the ratio of track labels to plot
#' @param nested_list boolean specifying whether plots are in a named list nested in another list
#' @return ggplotGrob object
#' @import gridExtra

plot_track <- function(..., bg_fill="black", text_fill="white", border="black", size=10, axis_align='none', width_ratio=c(1, 10), nested_list=F)
{
  # Grab all the tracks/data to be plotted as a named list, check and correct if list is within a list
  if(nested_list==T)
  {
    data <- list(...)
    data <- data[[1]]
  } else {
    data <- list(...)
  }
  
  # Build Track labels and store as a list containing grob objects, then convert to a single grob
  labels <- lapply(names(data), build_track_name, bg_fill=bg_fill, text_fill=text_fill, border=border, size=size)
  label_plot <- do.call(arrangeGrob, lapply(labels, ggplotGrob))
  
  # Convert the plots corresponding to track labels to a single grob
  data_plot <- align_plot(data, axis=axis_align)
  
  # arrange the label and data plot 
  p1 <- grid.arrange(label_plot, data_plot, ncol=2, widths=width_ratio)
  
  return(p1)
}