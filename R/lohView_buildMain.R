#' Plot LOH data
#' 
#' Build a ggplot2 object displaying calculated LOH data  
#' @name lohView_buildMain
#' @param 'data' object of class dataframe with loh difference 
#' calculations and column names "window_start", "window_stop", "chromosome", 
#' "sample", and "loh_diff"
#' @param 'gradient_midpoint' object of class numeric specifying the midpoint 
#' for legend's gradient scale
#' @param 'gradient_low' object of class character for hex color code for 
#' gradient's lower values
#' @param 'gradient_mid' object of class character for hex color code for 
#' gradient's middle values
#' @param 'graident_high' object of class character for hex color code for 
#' gradient's upper values
#' @param 'y_facet_lab.size' object of class numeric specifying the facet label 
#' size of y axis
#' @param 'y_facet_lab.angle' object of class numeric specifying the facet label
#' angle of y axis
#' @param 'x_facet_lab.size' object of class numeric specifying the facet label 
#' size of x axis
#' @param 'x_facet_lab.angle' object of class numeric specifying the facet label
#'  angle of x axis
#' @param 'background_color' object of class character specifying the background
#'  color of the facet labels
#' @param 'axis.title.size.x' object of class numeric specifying the x axis 
#' label size 
#' @param 'axis.title.size.y' object of class numeric specifying the y axis 
#' label size 
#' @param 'xlabel' object of class character specifying the title for the x 
#' axis label
#' @param 'ylabel' object of class character specifying the title for the y 
#' axis label
#' @param 'legend.text.size' object of class numeric specifying the text size 
#' for the legend
#' @param 'legend.title.size' object of class numeric specifying the legend 
#' title's size
#' @param 'legend.key.size' object of class numeric specifying the legend 
#' key's size
#' @param 'legend.key.size.unit' object of class cahracter specifying the units 
#' for the legend's parameters
#' @return object of class ggplot2
#' @import ggplot2
#' @export

lohView_buildMain <- function(data, gradient_midpoint, gradient_low, gradient_mid, 
                            gradient_high, y_facet_lab.size, y_facet_lab.angle, 
                            x_facet_lab.size, x_facet_lab.angle, 
                            background_color, axis.title.size.x, 
                            axis.title.size.y, xlabel, ylabel, 
                            legend.text.size, legend.title.size, 
                            legend.key.size.x, legend.key.size.unit)
{
    # Define the main plot
    p1 <- ggplot(data, aes(xmin=window_start,
                           xmax=window_stop,
                           ymin=-1,
                           ymax=1))
    
    loh_call <-geom_rect(aes(fill=loh_diff)) 
    
    # Define additional plot parameters
    facet <- facet_grid(sample~chromosome, scales="free", space="free")
    
    x_scale <- scale_x_continuous(expand = c(0, 0))
    y_scale <- scale_y_continuous(expand = c(0,0))
    
    lab_x <- xlab(xlabel)
    lab_y <- ylab(ylabel)
        
    # Define plot aesthetics
    BWscheme <- theme_bw()
    
    plot_theme <- theme(axis.ticks.x = element_blank(),
                        axis.text.x = element_blank(),
                        axis.ticks.y = element_blank(),
                        axis.text.y = element_blank(),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        strip.text.y = element_text(size = y_facet_lab.size,
                                                    angle = y_facet_lab.angle),
                        strip.text.x = element_text(size = x_facet_lab.size,
                                                    angle = x_facet_lab.angle),
                        strip.background = element_rect(fill=background_color),
                        axis.title.x=element_text(size = axis.title.size.x),
                        axis.title.y=element_text(size = axis.title.size.y),
                        legend.text=element_text(size = legend.text.size),
                        legend.title=element_text(size = legend.title.size),
                        legend.key.size=unit(legend.key.size.x,
                                             legend.key.size.unit))
    
    LOHgradient <-scale_fill_gradient2(midpoint = gradient_midpoint,
                                       guide = "colourbar",
                                       high = gradient_high,
                                       mid = gradient_mid,
                                       low = gradient_low)
    
    # Build the plot
    p1 <- p1 + loh_call + facet + x_scale + y_scale + lab_x + lab_y + BWscheme +
        plot_theme + LOHgradient
    
    return(p1)
}
