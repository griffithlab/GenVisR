#' Plot LOH data
#' 
#' Build a ggplot2 object displaying calculated LOH data  
#' @name lohView_buildMain
#' @param x object of class dataframe with loh difference 
#' calculations and column names "window_start", "window_stop", "chromosome", 
#' "sample", and "loh_diff"
#' @param dummyData object of class dataframe with column names "chromosome",
#' "start", "end" specifying chromosome boundaries
#' @param gradient_midpoint object of class numeric specifying the midpoint 
#' for legend's gradient scale
#' @param gradient_low object of class character for hex color code for 
#' gradient's lower values
#' @param gradient_mid object of class character for hex color code for 
#' gradient's middle values
#' @param graident_high object of class character for hex color code for 
#' gradient's upper values
#' @param theme_layer ggplot theme object specifying parameters for non data
#' elements
#' for the legend's parameters
#' @return object of class ggplot2
#' @import ggplot2

lohView_buildMain <- function(x, dummyData,
                              gradient_midpoint=gradient_midpoint,
                              gradient_low=gradient_low,
                              gradient_mid=gradient_mid,
                              gradient_high=gradient_high, xlabel=xlabel,
                              ylabel=ylabel, theme_layer=theme_layer)
{
    # define dummy data which will be chromosome boundaries, these are plotted
    # but are transparent and will not appear in the plot
    dummy_data <- geom_rect(data=dummyData, aes_string(xmin='start', xmax='end',
                                                       ymin=-1, ymax=1),
                            alpha=0)
    # Define the main plot
    data <- geom_rect(data=x, aes_string(xmin='window_start',
                                         xmax='window_stop',
                                         ymin=-1,
                                         ymax=1, fill='loh_diff'))

    
    # Define additional plot parameters
    facet <- facet_grid(sample ~ chromosome, scales="free", space="free")
    
    x_scale <- scale_x_continuous(expand = c(0, 0))
    y_scale <- scale_y_continuous(expand = c(0,0))
    
    lab_x <- xlab("Chromosome")
    lab_y <- ylab("Sample")
        
    # Define plot aesthetics
    BWscheme <- theme_bw()
    
    if(!is.null(theme_layer))
    {
        plot_theme <- theme_layer
    } else {
        plot_theme <- theme(axis.ticks.x=element_blank(),
                            axis.text.x=element_blank(),
                            axis.ticks.y=element_blank(),
                            axis.text.y=element_blank(),
                            panel.grid.major=element_blank(),
                            panel.grid.minor=element_blank())
    }

    
    LOHgradient <-scale_fill_gradient2(midpoint = gradient_midpoint,
                                       guide = "colourbar",
                                       high = gradient_high,
                                       mid = gradient_mid,
                                       low = gradient_low,
                                       space = 'Lab')
    
    # Build the plot
    p1 <- ggplot() + dummy_data + data + facet + x_scale + y_scale + lab_x + 
        lab_y + BWscheme + LOHgradient + plot_theme
 
    
    return(p1)
}
