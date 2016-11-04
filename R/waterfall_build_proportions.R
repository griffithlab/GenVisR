#'  @title Build mutational profile plot
#' 
#' @param data_frame input data.frame
#' @param plot_palette Color palette to use
#' @param file_type MAF, etc
#' @param layers layer(s) to add to this plot object
#' @param x_label label this plot?
#' 
#' @return a ggplot object
#' 
waterfall_build_proportions <- function(data_frame, plot_palette, file_type, layers,
    x_label) {

    # Declare the appropriate palette
    if(is.null(plot_palette)) {
        palette <- select_palette(file_type)
    } else {
        palette <- plot_palette
    }
    if (x_label) x_label_obj <- xlab(paste0('Sample (n=', nlevels(data_frame$sample), ')'))

    p5 <- ggplot(data_frame, aes_string(x = 'sample', fill = 'trv_type')) + 
        geom_bar(position = "fill", width = 0.95) + 
        scale_y_continuous(name = "% of total mutations", labels = scales::percent_format()) + 
        scale_fill_manual(name="Mutation Type", 
            values=palette, 
            breaks = rev(if (is.null(names(palette))) levels(data_frame[["trv_type"]]) 
                else names(palette))
            ) + 
        guides(fill = guide_legend(title = "Mutation type", ncol = 2)) +
        theme(
            axis.ticks.x = element_blank(),
            axis.text.x = if (x_label) element_text() else element_blank(),
            axis.title.x = element_blank(),
            axis.text.y = element_text(colour = "black"),
            axis.title.y = element_text(colour = "black"),
            panel.background = element_blank(),
            panel.border = element_blank(),
            panel.grid.minor = element_blank(),
            plot.background = element_blank()
        ) + layers + x_label_obj
    p5
}
