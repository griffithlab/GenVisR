#' @title Build mutational profile plot
#' 
#' @param data_frame input data.frame
#' @param plot_palette Color palette to use
#' @param file_type MAF, etc
#' @param layers layer(s) to add to this plot object
#' @param x_label label this plot?
#' 
#' @description Builds a ggplot object showing individuals' mutational profile
#' 
#' @return a ggplot object
#' 
waterfall_build_proportions <- function(data_frame, plot_palette, 
    file_type, layers, x_label) {

    # Declare the appropriate palette
    palette <- waterfall_select_palette(file_type, custom_palette = plot_palette)
    breaks_labels <- waterfall_palette_names(palette, file_type, data_frame)
    breaks <- breaks_labels[["breaks"]]
    labels <- breaks_labels[["labels"]]

    if (x_label) {
        x_label_obj <- xlab(paste0('Sample (n=', nlevels(data_frame$sample), ')'))
    } else {
        x_label_obj <- geom_blank()
    }

    p5 <- ggplot(data_frame, aes_string(x = 'sample', fill = 'trv_type')) + 
        geom_bar(position = "fill", width = 0.95) + 
        scale_y_continuous(
            name = "Proportion", 
            labels = scales::percent_format()) + 
        scale_fill_manual(name="Mutation Type", 
            values=palette, 
            breaks = breaks,
            labels = labels
            ) + 
        guides(fill = guide_legend(title = "Mutation Type", ncol = 2)) +
        theme(
            axis.ticks.x = element_blank(),
            axis.text.x =  element_blank(),
            axis.title.x =if (x_label) element_text() else element_blank(),
            axis.text.y = element_text(colour = "black"),
            axis.title.y = element_text(colour = "black"),
            panel.background = element_blank(),
            panel.border = element_blank(),
            panel.grid.minor = element_blank(),
            plot.background = element_blank()
        ) + layers + x_label_obj
    p5
}
