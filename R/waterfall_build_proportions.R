waterfall_build_proportions <- function(data_frame, layers) {
    p5 <- ggplot(data_frame, aes(x = sample, fill = trv_type)) + 
        geom_bar(position = "fill") + 
        scale_y_continuous(name = "% of total mutations", labels = percent_format()) + 
        guides(fill = guide_legend(title = "Mutation type", ncol = 2)) +
        theme(
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            panel.background = element_blank(),
            panel.grid.minor = element_blank(),
            plot.background = element_blank()) 
        # + layers
    p5
}
