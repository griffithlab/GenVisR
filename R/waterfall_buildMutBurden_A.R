#' plot mutation burden
#'
#' plot a barchart showing mutations per MB
#' @name waterfall_buildMutBurden_A
#' @param x a data frame in MAF format
#' @param coverage_space an integer specifying the coverage space in base pairs
#' from which a mutation could occur
#' @param layers Additional ggplot2 layers to plot
#' @return a ggplot object

waterfall_buildMutBurden_A <- function(x, coverage_space, layers=NULL)
{
    # Add in mutations per MB calculation
    x$mutation_per_MB <-
    x$mutation_total/coverage_space * 1000000

    # Alter GGplot2 Theme
    theme <- theme(axis.ticks.x = element_blank(),
                   axis.text.x = element_blank(),
                   axis.title.x = element_blank(),
                   legend.title = element_text(size=14),
                   axis.text.y = element_text(colour = "black"),
                   axis.title.y = element_text(colour = "black"),
                   panel.background = element_blank(),
                   # panel.grid.minor.y = element_line(colour = "black"),
                   panel.grid.major.y = element_line(colour = "grey80"),
                   panel.grid.minor.x = element_blank(),
                   panel.grid.major.x = element_blank(),
                   panel.border = element_rect(fill = NA)
            )

    # Add Legend
    legend <- scale_fill_manual(name="Translational Effect",
                                values=c("Synonymous"="red", "Non Synonymous"="blue"),
                                breaks=c("Synonymous", "Non Synonymous"),
                                drop=FALSE)

    # add y label
    y_label <- ylab('Mutations per MB')

    # ggplot2 call
    p1 <- ggplot(x, aes_string(x='sample', y='mutation_per_MB',
                               fill='trv_type')) +
        geom_bar(stat='identity', alpha=.75, width=1) +
        theme + y_label + legend + layers

    return(p1)
}
