#' build label for plot
#'
#' given a name create a label
#' @name genCov_buildTrack
#' @param name character string giving the name of the track
#' @param bg_fill character string giving the colour to fill the label
#' @param text_fill character string giving the colour to fill the text
#' @param border character string specifying the colour to fill the border of
#' the label
#' @param size integer specifying the size of the text within the label
#' @return ggplot object
#' @import ggplot2

genCov_buildTrack <- function(name, bg_fill="black", text_fill="white",
                              border="black", size=10)
{
    # Define various parameters of the plot
    background <- geom_rect(colour=border, fill=bg_fill)
    label <- geom_text(aes(x=.5, y=.5), label=name, colour=text_fill, angle=90,
                       size=size)
    scale_x <- scale_x_continuous(expand=c(0,0))
    scale_y <- scale_y_continuous(expand=c(0,0))
    labels <- labs(x=NULL, y=NULL)

    # Define theme of plot
    theme <- theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
                   axis.ticks.x=element_blank(), axis.title.y=element_blank(),
                   axis.text.y=element_blank(), axis.ticks.y=element_blank(),
                   plot.margin=unit(c(0, 0, 0, 0), "null"),
                   axis.ticks.length=unit(0,"null"),
                   panel.spacing=unit(0,"null"))

    # Define the main plot
    label <- ggplot(mapping=aes(xmin=0, xmax=1, ymin=0, ymax=1)) +
      background + label + theme + scale_x + scale_y + labels

    return(label)
}
