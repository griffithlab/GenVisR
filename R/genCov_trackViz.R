#' Overlay tracks with plots
#'
#' given a named list of plots, display them on tracks
#' @name genCov_trackViz
#' @param ... named list of ggplot2 plots
#' @param bgFill character string giving the colour to fill the label
#' @param textFill character string giving the colour to fill the text
#' @param border character string specifying the colour to fill the border of
#' the label
#' @param size integer specifying the size of the text within the label
#' @param axis_align character string specifying axis to align plotting space
#' on, one of 'both', 'height', 'width', 'none'
#' @param widthRatio vector of length 2 giving the ratio of track labels to
#' plots
#' @param list boolean specifying whether plots are in a named list or specified
#' individually via ...
#' @return ggplotGrob object
#' @importFrom gridExtra arrangeGrob
#' @importFrom gridExtra grid.arrange

genCov_trackViz <- function(..., bgFill="black", textFill="white", border="black",
                     size=10, axis_align='none', widthRatio=c(1, 10), list=TRUE)
{
    # Grab all the tracks/data to be plotted as a named list, check and
    # correct if list is within a list
    if(list==TRUE)
    {
        data <- list(...)
        data <- data[[1]]
    } else {
        data <- list(...)
    }

    # Build Track labels and store as a list containing grob objects,
    # then convert to a single grob
    labels <- lapply(names(data), genCov_buildTrack, bg_fill=bgFill,
                     text_fill=textFill, border=border, size=size)
    label_plot <- do.call(gridExtra::arrangeGrob,
                          c(lapply(labels, ggplot2::ggplotGrob), ncol=1, nrow=length(labels)))

  # Convert the plots corresponding to track labels to a single grob
    data_plot <- genCov_alignPlot(data, axis=axis_align)

  # arrange the label and data plot
    p1 <- gridExtra::grid.arrange(label_plot, data_plot, ncol=2,
                                  widths=widthRatio)

    return(p1)
}
