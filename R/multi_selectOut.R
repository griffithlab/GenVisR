#' Choose output
#' 
#' Selector for choosing output for GenVisR functions
#' 
#' @name multi_selectOut
#' @param data Data object to output
#' @param plot Plot object to output
#' @param out Character vector specifying the the object to output, one of
#' "data", "grob", or "plot".
#' @param draw Boolean specifying if the input to plot needs to be drawn
#' @return One of the following, a list of dataframes containing data to be
#' plotted, a grob object, or a plot.
#' @importFrom ggplot2 ggplotGrob
#' @importFrom grid grid.draw
#' @noRd

multi_selectOut <- function(data, plot, out="plot", draw="FALSE")
{
    # Decide what to output
    if(toupper(out) == "DATA")
    {
        return(data)
    } else if(toupper(out) == "PLOT" & isTRUE(draw)) {
        return(grid::grid.draw(plot))
    } else if(toupper(out) == "PLOT" & !isTRUE(draw)) {
        return(plot)
    } else if(toupper(out) == "GROB" & isTRUE(draw)) {
        return(plot)
    } else if(toupper(out) == "GROB" & !isTRUE(draw)) {
        return(ggplot2::ggplotGrob(plot))
    } else {
        warning("Did not recognize input to out...")
        if(isTRUE(draw))
        {
            return(grid::grid.draw(plot))
        } else {
            return(plot)
        }
    }
}
