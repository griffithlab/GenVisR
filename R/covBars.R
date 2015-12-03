#' Construct an overall coverage cohort plot
#'
#' Given a matrix construct a plot to display sequencing depth acheived
#' as percentage bars for a cohort of samples.
#' @name covBars
#' @param x Object of class matrix with rows representing coverage achieved
#' at bases and columns corresponding to each sample in the cohort.
#' @param colour Character vector specifying colours to represent sequencing
#' depth.
#' @param plot_title Character string specifying the title to display on the
#' plot.
#' @param x_title_size Integer specifying the size of the x-axis title.
#' @param y_title_size Integer specifying the size of the y-axis title.
#' @param facet_lab_size Integer specifying the size of the faceted labels
#'  plotted.
#' @param plotLayer Valid ggplot2 layer to be added to the plot.
#' @return ggplot object
#' @importFrom reshape2 melt
#' @examples
#' # Create data
#' x <- matrix(sample(100000,500), nrow=50, ncol=10, dimnames=list(0:49,paste0("Sample",1:10)))
#'
#' # Call plot function
#' covBars(x)
#' @export

covBars <- function(x, colour=NULL, plot_title=NULL, x_title_size=12,
                    y_title_size=12, facet_lab_size=10, plotLayer=NULL)
{
  # Perform quality check on input data
    dat <- covBars_qual(x, colour)
    x <- dat[[1]]
    colour <- dat[[2]]

    # resort the rows (increasing rowname as integer)
    x <- x[order(as.numeric(rownames(x))),]

    # normalize each sample (each sample should sum to 1)
    xnorm <- apply(x, 2, function(y){y/sum(as.numeric(y))})

    # get the cumulative sum of each sample
    xcs <- apply(xnorm, 2, cumsum)

    # melt the data for ggplot2 call
    xmelt <- reshape2::melt(xcs)
    colnames(xmelt) <- c('depth', 'sample', 'bp')

    # define the xmin to be used in the plot (xmax is bp)
    xmelt <- cbind(xmelt, xmin=rep(NA,nrow(xmelt)))
    for(i in unique(xmelt$sample)) 
    {
        tmpcs <- xmelt$bp[xmelt$sample==i]
        xmelt$xmin[xmelt$sample==i] <- c(0, tmpcs[0:(length(tmpcs)-1)])
    }
    xmelt <- as.data.frame(xmelt)

    # Maintain the order of samples
    xmelt$sample <- factor(xmelt$sample, levels=colnames(x))

    # Construct the plot
    p1 <- covBars_buildMain(xmelt, col=colour, plot_title=plot_title,
                            x_lab_size=x_title_size, y_lab_size=y_title_size,
                            facet_lab_size=facet_lab_size, layers=layers)

    return(p1)
}
