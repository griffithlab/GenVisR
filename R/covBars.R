#' Construct coverage cohort plot
#'
#' given a matrix construct a plot to display coverage as percentage bars for a group of samples
#' @name covBars
#' @param x object of class data frame containing rows for the coverage and columns the sample names
#' @param col vector of colors for the coverage bars
#' @param plot_title character string for title of plot
#' @param background character string specifying backround color of plot
#' @param x_lab_size integer specifying the size of the X label
#' @param y_lab_size integer specifying the size of the Y label
#' @param facet_lab_size integer specifying the size of the faceted labels
#' @param layers Additional layers to be plotted, can be a theme but must be a ggplot layer
#' @return ggplot object
#' @examples
#' # Create data
#' x <- matrix(sample(100000,500), nrow=50, ncol=10, dimnames=list(0:49,paste0("Sample",1:10)))
#'
#' # Call function
#' covBars(x)
#' @export
#' @import ggplot2
#' @import reshape2
#' @export

covBars <- function(x, col=NULL, plot_title=NULL, background='grey90', x_lab_size=12, y_lab_size=12, facet_lab_size=10, layers=NULL)
{
  # Perform quality check on input data
    dat <- covBars.qual(x, col)
    x <- dat[[1]]
    col <- dat[[2]]

  # resort the rows (increasing rowname as integer)
    x <- x[order(as.numeric(rownames(x))),]

  # normalize each sample (each sample should sum to 1)
xnorm <- apply(x, 2, function(y){y/sum(as.numeric(y))})

  # get the cumulative sum of each sample
    xcs <- apply(xnorm, 2, cumsum)

  # melt the data for ggplot2 call
    xmelt <- melt(xcs)
    colnames(xmelt) <- c('depth', 'sample', 'bp')

  # define the xmin to be used in the plot (xmax is bp)
    xmelt <- cbind(xmelt, xmin=rep(NA,nrow(xmelt)))
    for(i in unique(xmelt$sample)) {
        tmpcs <- xmelt$bp[xmelt$sample==i]
        xmelt$xmin[xmelt$sample==i] <- c(0, tmpcs[0:(length(tmpcs)-1)])
    }
    xmelt <- as.data.frame(xmelt)

  # Maintain the order of samples
    xmelt$sample <- factor(xmelt$sample, levels=colnames(x))

  # Construct the plot
    p1 <- build.covBars(xmelt, col=col, plot_title=plot_title, background=background, x_lab_size=x_lab_size, y_lab_size=y_lab_size, facet_lab_size=facet_lab_size, layers=layers)

    return(p1)
}
