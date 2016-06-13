#' Construct copy-number frequency plot
#'
#' Given a data frame construct a plot to display copy number changes across the
#' genome for a group of samples.
#' @name cnFreq
#' @param x Object of class data frame with rows representing the proportion of
#' CN losses/gains across the genome (default), or actual CN values. The former
#' option must contain columns with the following names "chromosome", "start",
#' "end", "gain", and "loss", and the latter option must contain column names
#' "chromosome", "start", "end", "segmean", and "sample". Windows supplied must
#' be consistent across samples!
#' @param CN_low_cutoff Numeric value representing the point at or below which
#' copy number alterations are considered losses. Only used if x represents CN
#' values.
#' @param CN_high_cutoff Numeric value representing the point at or above which
#' copy number alterations are considered gains. Only used if x represents CN
#' values.
#' @param plot_title Character string specifying the title to display on the
#' plot.
#' @param CN_Loss_colour Character string specifying the colour value for copy
#' number losses.
#' @param CN_Gain_colour Character string specifying the colour value for copy
#' number gains.
#' @param x_title_size Integer specifying the size of the x-axis title.
#' @param y_title_size Integer specifying the size of the y-axis title.
#' @param facet_lab_size Integer specifying the size of the faceted labels
#'  plotted.
#' @param plotLayer Valid ggplot2 layer to be added to the plot.
#' @param out Character vector specifying the the object to output, one of
#' "data", "grob", or "plot", defaults to "plot" (see returns).
#' @details cnFreq will detect the column names present in the data frame
#' supplied to x, and will perform one of the following actions. If "gain" and
#' "loss" columns are detected the raw data will be plotted, if "segmean" and 
#' "sample" columns are detected the frequency of copy-number gains and losses
#' present in the cohort will be calculated and plotted. The `plotLayer`
#' parameter can be used to add an additional layer to the ggplot2 graphic
#' (see vignette).
#' @return One of the following, a dataframe containing data to be
#' plotted, a grob object, or a plot.
#' @importFrom gtools mixedsort
#' @examples
#' # Create data
#' xstart <- seq(0,4990000,length.out=500)
#' xloss <- rep(runif(10,0,0.6),rep(50,10))/1.5
#' xloss <- xloss + jitter(xloss,amount=0.002)
#' x <- data.frame(chromosome=rep(paste0("chr",1:5),rep(500,5)), start=xstart,
#' end=xstart+10000, loss=xloss, gain=(1-xloss))
#' 
#' # Plot the data
#' cnFreq(x)
#' @export

cnFreq <- function(x, CN_low_cutoff=1.5, CN_high_cutoff=2.5, plot_title=NULL,
                   CN_Loss_colour='#002EB8', CN_Gain_colour='#A30000',
                   x_title_size=12, y_title_size=12, facet_lab_size=10,
                   plotLayer=NULL, out="plot")
{
    # Perform quality check on input data
    data <- cnFreq_qual(x)
    x <- data[[1]]
    plotType <- data[[2]]

    # If x contains actual CN values, transform into frequencies
    if(plotType=="freq")
    {
        xuniq <- unique(x[,c("chromosome","start","end")])
        gain = loss = obs = rep(NA, nrow(xuniq))
        
        for(i in 1:nrow(xuniq))
        {
            tmpind = Reduce(intersect,
                            list(which(x$chromosome==xuniq[i,1]),
                                 which(x$start==xuniq[i,2]),
                                 which(x$end==xuniq[i,3])))
            
            gain[i] = sum(x[tmpind,"segmean"] >= CN_high_cutoff, na.rm=TRUE)
            loss[i] = sum(x[tmpind,"segmean"] <= CN_low_cutoff, na.rm=TRUE)
            obs[i] = length(tmpind)
        }
        x <- data.frame(chromosome=xuniq$chromosome,
                        start=xuniq$start,
                        end=xuniq$end,
                        gain=gain,
                        loss=loss,
                        obs=obs)
    }

    # Transform losses to be negative values
    x$loss <- -1*x$loss
    x <- as.data.frame(x)

    # Maintain the order of chromosomes
    chr_order <- as.vector(unique(x$chromosome))
    chr_order <- gtools::mixedsort(chr_order)
    x$chromosome <- factor(x$chromosome, levels=chr_order)
    
    # build the plot
    p1 <- cnFreq_buildMain(x, plotType, plot_title=plot_title,
                           CN_low_colour=CN_Loss_colour,
                           CN_high_colour=CN_Gain_colour,
                           x_lab_size=x_title_size,
                           y_lab_size=y_title_size,
                           facet_lab_size=facet_lab_size,
                           plotLayer=plotLayer)
    
    # Decide what to output
    output <- multi_selectOut(data=list("data"=x), plot=p1, out=out)
    return(output)
}
