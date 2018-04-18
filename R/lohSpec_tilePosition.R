#' Obtain window information
#' 
#' Calculate window positions to perform LOH calculation
#' @name lohSpec_tilePosition
#' @param out object of class dataframe with columns 'chromosome', 
#' 'position', 'n_vaf', 't_vaf', and 'sample'
#' @param window_size integer with the size of the sliding window (bp) to be 
#' applied
#' @return list containing window start/stop positions for each chromosome 
#' from each sample to perform LOH calculations
#' @noRd


lohSpec_tilePosition <- function(out, window_size)
{
    ## Tiling function
    window <- lapply(out, function(x) {
        min <- as.integer(min(as.numeric(as.character(x$position))))
        max <- as.integer(max(as.numeric(as.character(x$position))))
        bin <- seq(min, max, by=window_size)
        num <- length(bin)
        if (bin[num] < max) {
            bin[num] <- max
        }    
        breaks <- cut(as.numeric(as.character(x$position)), breaks=bin, 
                      include.lowest=TRUE, dig.lab=11)
        x$bin <- breaks
            
        bins <- vector()
        bin_min <- vector()
        bin_max <- vector()
        for(i in 1:length(breaks)) {
            bins[i] <- strsplit(as.character(breaks[i]), 
                                split="[[:punct:]]")
            bin_min[i] <- as.numeric(as.character(bins[[i]][2]))
            bin_max[i] <- as.numeric(as.character(bins[[i]][3]))
        }
        x$min <- bin_min
        x$max <- bin_max
        colnames(x) <- c("chromosome", "position", "n_vaf", "t_vaf", "sample",
                         "bin", "window_start", "window_stop")
        return(x)
    })
    
    return(window)
}

