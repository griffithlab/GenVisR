#' Obtain window information
#' 
#' Calculate window positions to perform LOH calculation
#' @name lohView_windowPosition
#' @param "out" object of class dataframe with columns 'chromosome', 
#' 'position', 'n_freq', 't_freq', and 'sample'
#' @param step integer with the length of divisions (bp) in chromosomes
#' @param window_size integer with the size of the sliding window (bp) to be 
#' applied
#' @return list containing window start/stop positions for each chromosome 
#' from each sample
#' to perform LOH calculations
#' @import dplyr
#' @export


lohView_windowPosition <- function(out, step, window_size) {
    window <- lapply(out, function(x) {
        min <- integer()
        max <- integer()
        window_stop_1 <- integer()
        window_num <- integer()
        min <- as.integer(min(as.numeric(as.character(x$position))))
        max <- as.integer(max(as.numeric(as.character(x$position))))
        window_stop_1 <- min+window_size
        window_num <- as.integer((max-min)/step)
        window_data <- data.frame()
        
        num <- trunc(window_num)
        for (w in 1:num) {
            window_data[w,1] <- as.integer(min + (step*(w-1)))
            window_data[w,2] <- as.integer(window_stop_1 + (step*(w-1)))
        }
        colnames(window_data) <- c("window_start", "window_stop")
        window_final <- dplyr::filter(window_data, 
                                      max > window_data$window_stop)
        window_final[nrow(window_final), 2] <- max
        return(window_final)
    })
    return(window)
}
    
    