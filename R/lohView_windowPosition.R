#' Obtain window information
#' 
#' Calculate window positions to perform LOH calculation
#' @name lohView_windowPosition
#' @param out object of class dataframe with columns 'chromosome', 
#' 'position', 'n_vaf', 't_vaf', and 'sample'
#' @param step integer with the length of divisions (bp) in chromosomes
#' @param window_size integer with the size of the sliding window (bp) to be 
#' applied
#' @return list containing window start/stop positions for each chromosome 
#' from each sample to perform LOH calculations


lohView_windowPosition <- function(out, step, window_size)
{
    window <- lapply(out, function(x) {
        ## Calculate the number of windows necessary for each list
        min <- integer()
        max <- integer()
        window_stop_1 <- integer()
        window_num <- integer()
        min <- as.integer(min(as.numeric(as.character(x$position))))
        max <- as.integer(max(as.numeric(as.character(x$position))))
        window_stop_1 <- min+window_size
        num <- as.integer((max-min)/step)
        num <- as.vector(1:num)
        window_data_start <- vector()
        window_data_stop <- vector()
        
        ## Calculate exact window positions
        window_data <- lapply(num, function(x){
            window_data_start[x] <- as.integer(min+(step*(x-1)))
            window_data_stop[x] <- as.integer(window_stop_1+(step*(x-1)))
            window_data <- cbind(window_data_start[x], window_data_stop[x])
            return(window_data)
        })
        window_data <- plyr::ldply(window_data, data.frame)
        # Get window positions whose values are below max & set max as the 
        # final window position (end of the chromosome)
        colnames(window_data) <- c("window_start", "window_stop")
        window_final <- window_data[window_data$window_stop < max,]
        window_final[nrow(window_final), 2] <- max
        return(window_final)
    })
    return(window)
}
    
    