#' Calculate loh difference
#' 
#' Obtain LOH on an entire chromsomes from samples in a cohort
#' @name lohView_lohCalc
#' @param "window_data" object of class data frame with columns 
#' 'window_start' and 'window_stop
#' @param "out" object of class dataframe with columns 'chromosome', 
#' 'position', 'n_freq', 't_freq', and 'sample'
#' @param "normal" integer specifying the subtraction value from tumor VAF
#' @return object of class dataframe containing mean LOH difference calculations
#' and column names "window_start", "window_stop", "chromosome", "position", 
#' "n_freq", "t_freq", "sample", "loh_diff"
#' @import dplyr
#' @export

lohView_lohCalc <- function(window_data, out, normal) {
    window <- as.data.frame(window_data)
    sample_data <- as.data.frame(out)
    total <- data.frame()
    for (i in 1:nrow(window)) {
        filtered_data <- 
            dplyr::filter(sample_data, 
                          as.numeric(as.character(sample_data$position)) 
                                       >= window$window_start[i] & 
                              as.numeric(as.character(sample_data$position)) 
                                       <= window$window_stop[i])
        loh_calc <- abs(as.numeric(as.character(filtered_data$t_freq)) - normal)
        loh_avg <- mean(loh_calc)
        if (is.na(loh_avg)==TRUE) {
            loh_avg <- NULL
        }
        filtered_data$loh_diff <- loh_avg
        total <- rbind(total, filtered_data)
    }
    return(total)
}



