#' Calculate loh difference
#'
#' Obtain LOH on an entire chromsomes from samples in a cohort
#' @name lohView_lohCalc
#' @param window_data object of class data frame with columns
#' 'window_start' and 'window_stop
#' @param out object of class dataframe with columns 'chromosome',
#' 'position', 'n_vaf', 't_vaf', and 'sample'
#' @param normal integer specifying the subtraction value from tumor VAF
#' @return object of class dataframe containing mean LOH difference calculations
#' and column names "window_start", "window_stop", "chromosome", "position",
#' "n_vaf", "t_vaf", "sample", "loh_diff"

lohView_lohCalc <- function(window_data, out, normal)

{
    window <- data.frame(window_data)
    sample_data <- data.frame(out)
    sample_data$position <- as.numeric(as.character(sample_data$position))
    total <- data.frame()

    for (i in 1:nrow(window))
    {
        filtered_data <-
            sample_data[sample_data$position >= window$window_start[i] &
                            sample_data$position <= window$window_stop[i],]
        loh_calc <- abs(as.numeric(as.character(filtered_data$t_vaf)) - normal)
        loh_avg <- mean(loh_calc)
        if (is.na(loh_avg)==TRUE)
        {
            loh_avg <- NULL
        }

        filtered_data$loh_diff_avg <- loh_avg
        total <- rbind(total, filtered_data)
    }
    return(total)

} 
