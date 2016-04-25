#' Obtain LOH data
#' 
#' Obtain LOH heatmap on entire chromsomes from samples in a cohort
#' @name lohSpec_slidingWindow
#' @param loh_data data frame with columns "chromosome", "position", "n_vaf",
#' "t_vaf", "sample" giving raw vaf calls for germline variants
#' @param step integer with the length of divisions (bp) in chromosomes
#' @param window_size integer with the size of the sliding window (bp) to be 
#' applied
#' @param normal integer specifying the normal VAF frequency used in LOH 
#' calculations
#' @return object of class dataframe containing LOH data
#' @importFrom plyr adply
#' @importFrom plyr ldply

lohSpec_slidingWindow <- function(loh_data, step, window_size, normal)
{     
    ## Obtain lists for each sample and chromosome
    out <- split(loh_data, list(as.character(loh_data$chromosome),
                                      as.character(loh_data$sample)))

    ## Obtain the window position values 
    window_data <- lohSpec_windowPosition(out, step, window_size)

    total <- data.frame()
    final_dataset <- list()
    final_df <- list()
    loh_df <- list()
    ## Perform loh Calculations on each chromosome and sample within each window
    for (i in 1:length(out))
    {
        final_dataset[[i]] <- lohSpec_lohCalc(window_data[[i]], out[[i]], 
                                              normal)
        #final_dataset[[i]] <- plyr::ldply(final_df[[i]], data.frame)
    }
    
    ## Calculate avg loh for overlapping regions 
    df <- lohSpec_stepCalc(final_dataset, step = step) 
    
    ## Combine the lists into a single dataframe
    loh_dataset <- plyr::ldply(df, data.frame)
    colnames(loh_dataset) <- c("window_start", "window_stop", "chromosome", 
                               "sample", "loh_diff_avg")
    loh_dataset$loh_diff_avg <- loh_dataset$loh_diff_avg
    
    return(loh_dataset)
}

