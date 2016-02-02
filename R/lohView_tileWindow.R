#' Obtain LOH data
#' 
#' Obtain LOH heatmap on entire chromsomes from samples in a cohort
#' @name lohView_tileWindow
#' @param loh_data data frame with columns "chromosome", "position", "n_vaf",
#' "t_vaf", "sample" giving raw vaf calls for germline variants
#' @param tile object specifying whether the loh calculations should be tiled
#' @param window_size integer with the size of the sliding window (bp) to be 
#' applied
#' @param normal integer specifying the normal VAF frequency used in LOH 
#' calculations
#' @return object of class dataframe containing LOH data
#' @importFrom plyr adply
#' @importFrom plyr ldply

lohView_tileWindow <- function(loh_data, window_size, normal)
{     
    ## Obtain lists for each sample and chromosome
    out <- split(loh_data, list(as.character(loh_data$chromosome),
                                as.character(loh_data$sample)))
    
    ## Obtain the window position values 
    window_data <- lohView_tilePosition(out, window_size)
    
    ## Perform loh difference calculation
    total <- data.frame()
    final_dataset <- list()
    final_dataset <- lohView_tileCalc(window_data, normal=normal)
    
    ## Combine the lists into a single dataframe
    loh_dataset <- plyr::ldply(final_dataset, data.frame)
    loh_dataset <- loh_dataset[!duplicated(loh_dataset),]
    
    return(loh_dataset)
}

