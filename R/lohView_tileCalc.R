#' Calculate loh difference
#' 
#' Obtain LOH on an entire chromsomes from samples in a cohort
#' @name lohView_tileCalc
#' @param window_data object of class data frame with columns "chromosome", 
#' "position", "n_vaf", "t_vaf", "sample", "bin", "window_start", "window_stop"
#' @param normal integer specifying the subtraction value from tumor VAF
#' @return object of class dataframe containing mean LOH difference calculations
#' and column names "window_start", "window_stop", "chromosome", "position", 
#' "n_vaf", "t_vaf", "sample", "loh_diff"

lohView_tileCalc <- function(window_data, normal)
{
    
    ## Tile function to calculate LOH    
    ## Calculate absolute loh difference from normal for each call
    total <- lapply(window_data, function(x) {
        ## Get the tumor vaf difference from 50
        t_vaf <- x$t_vaf
        diff <- abs(t_vaf-50)
        x$loh_diff <- diff
        
        ## Obtain the bin values 
        breaks <- as.character(unique(x$bin))
        
        ## Calculate the average loh difference within each bin
        for (i in 1:length(breaks)) {
            data <- subset(x, as.character(x$bin)==breaks[i])
            data$loh_diff_avg <- mean(data$loh_diff)
            
            ## Merge the datasets for each sample
            if(!exists("new_loh_data")) {
                new_loh_data <- data
            }
            if(exists("new_loh_data")) {
                temp <- data
                new_loh_data <- rbind(new_loh_data, temp)
                rm(temp)
            }
            rm(data)
        }
        return(new_loh_data)
    })
    return(total)
}


