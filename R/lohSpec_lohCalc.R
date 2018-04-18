#' Calculate loh difference
#' 
#' Obtain LOH on an entire chromsomes from samples in a cohort
#' @name lohSpec_lohCalc
#' @param window_data object of class data frame with columns 
#' 'window_start' and 'window_stop
#' @param out object of class dataframe with columns 'chromosome', 
#' 'position', 'n_vaf', 't_vaf', and 'sample'
#' @param normal integer specifying the subtraction value from tumor VAF
#' @return object of class dataframe containing mean LOH difference calculations
#' and column names "window_start", "window_stop", "chromosome", "position", 
#' "n_vaf", "t_vaf", "sample", "loh_diff"
#' @noRd

lohSpec_lohCalc <- function(window_data, out, normal)
    
{
    ## Set each list as a dataframe
    window <- data.frame(window_data)   
    sample_data <- data.frame(out)    
    sample_data$position <- as.numeric(as.character(sample_data$position))
    total <- data.frame()
    window_start <- vector()
    window_stop <- vector()
    num <- as.vector(1:nrow(window))
    
    df <- do.call("rbind", lapply(num, function(i){
        ## filter sample_data for loh calls whose positions fall in the
        ## window frame
        filtered_data <- 
            sample_data[sample_data$position >= window$window_start[i] &                           
                            sample_data$position <= window$window_stop[i],] 
        
        ## Perform the loh calculation to obtain avg loh in the window's frame
        loh_calc <- abs(as.numeric(as.character(filtered_data$t_vaf)) - normal)       
        loh_avg <- mean(loh_calc)
        window_start <- window$window_start[i]
        window_stop <- window$window_stop[i]
        if (is.na(loh_avg)==TRUE)    
        {            
            loh_avg <- NULL
            window_start <- NULL
            window_stop <- NULL
        }
        
        filtered_data$loh_diff_avg <- loh_avg
        filtered_data$window_start <- window_start
        filtered_data$window_stop <- window_stop
        return(filtered_data)
    }))

    return(df)
} 



