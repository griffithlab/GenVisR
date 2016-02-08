#' Obtain average loh within each step 
#' 
#' Calculate avverage LOH within each step
#' @name lohView_stepCalc
#' @param x object of class dataframe with columns 'window_start', 
#' 'window_stop', 'chromosome', 'position', 'n_vaf', 't_vaf', 'sample', and 
#' 'loh_diff_avg'
#' @param step integer with the length of divisions (bp) in chromosomes
#' @return list containing avg loh calculations for each step interval 

lohView_stepCalc <- function(final_dataset, step) {
    ## Remove irrelevant columns
    loh_df1 <- lapply(final_dataset, function(x) 
        x[!(names(x) %in% c("n_vaf", "t_vaf"))])
    
    step_final <- lapply(loh_df1, function(x){
        
        ## Obtain dataset with out repeating windowsloh
        loh_df <- x
        sample <- unique(loh_df$sample)
        chromosome <- unique(loh_df$chromosome)
        ## Obtain boundaries for each step sized window
        step_startMin <- min(as.integer(as.character(loh_df$window_start)))
        step_startMax <- max(as.integer(as.character(loh_df$window_start)))
        step_stopMin <- min(as.integer(as.character(loh_df$window_stop)))
        step_stopMax <- max(as.integer(as.character(loh_df$window_stop)))
        
        ## Calculate number of step-sized windows needed 
        intervals <- ((step_startMax + step) - step_startMin)/step
        
        step_start <- vector()
        step_stop <- vector()
        
        for (i in 1:intervals) {
            ## Obtain positions for the first step-sized window
            if (i == 1) {
                step_start[i] <- step_startMin
                step_stop[i] <- step_start[1] + step 
            }
            if (i > 1) {
                step_start[i] <- step_start[i-1] + step
                step_stop[i] <- step_stop[i-1] + step
            }
        }
        step_df <- data.frame(cbind(step_start, step_stop))
        
        ## Set the end position of the last step-sized window to be the 
        ## end of the chromosome
        step_df[nrow(step_df),2] <- as.numeric(as.character(step_stopMax))
        step_df$chromosome <- chromosome
        step_df$sample <- sample
        
        ## Get the avg loh in overlapping windows for each step
        for (w in 1:nrow(step_df)) {
            start <- as.numeric(as.character(step_df$step_start[w]))
            stop1 <- as.numeric(as.character(step_df$step_stop[w]))
            
            ## Obtain dataset with loh call positions are within the 
            ## step boundaries
            df <- loh_df[loh_df$position >= start & 
                             loh_df$position <= stop1,]
            if (nrow(df) == 0) {
                step_df$loh_avg[w] <- 0
            }
            if (nrow(df) >= 1) {
                step_df$loh_avg[w] <- mean(df$loh_diff_avg)
            }
        }
        
        ## Remove rows without loh calls made in the step-sized window
        final_step_df <- step_df[step_df$loh_avg > 0,]
        return(final_step_df)  
    })
    return(step_final)
}
