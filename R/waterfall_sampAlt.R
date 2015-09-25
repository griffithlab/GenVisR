#' mutation sample subset sample based
#' 
#' Alter a mutSpec input file keeping/adding entries in a selection of samples
#' @name waterfall_sampAlt
#' @param x a data frame in long format with columns 'sample', 'trv_type'
#' @param samples character vector giving samples to plot
#' @return a subset data frame

waterfall_sampAlt <- function(x, samples)
{
    message("Retrieving requested samples from supplied data...")
    # Perform quality check
    if(typeof(samples) != 'character' & class(samples) != 'character')
    {
        memo <- paste0("argument supplied to main.samples is not a ",
                       "character vector, attempting to coerce")
        warning(memo)
        samples <- as.character(samples)
    }
    
    # add in samples not in the original data frame x
    if(!all(toupper(samples) %in% toupper(x$sample)))
    {
        memo <- paste0("Detected one or more samples supplied to main.samples ",
                       "not found in x ... adding samples to plot")
        message(memo)
        
        samp_not_in_x <- as.data.frame(samples[which(!(samples %in% x$sample))])
        colnames(samp_not_in_x) <- "sample"
        x <- plyr::rbind.fill(samp_not_in_x, x)
    }
    
    # Subset the data based on the arguments supplied in main.samples,
    # then relevel the data frame
    x <- x[(toupper(x$sample) %in% toupper(samples)), ]
    x <- droplevels(x)
    x$sample <- factor(x$sample, levels=samples)
    
    return(x)
}