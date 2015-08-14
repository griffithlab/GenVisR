#' mutation sample subset sample based
#' 
#' Alter a mutSpec input file keeping/adding entries in a selection of samples
#' @name mutSpec.mutation_sample_subset
#' @param x a data frame in long format with columns 'sample', 'trv_type'
#' @param samples character vector giving samples to plot
#' @return a subset data frame

mutSpec.mutation_sample_subset <- function(x, samples)
{
    # Perform quality check
    if(typeof(samples) != 'character' & class(samples) != 'character')
    {
        warning("argument supplied to main.samples is not a character vector,
                attempting to coerce")
        samples <- as.character(samples)
    }
    
    # add in samples not in the original data frame x
    if(!all(toupper(samples) %in% toupper(x$sample)))
    {
        message("Detected one or more samples supplied to main.samples not
                found in x ... adding samples to plot")
        
        samp_not_in_x <- as.data.frame(samples[which(!(samples %in% x$sample))])
        colnames(samp_not_in_x) <- "sample"
        x <- plyr::rbind.fill(samp_not_in_x, x)
    }
    
    # Subset the data based on the arguments supplied in main.samples,
    # then relevel the data frame
    x <- x[(toupper(x$sample) %in% toupper(samples)), ]
    x$sample <- factor(x$sample, levels=samples)
    
    return(x)
}