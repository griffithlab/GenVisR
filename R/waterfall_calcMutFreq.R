#' Calculate Synonymous/Nonsynonymous mutation frequency
#' 
#' Creates a data frame giving synonymous/nonsynonymous counts on a sample level
#' @name waterfall_calcMutFreq
#' @param x data frame in long format with columns sample, trv_type
#' @return a data frame with synonymous/nonsynonymous counts appended 

waterfall_calcMutFreq <- function(x)
{
    message("Calculating frequency of mutations...")
    # Change trv_type calls to either synonymous or non synonymous,
    # for use in the mutation per Mb plot
    x$trv_type <- as.character(x$trv_type)
    x$trv_type[toupper(x$trv_type) != toupper('silent')] <- 'Non Synonymous'
    x$trv_type[toupper(x$trv_type) == toupper('silent')] <- 'Synonymous'
    x$trv_type <- factor(x$trv_type, levels=c('Synonymous', 'Non Synonymous'))
    
    # Obtain a data frame of mutation counts on the sample level
    mutation_counts <- table(x[,c('sample', 'trv_type')])
    mutation_counts <- as.data.frame(reshape2::melt(mutation_counts))
    colnames(mutation_counts) <- c('sample', 'trv_type', 'mutation_total')
    
    return(mutation_counts)
}