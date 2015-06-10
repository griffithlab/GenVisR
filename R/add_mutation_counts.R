#' Calculate Synonymous/Nonsynonymous mutation frequency
#' 
#' Creates a data frame giving synonymous/nonsynonymous counts
#' @name add_mutation_counts
#' @param data_frame a data frame in long format with columns sample, trv_type
#' @return a data frame with synonymous/nonsynonymous counts appended 
#' @import reshape2

add_mutation_counts <- function(data_frame)
{
  #############################################################################################################################
  ###################### Function to obtain Mutation Frequency counts on a sample level as a data frame #######################
  #############################################################################################################################

  # Change trv_type calls to either synonymous or non synonymous, for use in the mutation per Mb plot
  data_frame$trv_type <- as.character(data_frame$trv_type)
  data_frame$trv_type[toupper(data_frame$trv_type) != toupper('silent')] <- 'Non Synonymous'
  data_frame$trv_type[toupper(data_frame$trv_type) == toupper('silent')] <- 'Synonymous'
  data_frame$trv_type <- factor(data_frame$trv_type, levels=c('Synonymous', 'Non Synonymous'))
  
  # Obtain a data frame of mutation counts on the sample level
  mutation_counts <- table(data_frame[,c('sample', 'trv_type')])
  mutation_counts <- as.data.frame(melt(mutation_counts))
  colnames(mutation_counts) <- c('sample', 'trv_type', 'mutation_total')
  
  return(mutation_counts)
}