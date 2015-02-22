#' format mutation observations
#' 
#' Create a data frame of mutation observations
#' @name mutationObs
#' @param data object of class data frame with columns trv_type and amino acid change
#' @return object of class data frame giving mutation observations

mutationObs <- function(data)
{
  ###################################################################
  ####### Function to extract mutations and their coordinates #######
  ###################################################################
  
  # extract the mutation types
  mutation <- as.character(data$trv_type)
  
  # extract the mutation coordinates
  mutation_coord <- data$amino_acid_change
  mutation_coord <- as.numeric(gsub("[\\D]+", "", mutation_coord, perl=TRUE))
  
  # combine mutation type and mutation coord into a data frame
  mutation_data <- as.data.frame(cbind(mutation_coord, mutation))
  mutation_data$mutation_coord <- as.numeric(as.character(mutation_data$mutation_coord))
  
  # add extra column giving height of Y axis for points to be plotted
  mutation_data$height_max <- 1.5
  
  return(mutation_data)
}