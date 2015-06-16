#' format mutation observations
#' 
#' Create a data frame of mutation observations
#' @name mutationObs
#' @param data object of class data frame with columns trv_type and amino acid change
#' @param track character string specifying one to 'top', 'bottom' to specify proper track
#' @return object of class data frame giving mutation observations

mutationObs <- function(data, track, fill_value, label_column, rep.fact, rep.dist.lmt, attr.fact, adj.max, adj.lmt, iter.max)
{
  ###################################################################
  ####### Function to extract mutations and their coordinates #######
  ###################################################################
  
  # extract the mutation types and set a flag specifying they are present
  if(any(colnames(data) %in% fill_value))
  {
    fill_value_flag <- TRUE
    fill <- as.character(data[,eval(fill_value)])
  } else {
    fill_value_flag <- FALSE
  }
  
  # extract the mutation coordinates
  mutation_coord <- data$amino_acid_change
  if(all(grepl("p\\.", mutation_coord)))
  {
    message("Detected p. notation for amino_acid_change")
    mutation_coord <- as.numeric(gsub("p\\.[a-zA-z]*(\\d+).*?$", "\\1", mutation_coord, perl=TRUE))
  } else if(all(grepl("c\\.", mutation_coord)))
  {
    stop("C. notation is not currently supported please specify amino acid change in P. notation")
    message("Detected c. notation for amino_acid_change, converting to p. notation")
    mutation_coord <- as.numeric(gsub("c\\.[a-zA-z]*(\\d+).*?$", "\\1", mutation_coord, perl=TRUE))
    fiveUTR_len <- fetchfiveUTRlen(data$transcript_name[1], dataset=ensembl.dataset)
    mutation_coord <- ceiling((mutation_coord - fiveUTR_len)/3)
  } else {
    stop("Could not determine notation type for amino_acid_change, check input")
  }
  
  # combine mutation type and mutation coord into a data frame
  if(fill_value_flag)
  {
    mutation_data <- as.data.frame(cbind(mutation_coord, fill))
    colnames(mutation_data) <- c('mutation_coord', eval(fill_value))
  } else {
    mutation_data <- as.data.frame(mutation_coord)
    colnames(mutation_data) <- c('mutation_coord')
  }
  mutation_data$mutation_coord <- as.numeric(as.character(mutation_data$mutation_coord))
  
  
  # add extra column giving height of Y axis for points to be plotted
  if(track == 'top')
  {
    mutation_data$height_max <- 2
  } else if (track == 'bottom')
  {
    mutation_data$height_min <- -2
  } else {
    stop("Fatal error: incorrect track type specified in mutationObs")
  }
  
  # extract the mutation types and set a flag specifying they are present
  if(any(colnames(data) %in% label_column))
  {
    label_column_flag <- TRUE
    mutation_data$labels <- as.character(data[,eval(label_column)])
  } else {
    label_column_flag <- FALSE
  }
  
  # extract optional labels for points to be plotted
  #if(!is.null(label_column))
  #{
  #  mutation_data$labels <- as.character(data[,eval(label_column)])
  #}
  
  # Dodge mutation coordinates on the x axis
  if(track == 'top')
  {
    message("applying force field to observed mutations for top track")
  } else if (track == 'bottom')
  {
    message("applying force field to observed mutations for bottom track")    
  }
  mutation_data <- mutation_data[order(mutation_coord),] 
  mutation_data$coord_x_dodge <- dodge_coord_x(as.vector(mutation_data$mutation_coord), rep.fact=rep.fact, rep.dist.lmt=rep.dist.lmt, attr.fact=attr.fact, adj.max=adj.max, adj.lmt=adj.lmt, iter.max=iter.max)
  
  # Redefine and return grouping information and then dodge y coordinates
  mutation_data$group <- group_mutation_coord(as.vector(mutation_data$mutation_coord))
  if(track == 'top')
  {
    mutation_data$coord_y_dodge <- dodge_coord_y(mutation_data, track='top')
  } else if(track == 'bottom')
  {
    mutation_data$coord_y_dodge <- dodge_coord_y(mutation_data, track='bottom')
  }
   
  return(mutation_data)
}