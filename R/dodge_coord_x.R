#' dodge coordinates
#' 
#' given amino acid position coordinates on, before grouping and dodge on x axis
#' @name dodge_coord_x
#' @param x numeric vector of position coordinates on x axis
#' @return numeric vector of dodged position coordinates on x axis

dodge_coord_x <- function(x, spread_degree)
{
  # Sort the vector, Then group coordinates meeting pre-defined threshold
  x <- sort(x)
  for(i in 1:length(x))
  {
    # define the first element as belonging to group 1
    if(i == 1)
    {
      group <- 1
      begin <- x[i]
      group_vec <- c(group)
      next
    }
    
    # assign subsequent coordinates to the same group if threshold is met, else start new group
    if(x[i] - x[i - 1] <= 50 & x[i] - begin <= 200)
    {
      group_vec <- c(group_vec, group)
    } else if(x[i] - x[i - 1] > 50 | x[i] - begin > 200) {
      begin <- x[i]
      group <- group + 1
      group_vec <- c(group_vec, group)
    }
  }
  
  # bind group information with vector of coordinates
  coord_data <- as.data.frame(cbind(x, group_vec))
  colnames(coord_data) <- c("coord", "group")
  
  # Obtain median of coordinates for each group and merge onto coordinate data
  median_data <- as.data.frame(aggregate(coord_data, by=list(group_vec), FUN=median))
  median_data <- median_data[,2:3]
  colnames(median_data) <- c("median", "group")
  coord_data <- merge(coord_data, median_data, by="group")
  
  # Obtain a count of coordinates in each group and add to coord data frame
  freq_data <- as.data.frame(table(coord_data$group))
  colnames(freq_data) <- c("group", "Freq")
  coord_data <- merge(coord_data, freq_data, by="group")
  
  # Obtain artificially spread coordinate data and add it to original coordinate data
  spread_coord_data <- apply(coord_data, 1, spread_coord_points, spread_degree)
  spread_coord_data <- unique(spread_coord_data)
  spread_coord_data <- rle(unlist(spread_coord_data))$values
  coord_data <- cbind(coord_data, spread_coord_data)
  
  # correct the skewing for single points introduced by the seq function in spread_coord_points
  spread_coord_data <- apply(coord_data, 1, correct_spread_coord_points)
  coord_data$spread_coord_data <- spread_coord_data
  
  # Re-classify duplicate coordinates as the same coordinate (for stacking)
  coord_dup <- duplicated(coord_data$coord) | duplicated(coord_data$coord, fromLast=TRUE)
  coord_data <- cbind(coord_data, coord_dup)
  # problem below, need to resolve
  coord_data$spread_coord_data[coord_data$coord_dup == TRUE] <- coord_data$coord[coord_data$coord_dup == TRUE]
  
  
  
  # Return the now artificially spread coordinate data
  return(as.vector(coord_data$spread_coord_data))	
}