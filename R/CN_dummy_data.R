#' retrieve and format CN_cohort plot supplemental data
#' 
#' given a genome obtain Start and Stop positions for all chromosomes in the genome
#' @name CN_dummy_data
#' @param genome character string specifying UCSC genome from which data is derived
#' @return object of class data frame

CN_dummy_data <- function(genome)
{
  
  # Obtain data for UCSC genome and extract relevant columns
  data <- get_cytobands(genome)
  data <- data[,c(1,2,3)]
  
  # Obtain max for each chromosome
  maxChrom <- aggregate(chromStart ~ chrom, data=data, max)
  maxChrom <- cbind(maxChrom, maxChrom[,2])
  colnames(maxChrom) <- c('Chromosome', 'Start', 'Stop')
  
  # Obtain max for each chromosome
  minChrom <- aggregate(chromStart ~ chrom, data=data, min)
  minChrom <- cbind(minChrom, minChrom[,2])
  colnames(minChrom) <- c('Chromosome', 'Start', 'Stop')
  
  # bind all the data together
  data <- rbind(maxChrom, minChrom)
  
  return(data)
}