#' Compare the identities of multiple samples
#' 
#' Given the bam file path, count the number of reads at the 24 SNP locations
#' @name comparemysamples
#' @param x data frame of bam files and associated sample names
#' @return gridExtra object
#' @export 

comparemysamples <- function(x){
  count_tables <- lapply(as.character(x[,1]), bamreadcount)
  plot <- compare_identities_plot(count_tables,as.character(x[,2]))
  return(plot)
}
