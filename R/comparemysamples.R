#' Compare the identities of multiple samples
#' 
#' Given the bam file path, count the number of reads at the 24 SNP locations
#' @name comparemysamples
#' @param x data frame with colnames sample_name, bamfile
#' @param genome Object of class BSgenome specifying the genome
#' @return gridExtra object
#' @export 
## x= data frame, column 1: Names of samples, column 2: bam file paths
comparemysamples <- function(x, genome){
  bams <- as.character(x$bamfile)
  samplenames <- as.character(x$sample_name)
  count_tables <- lapply(bams, bamreadcount, genome=genome)
  plot <- compare_identities_plot(count_tables,samplenames)
  return(plot)
}
