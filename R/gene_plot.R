#' plot gene
#' 
#' given a Granges object plot genomic features within the Granges object
#' @name gene_plot
#' @param txdb A TxDb object for a genome
#' @param gr A Granges object specifying a region of interest
#' @param genome Object of class BSgenome specifying the genome
#' @param reduce Boolean specifying whether to collapse isoforms in the ROI
#' @param transformIntronic Boolean specifying whether to perform a log transform on intronic space
#' @import GenomicRanges
#' @return ggplot object
#' @export

gene_plot <- function(txdb, gr, genome, reduce=FALSE, transformIntronic=FALSE)
{
  # extract a data frame for each type of gene feature given a transcript database and Granges object as a list
  cds <- formatcds(txdb, gr, genome=genome, reduce=reduce)
  threeUTR <- formatthreeUTR(txdb, gr, genome=genome, reduce=reduce)
  fiveUTR <- formatfiveUTR(txdb, gr, genome=genome, reduce=reduce)
  
  # bind together a data frame of gene features as a list and remove erroneous NA values
  # Here be dragons !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  gene_features <- mapply(rbind, cds, threeUTR, fiveUTR, SIMPLIFY=FALSE)
  gene_features <- lapply(gene_features, na.omit)
  
  #!!!!!!!!!! Test to  grab intronic regions and transform intronic space
  if(transformIntronic == TRUE)
  {
#     introns <- formatIntronic(txdb, gr)
#     gene_features <- format_intronic_space(introns, gene_features, gr)
    master <- mergeRegions(gene_features, gr)
  }
  
  return(gene_features)
  
  #!!!!!!!!!!!!! End of the dragons
  
  # Adjust the Y axis gene locations based on the presense of isoforms
  increment <- 0
  for(i in 1:length(gene_features))
  {
    gene_features[[i]]$Upper <- gene_features[[i]]$Upper + increment
    gene_features[[i]]$Lower <- gene_features[[i]]$Lower + increment
    gene_features[[i]]$Mid <- gene_features[[i]]$Mid + increment
    increment <- increment + 2.2
  }	
  
  # Convert the list object of gene_features into a data frame
  gene_features <- as.data.frame(do.call("rbind", gene_features))
  
  # construct the gene in gplot
  gene_plot <- build_gene(gene_features)
  
  return(gene_plot)
}