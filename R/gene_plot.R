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
  
  # Create a master table based on an intronic log transform then use the master table as a map for mapping coordinates to transformed space
  if(transformIntronic == TRUE)
  {
    message("transforming intronic space")
    master <- mergeRegions(gene_features, gr)
    gene_features <- lapply(gene_features, map_coord_space, master=master)
  }
  
  #!!!!!!!!!!!!! End of the dragons
  
  # Adjust the Y axis gene locations based on the presense of isoforms
  increment <- 0
  for(i in 1:length(gene_features))
  {
    gene_features[[i]]$Upper <- gene_features[[i]]$Upper + increment
    gene_features[[i]]$Lower <- gene_features[[i]]$Lower + increment
    gene_features[[i]]$Mid <- gene_features[[i]]$Mid + increment
    if(transformIntronic == TRUE)
    {
      gene_features[[i]]$trans_segStart <- min(gene_features[[i]]$trans_start_vec)
      gene_features[[i]]$trans_segEnd <- max(gene_features[[i]]$trans_end_vec)
    }
    increment <- increment + 2.2
  }	
  
  # Convert the list object of gene_features into a data frame
  gene_features <- as.data.frame(do.call("rbind", gene_features))
  
  # If it is requested to transform intronic, dump the original values in the data frame and rename the transform data
  if(transformIntronic==T)
  {
    gene_features <- gene_features[,c('trans_start_vec', 'trans_end_vec', 'GC', 'width', 'Type', 'Upper', 'Lower', 'Mid', 'trans_segStart', 'trans_segEnd')]
    colnames(gene_features) <- c('start', 'end', 'GC', 'width', 'Type', 'Upper', 'Lower', 'Mid', 'segStart', 'segEnd')
  }
  
  # construct the gene in gplot
  gene_plot <- build_gene(gene_features)
  
  return(gene_plot)
}