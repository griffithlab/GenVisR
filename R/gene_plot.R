#' plot gene
#' 
#' given a Granges object plot genomic features within the Granges object
#' @name gene_plot
#' @param txdb A TxDb object for a genome
#' @param gr A Granges object specifying a region of interest
#' @param genome Object of class BSgenome specifying the genome
#' @param reduce Boolean specifying whether to collapse isoforms in the ROI
#' @param transformIntronic Boolean specifying whether to perform a log transform on intronic space
#' @param output_transInt_table Boolean specifying whether to output a master gene features table instead of a plot when transformIntronic is TRUE
#' @import GenomicRanges
#' @return ggplot object
#' @export

gene_plot <- function(txdb, gr, genome, reduce=FALSE, transformIntronic=FALSE, output_transInt_table=FALSE)
{
  require(plyr)
  library(doMC)
  doMC::registerDoMC(cores=8)
  
  # extract a data frame for each type of gene feature given a transcript database and Granges object as a list
  cds <- formatcds(txdb, gr, genome=genome, reduce=reduce)
  threeUTR <- formatthreeUTR(txdb, gr, genome=genome, reduce=reduce)
  fiveUTR <- formatfiveUTR(txdb, gr, genome=genome, reduce=reduce)
  
  # bind together a data frame of gene features as a list and remove erroneous NA values
  gene_features <- mapply(rbind, cds, threeUTR, fiveUTR, SIMPLIFY=FALSE)
  gene_features <- lapply(gene_features, na.omit)
  
  # obtain xlimits for gene plot, this is overwritten of transformIntronic == TRUE
  xlimits <- c(start(gr), end(gr))
  
  # Create a master table based on an intronic log transform then use the master table as a map for mapping coordinates to transformed space
  if(transformIntronic == TRUE)
  {
    # status message
    message("transforming intronic space")
    
    # Create Master table and return it instead of plot if requested
    master <- mergeRegions(gene_features, gr)
    if(output_transInt_table == TRUE)
    {
      return(master)
    }
    
    # Map the original coordinates into transformed space
    gene_features <- lapply(gene_features, function(x, master) adply(x, 1, map_coord_space, master=master, .parallel=TRUE), master=master)
  }
  
  # Adjust the Y axis gene locations based on the presense of isoforms
  increment <- 0
  for(i in 1:length(gene_features))
  {
    gene_features[[i]]$Upper <- gene_features[[i]]$Upper + increment
    gene_features[[i]]$Lower <- gene_features[[i]]$Lower + increment
    gene_features[[i]]$Mid <- gene_features[[i]]$Mid + increment
    if(transformIntronic == TRUE)
    {
      gene_features[[i]]$trans_segStart <- min(gene_features[[i]]$trans_start)
      gene_features[[i]]$trans_segEnd <- max(gene_features[[i]]$trans_end)
    }
    increment <- increment + 2.2
  }	
  
  # Convert the list object of gene_features into a single data frame and set flag to display x axis in plot (this flag is overwritten if transformIntronic == T)
  gene_features <- as.data.frame(do.call("rbind", gene_features))
  display_x_axis <- TRUE
  
  # Replace the original coordinates with the transformed coordinates
  if(transformIntronic==T)
  { 
    # replace original coordinates with transformed coordinates
    gene_features <- gene_features[,c('trans_start', 'trans_end', 'GC', 'width', 'Type', 'Upper', 'Lower', 'Mid', 'trans_segStart', 'trans_segEnd')]
    colnames(gene_features) <- c('start', 'end', 'GC', 'width', 'Type', 'Upper', 'Lower', 'Mid', 'segStart', 'segEnd')
    
    # set flag to not display x axis values if plot is transformed
    display_x_axis <- FALSE
    
    # Obtain x limits for gene plot based on granges object
    start <- cbind(start(gr), start(gr))
    end <- cbind(end(gr), end(gr))
    temp <- as.data.frame(rbind(start, end))
    colnames(temp) <- c('start', 'end')
    temp <- adply(temp, 1, map_coord_space, master=master)
    xlimits <- c(min(temp$trans_start), max(temp$trans_end))
  }
  
  # construct the gene in gplot
  gene_plot <- build_gene(gene_features, display_x_axis=display_x_axis, x_limits=xlimits)
  
  return(gene_plot)
}