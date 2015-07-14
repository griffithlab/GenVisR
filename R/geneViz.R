#' plot gene
#' 
#' given a Granges object plot genomic features within the Granges object
#' @name geneViz
#' @param txdb A TxDb object for a genome
#' @param gr A Granges object specifying a region of interest
#' @param genome Object of class BSgenome specifying the genome
#' @param reduce Boolean specifying whether to collapse isoforms in the ROI
#' @param base  base A vector of log bases to transform the data, corresponding to the elements of transform 
#' @param transform A vector of strings designating what objects to log transform
#' @param plot_transcript_name Boolean specifying whether to plot the transcript name
#' @param transcript_name_size Integer specifying the size of the transcript name text
#' @param gene_colour character string specifying the colour of genomic features plotted
#' @return ggplot object
#' @export
#' @import GenomicRanges
#' @import plyr
#' @import GenomicFeatures
#' @importFrom "IRanges" IRanges

geneViz <- function(txdb, gr, genome, reduce=FALSE, gene_colour=NULL, base=c(10,2,2), transform=c('Intron','CDS','UTR'), plot_transcript_name=TRUE, transcript_name_size=6){

  # extract a data frame for each type of gene feature given a transcript database and Granges object as a list
  cds <- formatCDS(txdb, gr, genome=genome, reduce=reduce)
  UTR <- formatUTR(txdb, gr, genome=genome, reduce=reduce)
  
  # bind together a data frame of gene features as a list and remove erroneous NA values
  if(is.null(cds)){
    if(is.null(UTR)){
      stop('No exons! There should be at least one gene present, otherwise use coverage_plot()')
    }else{
      gene_features <- UTR
    }
  }else if(is.null(UTR)){
    gene_features <- cds          #AHW: I don't think this would ever happen. But, it's a jungle out there.
  }else{
    keys <- unique(c(names(cds), names(UTR)))
    gene_features <- setNames(mapply(rbind, cds[keys], UTR[keys], SIMPLIFY=FALSE), keys)
  }
  gene_features <- lapply(gene_features, na.omit)
  if(reduce){
    gene_features <- lapply(gene_features, mergeTypes)
    chr <- as.character(seqnames(gr))
    x <- gene_features[[1]][,c('start','end')]
    x$chr <- chr
    gc <- function(chr, s, e){
      gr.temp <- GRanges(seqnames=c(chr), ranges=IRanges(start=c(as.integer(s)), end=c(as.integer(e))), strand=strand(c("+")))
      return(calcGC(gr=gr.temp, genome=genome))
    }
    l <- apply(x,1,function(x) gc(chr=x[3], s=x[1], e=x[2]))
    l <- lapply(l, Granges2dataframe)
    gene_features[[1]]$GC <- rbind.fill(l)$GC
  }
  
  # obtain xlimits for gene plot, this is overwritten if transform.
  xlimits <- c(start(gr), end(gr))
  
  # Create a master table based on log transform(s) then use the master table as a map for mapping coordinates to transformed space
  if(!is.null(transform) && length(transform) > 0)
  {
    # status message
    message("Calculating transform")
    
    # Create Master table and return it instead of plot if requested
    master <- mergeRegions(gene_features, gr, transform=transform, base=base)
    
    # Map the original coordinates into transformed space
    gene_features <- lapply(gene_features, function(x, master) adply(x, 1, map_coord_space, master=master), master=master)
  } else {
    master <- NULL
  }
  
  # Adjust the Y axis gene locations based on the presense of isoforms
  increment <- 0
  for(i in 1:length(gene_features))
  {
    gene_features[[i]]$Upper <- gene_features[[i]]$Upper + increment
    gene_features[[i]]$Lower <- gene_features[[i]]$Lower + increment
    gene_features[[i]]$Mid <- gene_features[[i]]$Mid + increment
    if(length(transform) > 0)
    {
      gene_features[[i]]$trans_segStart <- min(gene_features[[i]]$trans_start)
      gene_features[[i]]$trans_segEnd <- max(gene_features[[i]]$trans_end)
    }
    increment <- increment + 2.2
  }	
  
  # Convert the list object of gene_features into a single data frame and set flag to display x axis in plot (this flag is overwritten if space is transformed)
  gene_features <- as.data.frame(do.call("rbind", gene_features))
  display_x_axis <- TRUE
  
  # Replace the original coordinates with the transformed coordinates
  if(length(transform) > 0)
  { 
    # replace original coordinates with transformed coordinates
    gene_features <- gene_features[,c('trans_start', 'trans_end', 'GC', 'width', 'Type', 'Upper', 'Lower', 'Mid', 'trans_segStart', 'trans_segEnd', 'txname')]
    colnames(gene_features) <- c('start', 'end', 'GC', 'width', 'Type', 'Upper', 'Lower', 'Mid', 'segStart', 'segEnd', 'txname')
    
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
  if(reduce == TRUE || plot_transcript_name == FALSE)
  {
    gene_plot <- build.gene(gene_features, display_x_axis=display_x_axis, x_limits=xlimits, gene_colour=gene_colour, transcript_name=FALSE)
  } else if(reduce == FALSE && plot_transcript_name == TRUE)
  {
    gene_plot <- build.gene(gene_features, display_x_axis=display_x_axis, x_limits=xlimits, gene_colour=gene_colour, transcript_name=TRUE, transcript_name_size=transcript_name_size)
  }
  
  out <- list('plot' = gene_plot, 'features' = gene_features, 'master' = master)
  return(out)
}
