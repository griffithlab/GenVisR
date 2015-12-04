#' Construct a gene-features plot
#'
#' Given a GRanges object specifying a region of interest, plot genomic features
#' within that region.
#' @name geneViz
#' @param txdb Object of class TxDb giving transcription meta data for a genome
#' assembly. See Bioconductor annotation packages.
#' @param gr Object of class GRanges specifying the region of interest and 
#' corresponding to a single gene. See Bioconductor package GRanges.
#' @param genome Object of class BSgenome specifying the genome sequence of
#' interest. See Bioconductor annotation packages.
#' @param reduce Boolean specifying whether to collapse gene isoforms within the
#' region of interest into one representative transcript. Experimental use with
#' caution!
#' @param base Numeric vector of log bases to transform the data
#' corresponding to the elements supplied to the variable transform See details.
#' @param transform Character vector specifying what objects to log transform,
#' accepts "Intron", "CDS", and "UTR" See details.
#' @param isoformSel Character vector specifying the names
#' (from the txdb object) of isoforms within the region of interest to display. 
#' @param labelTranscript Boolean specifying whether to plot the transcript
#' names in the gene plot.
#' @param labelTranscriptSize Integer specifying the size of the transcript
#' name text in the gene plot.
#' @param gene_colour Character string specifying the colour of the gene to be
#' plotted.
#' @param plotLayer Valid ggplot2 layer to be added to the gene plot.
#' @details geneViz is an internal function which will output a list of three
#' elements. As a convenience the function is exported however to obtain the
#' plot from geneViz the user must call the first element of the list. geneViz
#' is intended to plot gene features within a single gene with boundaries
#' specified by the GRanges object, plotting more that one gene is advised
#' against.
#' 
#' GET ALEX TO ADD TRANSFORM DETAILS HERE
#' @return object of class list with list elements containing a ggplot object, 
#' the gene features within the plot as a data frame, and mapping information
#' of the gene features within the ggplot object.
#' @importFrom GenomicRanges seqnames
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomicRanges strand
#' @importFrom IRanges IRanges
#' @importFrom plyr adply
#' @importFrom GenomicRanges mcols
#' @importFrom plyr rbind.fill
#' @importFrom GenomicRanges start
#' @importFrom GenomicRanges end
#' @examples
#' # need transcript data for reference
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#'
#' # need a biostrings object for reference
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' genome <- BSgenome.Hsapiens.UCSC.hg19
#'
#' # need Granges object
#' gr <- GRanges(seqnames=c("chr10"),
#' ranges=IRanges(start=c(89622195), end=c(89729532)), strand=strand(c("+")))
#'
#' # Plot the graphic
#' geneViz(txdb, gr, genome)
#' @export

geneViz <- function(txdb, gr, genome, reduce=FALSE, gene_colour=NULL,
                    base=c(10,2,2), transform=c('Intron','CDS','UTR'),
                    isoformSel=NULL, labelTranscript=TRUE,
                    labelTranscriptSize=4, plotLayer=NULL)
{

    # extract a data frame for each type of gene feature given a transcript
    # database and Granges object as a list
    cds <- geneViz_formatCDS(txdb=txdb, gr=gr, genome=genome, reduce=reduce)
    UTR <- geneViz_formatUTR(txdb=txdb, gr=gr, genome=genome, reduce=reduce)

    # bind together a data frame of gene features as a list and remove erroneous
    # NA values
    if(is.null(cds))
    {
        if(is.null(UTR))
        {
            memo <- paste0("No gene features found!, there should be at least", 
                           "one gene present")
            stop(memo)
        }else{
            gene_features <- UTR
        }
    }else if(is.null(UTR)){
        #AHW: I don't think this would ever happen. But, it's a jungle out there
        gene_features <- cds         
    }else{
        keys <- unique(c(names(cds), names(UTR)))
        gene_features <- setNames(mapply(rbind,
                                         cds[keys],
                                         UTR[keys],
                                         SIMPLIFY=FALSE), keys)
    }
    gene_features <- lapply(gene_features, na.omit)
    
    # give the user the ability to select specific isoforms
    if(!is.null(isoformSel))
    {
        # Check that the isoforms requested for subset exist
        
        # rename list with isoform identifications
        isoformNames <- lapply(gene_features, function(x) unique(x$txname))
        names(gene_features) <- isoformNames
        
        if(!any(names(gene_features) %in% isoformSel))
        {
            memo <- paste0("no isoforms supplied to isoformSel were found", 
                           " within the genomic range supplied to gr.", 
                           "Valid isoform names within this range are: ",
                           toString(names(gene_features)))
            stop(memo)
        } else {
            memo <- paste0("Subsetting genomic features based on",
                           " supplied isoforms")
            message(memo)
        }
        
        # extract the isoforms specified
        gene_features <- gene_features[which(names(gene_features)
                                             %in% isoformSel)]
    }
    
    if(reduce)
    {
        gene_features <- lapply(gene_features, geneViz_mergeTypes)
        gene_features[[1]] <- gene_features[[1]][gene_features[[1]]$width > 0,]
        chr <- as.character(GenomicRanges::seqnames(gr))
        x <- gene_features[[1]][,c('start','end')]
        x$chr <- chr
        gc <- function(chr, s, e){
            gr_temp <- 
                GenomicRanges::GRanges(seqnames=c(chr),
                                       ranges=IRanges::IRanges(start=c(as.integer(s)),
                                                               end=c(as.integer(e))),
                                       strand=GenomicRanges::strand(c("+")))
            
            return(geneViz_calcGC(gr=gr_temp, genome=genome))
        }
        
        l <- apply(x,1,function(x) gc(chr=x[3], s=x[1], e=x[2]))
        # AHW: Ugh. A for loop. I know. If you want to make this better,
        # I would welcome it.
        l_new <- l[[1]]
        for (i in 2:length(l))
        {
            l_new <- c(l_new, l[[i]])
        }
        l <- l_new
        gc <- as.data.frame(GenomicRanges::mcols(l))
        gene_features[[1]]$GC <- plyr::rbind.fill(gc)$GC
    }

    # obtain xlimits for gene plot, this is overwritten if transform.
    xlimits <- c(GenomicRanges::start(gr), GenomicRanges::end(gr))

    # Create a master table based on log transform(s) then use the master table
    # as a map for mapping coordinates to transformed space
    if(!is.null(transform) && length(transform) > 0)
    {
        # status message
        message("Calculating transform")

        # Create Master table and return it instead of plot if requested
        master <- geneViz_mergeRegions(gene_features,
                               gr,
                               transform=transform,
                               base=base)

        # Map the original coordinates into transformed space
        gene_features <- lapply(gene_features,
                                function(x, master) plyr::adply(x, 1,
                                                                geneViz_mapCoordSpace,
                                                                master=master),
                                master=master)
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
            gene_features[[i]]$trans_segStart <-
                min(gene_features[[i]]$trans_start)
            gene_features[[i]]$trans_segEnd <-
                max(gene_features[[i]]$trans_end)
        }
        
        increment <- increment + 2.2
    }

    # Convert the list object of gene_features into a single data frame and set
    # flag to display x axis in plot
    # (this flag is overwritten if space is transformed)
    gene_features <- as.data.frame(do.call("rbind", gene_features))
    display_x_axis <- TRUE

    # Replace the original coordinates with the transformed coordinates
    if(length(transform) > 0)
    {
        # replace original coordinates with transformed coordinates
        gene_features <- gene_features[,c('trans_start', 'trans_end', 'GC',
                                          'width', 'Type', 'Upper', 'Lower',
                                          'Mid', 'trans_segStart',
                                          'trans_segEnd', 'txname')]
        
        colnames(gene_features) <- c('start', 'end', 'GC', 'width', 'Type',
                                     'Upper', 'Lower', 'Mid', 'segStart',
                                     'segEnd', 'txname')

        # set flag to not display x axis values if plot is transformed
        display_x_axis <- FALSE

        # Obtain x limits for gene plot based on granges object
        start <- cbind(GenomicRanges::start(gr), GenomicRanges::start(gr))
        end <- cbind(GenomicRanges::end(gr), GenomicRanges::end(gr))
        temp <- as.data.frame(rbind(start, end))
        colnames(temp) <- c('start', 'end')
        temp <- plyr::adply(temp, 1, geneViz_mapCoordSpace, master=master)
        xlimits <- c(min(temp$trans_start), max(temp$trans_end))
    }

    # construct the gene in gplot
    if(reduce == TRUE || labelTranscript == FALSE)
    {
        gene_plot <- geneViz_buildGene(gene_features,
                                       display_x_axis=display_x_axis,
                                       x_limits=xlimits,
                                       gene_colour=gene_colour,
                                       transcript_name=FALSE, layers=plotLayer)
    } else if(reduce == FALSE && labelTranscript == TRUE) {
        
        gene_plot <- geneViz_buildGene(gene_features,
                                       display_x_axis=display_x_axis,
                                       x_limits=xlimits,
                                       gene_colour=gene_colour,
                                       transcript_name=TRUE,
                                       transcript_name_size=labelTranscriptSize,
                                       layers=plotLayer)
    }

    out <- list('plot' = gene_plot,
                'features' = gene_features,
                'master' = master)
    
    return(out)
}
