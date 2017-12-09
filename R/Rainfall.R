################################################################################
##################### Public/Private Class Definitions #########################

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Public Class !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

#' Class Rainfall
#' 
#' An S4 class for the Rainfall plot object
#' @name Rainfall-class
#' @rdname Rainfall-class
#' @slot PlotA gtable object for the rainfall plot
#' @slot PlotB gtable object for density plots based on the rainfall plot
#' @slot Grob gtable object for the arranged plot
#' @slot primaryData data.table object storing the primary data used for plotting.
#' @import methods
#' @importFrom gtable gtable
#' @importFrom data.table data.table
#' @exportClass Rainfall
methods::setOldClass("gtable")
setClass("Rainfall",
         representation=representation(PlotA="gtable",
                                       PlotB="gtable",
                                       Grob="gtable",
                                       primaryData="gtable"),
         validity=function(object){
             
         }
)

#' Constructor for the Rainfall class
#' 
#' @name Rainfall
#' @rdname Rainfall-class
#' @param object Object of class MutationAnnotationFormat, GMS, VEP.
#' @param BSgenome Object of class BSgenome to extract genome wide chromosome coordinates
#' @param palette Character vector specifying colors used for encoding transitions and transversions
#' , should be of length 6. If NULL a default palette will be used.
#' @param sectionHeights Numeric vector specifying relative heights of each plot section,
#' should sum to one. Expects a value for each section.
#' @param chromosomes Character vector specifying chromosomes for which to plot
#' @param highlightSample Character vector specifying the sample for which to highlight, if
#' NULL (the default) no differentiation will be made between samples.
#' @param verbose Boolean specifying if status messages should be reported.
#' @param plotALayers list of ggplot2 layers to be passed to the rainfall plot.
#' @param plotBLayers list of ggplot2 layers to be passed to the density plot.
#' @export
Rainfall <- function(object, BSgenome=NULL, palette=NULL, sectionHeights=NULL, chromosomes=NULL,
                     highlightSample=NULL, verbose=FALSE, plotALayers=NULL, plotBLayers=NULL){
    
    # initalize the RainfallPrimaryData object
    primaryData <- RainfallPrimaryData(object, BSgenome=BSgenome, chromosomes=chromosomes,
                                       highlightSample=highlightSample, verbose=verbose)
    
}

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Private Classes !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

#' Private Class RainfallPrimaryData object
#' 
#' An S4 class for the Primary Data of the Rainfall plot Object
#' @name RainfallPrimaryData-class
#' @rdname RainfallPrimaryData-Class
#' @slot primaryData data.table object storing inter-mutation distances and their classifications
#' @import methods
#' @importFrom data.table data.table
#' @noRd
setClass("RainfallPrimaryData",
         representation=representation(primaryData="data.table"),
         validity=function(object){
             
         }
)

#' Constructor for the RainfallPrimaryData class
#' 
#' @name RainfallPrimaryData
#' @rdname RainfallPrimaryData-class
#' @param object Object of class MutationAnnotationFormat, GMS, VEP
#' @noRd
RainfallPrimaryData <- function(object, BSgenome, chromosomes, highlightSample, verbose){
    
    # convert object to Rainfall format
    primaryData <- toRainfall(object, verbose=verbose)
    
    # annotate each variant with the transition/transversion type
    primaryData <- annoRainfall(primaryData, verbose=verbose)
    
    # calculate the distances between each mutation
    primaryData <- calcMutDist(primaryData, verbose=verbose)
    
    # subset data to only the chromosomes desired to be plotted
    primaryData <- chrSubset(primaryData, chromosomes=chromosomes, verbose=verbose)
    
    # add in chromosome boundaries from a BSgenome object
    primaryData <- annoGenomeCoord(primaryData, BSgenome=BSgenome, verbose=verbose)
    
    # alter the data to highlight the sample
    primaryData <- formatSample(primaryData, highlightSample=highlightSample, verbose=verbose)
    
    new("RainfallPrimaryData", primaryData=primaryData)
}

################################################################################
#################### Accessor function definitions #############################

################################################################################
#################### Method function definitions ###############################

#' @rdname RainfallPrimaryData-methods
#' @aliases RainfallPrimaryData
#' @param object Object of class data.table
#' @param verbose Boolean for status updates
#' @return data.table object with transversions and transitions annotated
#' @noRd
setMethod(f="annoRainfall",
          signature="data.table",
          definition=function(object, verbose, ...){
              
              # print status message
              if(verbose){
                  memo <- paste("annotating mutations as transitions and transversions")
                  message(memo)
              }
              
              # add an extra temporary column acting as a key for transition/transversion type
              object$base_change <- paste0(object$reference, "2", object$variant)
              
              # annotate the variant transition/transversion type
              # annotate the grouping of the base change
              getMutChange <- function(x){
                  switch(x, A2C="A->C or T->G (TV)",
                         T2G="A->C or T->G (TV)", A2G="A->G or T->C (TI)",
                         T2C="A->G or T->C (TI)", A2T="A->T or T->A (TV)",
                         T2A="A->T or T->A (TV)", G2A="G->A or C->T (TI)",
                         C2T="G->A or C->T (TI)", G2C="G->C or C->G (TV)",
                         C2G="G->C or C->G (TV)", G2T="G->T or C->A (TV)",
                         C2A="G->T or C->A (TV)", NA)
              }
              object$trans_tranv <- sapply(object$base_change, getMutChange)
              
              # anything that wasn't annotated above should be an insertion or deletion annotate this
              object[is.na(object$trans_tranv),"trans_tranv"] <- "INDEL"
              
              # set the order of transitions and transversions
              trans_tranv_order <- c("A->C or T->G (TV)", "A->G or T->C (TI)",
                                     "A->T or T->A (TV)", "G->A or C->T (TI)",
                                     "G->C or C->G (TV)", "G->T or C->A (TV)",
                                     "INDEL")
              object$trans_tranv <- factor(object$trans_tranv, levels=trans_tranv_order)
              
              # remove the temporary key column created above
              object$base_change <- NULL
              
              # return the final object
              return(object)
          })

#' @rdname RainfallPrimaryData-methods
#' @aliases RainfallPrimaryData
#' @param object Object of class data.table
#' @param verbose Boolean for status updates
#' @return data.table object with calculated mutation distances
#' @noRd
setMethod(f="calcMutDist",
          signature="data.table",
          definition=function(object, verbose, ...){
              
              # print status message
              if(verbose){
                  memo <- paste("Calculating intra-mutational distances")
                  message(memo)
              }
              
              # split the data.table into lists by sample/chromosome to avoid
              # calculating inter mutation differences among these variables
              object <- split(object, by=c("sample", "chromosome"))
              
              # caluclate log10 of distance between mutations
              a <- function(x){
                  difference <- log10(diff(x$start))
                  difference <- c(NA, difference)
                  x$log10_difference <- difference
                  return(x)
              }
              object <- lapply(object, a)
              
              # recombine dataframe
              object <- data.table::rbindlist(object, use.names=TRUE, fill=TRUE)
              
              return(object)
          })

#' @rdname RainfallPrimaryData-methods
#' @aliases RainfallPrimaryData
#' @param object Object of class data.table
#' @param chromosomes character vector of chromosomes to retain
#' @param verbose Boolean for status updates
#' @return data.table object with calculated mutation distances
#' @noRd
setMethod(f="chrSubset",
          signature="data.table",
          definition=function(object, chromosomes, verbose, ...){
              
              # print status message
              if(verbose){
                  memo <- paste("Performing chromosome subsets")
                  message(memo)
              }
              
              # if chromosomes is null we dont want to do anything just return the object back
              if(is.null(chromosomes)){
                  return(object)
              }
              
              # perform quality checks on the chromosome parameter arguments
              
              # check for character vector
              if(!is.character(chromosomes)){
                  memo <- paste("Input to chromosomes should be a character vector, attempting to coerce...")
                  warning(memo)
              }
              
              # check for specified chromosomes not in the original input
              missingChr <- chromosomes[!chromosomes %in% unique(object$chromosome)]
              if(length(missingChr) != 0){
                  memo <- paste("The following chromosomes were designated to be kept but were not found:",
                                toString(missingChr), "\nValid chromosomes are", toString(unique(object$chromosome)))
                  warning(memo)
              }
              
              # perform the subset
              object <- object[object$chromosome %in% chromosomes,]
              object$chromosome <- factor(object$chromosome)
              
              # check that the object has a size after subsets
              if(nrow(object) < 1){
                  memo <- paste("no entries left to plot after chromosome subsets")
                  stop(memo)
              }
              
              return(object)
          })

#' @rdname RainfallPrimaryData-methods
#' @aliases RainfallPrimaryData
#' @param object Object of class data.table
#' @param BSgenome Object of class BSgenome, used for extracting chromosome boundaries
#' @param verbose Boolean for status updates
#' @return data.table object with calculated mutation distances
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom data.table as.data.table
#' @importFrom data.table rbindlist
#' @importFrom gtools mixedsort
#' @noRd
setMethod(f="annoGenomeCoord",
          signature="data.table",
          definition=function(object, BSgenome, verbose, ...){
              
              # print status message
              if(verbose){
                  memo <- paste("adding chromosome boundaries from BSgenome object")
              }
              
              # perform quality check on the BSgenome object
              if(is.null(BSgenome)){
                  memo <- paste("BSgenome object is not specified, whole chromosomes",
                                "will not be plotted, this is not recommended!")
                  warning(memo)
              } else if(is(BSgenome, "BSgenome")) {
                  if(verbose){
                      memo <- paste("BSgenome passed object validity checks")
                  }
              } else {
                  memo <- paste("class of the BSgenome object is", class(BSgenome),
                                "should either be of class BSgenome or NULL",
                                "setting this to param to NULL")
                  warning(memo)
                  BSgenome <- NULL
              }
              
              # create a data table of genomic coordinates end positions
              genomeCoord_a <- data.table::as.data.table(seqlengths(BSgenome))
              colnames(genomeCoord_a) <- c("start")
              genomeCoord_a$chromosome <- names(seqlengths(BSgenome))
              genomeCoord_a$end <- genomeCoord_a$start
              genomeCoord_a$origin <- "chrStop"
              
              # create a data table of genomic coordinate start positions
              genomeCoord_b <- data.table::data.table("start"=1, "end"=1,
                                                      "chromosome"=genomeCoord_a$chromosome,
                                                      "origin"="chrStart")
              
              # bind the master genome coord data table together
              genomeCoord_a <- data.table::rbindlist(list(genomeCoord_a, genomeCoord_b), use.names=TRUE, fill=TRUE)
              genomeCoord_a <- genomeCoord_a[,c("chromosome", "start", "end", "origin")]
              
              # check that chromosomes between BSgenome and original input match
              chrMismatch <- as.character(unique(object[!object$chromosome %in% genomeCoord_a$chromosome,]$chromosome))
              
              if(length(chrMismatch) >= 1){
                  
                  memo <- paste("The following chromosomes do not match the supplied BSgenome object",
                                toString(chrMismatch))
                  warning(memo)
                  
                  # test if the chr mismatch is fixed by appending chr to chromosomes
                  chrMismatch_appendChr <- length(as.character(unique(object[!paste0("chr", object$chromosome) %in% genomeCoord_a$chromosome,]$chromosome)))
                  if(chrMismatch_appendChr < length(chrMismatch)){
                      memo <- paste("appending \"chr\" to chromosomes in attempt to fix mismatch with the BSgenome")
                      warning(memo)
                      object$chromosome <- paste0("chr", object$chromosome)
                  }
              }
              
              # check if any chromosomes in the origin input lack genomic coordinates
              if(any(!unique(object$chromosome) %in% unique(genomeCoord_a$chromosome))){
                  missingGenomeCoord <- unique(object$chromosome)
                  missingGenomeCoord <- missingGenomeCoord[!missingGenomeCoord %in% unique(genomeCoord_a$chromosome)]
                  memo <- paste("The following chromosomes are missing genomic coordinates", toString(missingGenomeCoord),
                                "Full genomic coordinates will not be plotted for these chromosomes")
                  warning(memo)
              }
              
              # filter the genomeCoord object to only include chromosomes in the input data
              genomeCoord_a <- genomeCoord_a[genomeCoord_a$chromosome %in% unique(object$chromosome),]
              
              # add all genomeCoords to the original input for each sample
              object <- split(object, by=c("sample"))
              combine_Input_Coord <- function(x, y){
                  y$sample <- as.character(unique(x$sample))
                  x <- data.table::rbindlist(list(x, y), use.names=TRUE, fill=TRUE)
              }
              object <- lapply(object, combine_Input_Coord, genomeCoord_a)
              object <- data.table::rbindlist(object, use.names=TRUE, fill=TRUE)
              
              # re-factor such that chromosomes are naturally ordered
              object$chromosome <- factor(object$chromosome, levels=unique(gtools::mixedsort(object$chromosome)))
              
              # return the new object
              return(object)
          })

#' @rdname RainfallPrimaryData-methods
#' @aliases RainfallPrimaryData
#' @param object Object of class data.table
#' @param highlightSample Character string of sample to emphasize
#' @param verbose Boolean for status updates
#' @return data.table object with calculated mutation distances
#' @importFrom gtools mixedsort
#' @noRd
setMethod(f="formatSample",
          signature="data.table",
          definition=function(object, highlightSample, verbose, ...){
              
              # sort samples in a smart sort
              object$sample <- factor(object$sample, levels=gtools::mixedsort(unique(object$sample)))
              
              # if highlightSample is null do nothing further
              if(is.null(highlightSample)){
                  return(object)
              }
              
              # print status message
              if(verbose){
                  memo <- paste("setting data for sample highlighting")
                  message(memo)
              }
              
              # perform some quality checks
              
              # check that input is a character vector
              if(!is.character(highlightSample)){
                  memo <- paste("input to highlightSample is not a character vector",
                                "attempting to coerce.")
                  warning(memo)
                  highlightSample <- as.character(highlightSample)
              }
              
              # check that input is of length 1
              if(length(highlightSample) != 1){
                  memo <- paste("input to highlightSample is not of length 1,",
                                "using only the first element")
                  warning(memo)
                  highlightSample <- highlightSample[1]
              }
              
              # check that the sample to highlight is in the data
              if(!highlightSample %in% unique(object$sample)){
                  memo <- paste("Could not find the sample", toString(highlightSample),
                                "in the data, valid entries are:", toString(unique(object$sample)))
                  stop(memo)
              }
              
              # perform the neccessary alterations to the data for sample highlighting
              object[!object$sample %in% highlightSample,"trans_tranv"] <- NA
              
              # return the formated object
              return(object)
          })
