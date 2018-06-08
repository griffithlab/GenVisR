################################################################################
##################### Public/Private Class Definitions #########################

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Public Class !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

#' Class Rainfall
#' 
#' An S4 class for the Rainfall plot object, under development!!!
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
                                       primaryData="data.table"),
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
#' , should be of length 7. If NULL a default palette will be used.
#' @param sectionHeights Numeric vector specifying relative heights of each plot section,
#' should sum to one. Expects a value for each section.
#' @param chromosomes Character vector specifying chromosomes for which to plot
#' @param sample Character vector specifying the samples for which to plot.
#' @param pointSize numeric value giving the size of points to plot (defaults to 2)
#' @param verbose Boolean specifying if status messages should be reported.
#' @param plotALayers list of ggplot2 layers to be passed to the rainfall plot.
#' @param plotBLayers list of ggplot2 layers to be passed to the density plot.
#' @export
Rainfall <- function(object, BSgenome=NULL, palette=NULL, sectionHeights=NULL, chromosomes=NULL,
                     sample=NULL, pointSize=NULL, verbose=FALSE, plotALayers=NULL, plotBLayers=NULL){
    
    message("This function is part of the new S4 feature and is under active development")
    
    # construct the RainfallPrimaryData object
    primaryData <- RainfallPrimaryData(object, BSgenome=BSgenome, chromosomes=chromosomes,
                                       sample=sample, verbose=verbose)
    
    # construct the Rainfallplots object
    rainfallPlots <- RainfallPlots(primaryData, palette=palette, pointSize=pointSize, plotALayers=plotALayers, plotBLayers=plotBLayers, verbose=verbose)
    
    # align the plots
    rainfallGrob <- arrangeRainfallPlot(rainfallPlots, sectionHeights=sectionHeights, verbose=verbose)
    
    new("Rainfall", PlotA=getGrob(rainfallPlots, index=1), PlotB=getGrob(rainfallPlots, index=2), Grob=rainfallGrob, primaryData=getData(primaryData, name="primaryData"))
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
RainfallPrimaryData <- function(object, BSgenome, chromosomes, sample, verbose){
    
    # convert object to Rainfall format
    primaryData <- toRainfall(object, BSgenome=BSgenome, verbose=verbose)
    
    # annotate each variant with the transition/transversion type
    primaryData <- annoRainfall(primaryData, verbose=verbose)
    
    # calculate the distances between each mutation
    primaryData <- calcMutDist(primaryData, verbose=verbose)
    
    # subset data to only the chromosomes desired to be plotted
    primaryData <- chrSubset(primaryData, chromosomes=chromosomes, verbose=verbose)
    
    # add in chromosome boundaries from a BSgenome object
    primaryData <- annoGenomeCoord(primaryData, BSgenome=BSgenome, verbose=verbose)
    
    # alter the data to highlight the sample
    primaryData <- formatSample(primaryData, sample=sample, verbose=verbose)
    
    new("RainfallPrimaryData", primaryData=primaryData)
}

#' Private Class RainfallPlots
#' 
#' An S4 class for the grobs of the Rainfall plot object
#' @name RainfallPlots-class
#' @rdname RainfallPlots-Class
#' @slot PlotA gtable object for the rainfall plot
#' @slot PlotB gtable object for density plots based on the rainfall plot
#' @import methods
#' @importFrom gtable gtable
#' @noRd
setClass("RainfallPlots",
         representation=representation(PlotA="gtable",
                                       PlotB="gtable"),
         validity=function(object){
             
         }
)

#' Constructor for the RainfallPlots class
#' 
#' @name RainfallPlots
#' @rdname RainfallPlots-class
#' @param object Object of class RainfallPrimaryData
#' @noRd
RainfallPlots <- function(object, palette=palette, pointSize=pointSize, plotALayers=plotALayers, plotBLayers=plotBLayers, verbose=verbose){
    
    # build the rainfall plot
    PlotA <- buildRainfallPlot(object, palette=palette, pointSize=pointSize, plotALayers=plotALayers, verbose=verbose)
    
    # build the corresponding density plot
    PlotB <- buildDensityPlot(object, palette=palette, plotBLayers=plotBLayers, verbose=verbose)
    
    new("RainfallPlots", PlotA=PlotA, PlotB=PlotB)
}

################################################################################
#################### Accessor function definitions #############################

#' Helper function to get data from classes
#' 
#' @rdname getData-methods
#' @aliases getData
.getData_Rainfall <- function(object, name=NULL, index=NULL, ...){
    
    if(is.null(name) & is.null(index)){
        memo <- paste("Both name and index are NULL, one must be specified!")
        stop(memo)
    }
    
    if(is.null(index)){
        index <- 0
    } else {
        if(index > 1){
            memo <- paste("index out of bounds")
            stop(memo)
        }
    }
    
    if(is.null(name)){
        name <- "noMatch"
    } else {
        slotAvailableName <- c("primaryData")
        if(!(name %in% slotAvailableName)){
            memo <- paste("slot name not found, specify one of:", toString(slotAvailableName))
            stop(memo)
        }
    }
    
    if(name == "primaryData" | index == 1){
        data <- object@primaryData
    }
    
    return(data)
}

#' @rdname getData-methods
#' @aliases getData
setMethod(f="getData",
          signature="RainfallPrimaryData",
          definition=.getData_Rainfall)

#' @rdname getData-methods
#' @aliases getData
setMethod(f="getData",
          signature="Rainfall",
          definition=.getData_Rainfall)

#' Helper function to extract grobs from objects
#'
#' @rdname getGrob-methods
#' @aliases getGrob
#' @noRd
.getGrob_Rainfall <- function(object, index=1, ...){
    if(index == 1){
        grob <- object@PlotA
    } else if(index == 2) {
        grob <- object@PlotB
    } else if(index == 3) {
        grob <- object@Grob
    } else {
        stop("Subscript out of bounds") 
    }
    return(grob)
}

#' @rdname getGrob-methods
#' @aliases getGrob
setMethod(f="getGrob",
          signature="RainfallPlots",
          definition=.getGrob_Rainfall)

#' @rdname getGrob-methods
#' @aliases getGrob
setMethod(f="getGrob",
          signature="Rainfall",
          definition=.getGrob_Rainfall)

#' @rdname drawPlot-methods
#' @aliases drawPlot
#' @importFrom grid grid.draw
#' @importFrom grid grid.newpage
#' @exportMethod drawPlot
setMethod(
    f="drawPlot",
    signature="Rainfall",
    definition=function(object, ...){
        mainPlot <- getGrob(object, index=3)
        grid::grid.newpage()
        grid::grid.draw(mainPlot)
    }
)

################################################################################
#################### Method function definitions ###############################

#' @rdname toRainfall-methods
#' @aliases toRainfall
#' @param object Object of class data.table
#' @param verbose Boolean specifying if status messages should be reported
#' @noRd
setMethod(f="toRainfall",
          signature="data.table",
          definition=function(object, verbose, ...){
              
              # print status message
              if(verbose){
                  memo <- paste("converting", class(object)[1], "to expected Rainfall format")
                  message(memo)
              }
              
              # check for appropriate columns
              expectCol <- c("sample", "chromosome", "start", "stop", "refAllele", "variantAllele")
              if(!all(expectCol %in% colnames(object))){
                  missingCol <- expectCol[!expectCol %in% colnames(object)]
                  memo <- paste("Could not find the following required columns in the input:",
                                toString(missingCol))
                  stop(memo)
              }
              rainfallFormat <- object
              
              # remove cases where a mutation does not exist
              rowCountOrig <- nrow(rainfallFormat)
              rainfallFormat <- rainfallFormat[rainfallFormat$refAllele != rainfallFormat$variantAllele,]
              if(rowCountOrig != nrow(rainfallFormat)){
                  memo <- paste("Removed", rowCountOrig - nrow(rainfallFormat),
                                "entries where the reference matches the variant")
                  warning(memo)
              }
              
              # remove mutations at duplicate genomic mutation as this could artifically increase
              # the density of mutations
              rowCountOrig <- nrow(rainfallFormat)
              rainfallFormat <- rainfallFormat[!duplicated(rainfallFormat[,c("sample", "chromosome", "start", "stop")]),]
              if(rowCountOrig != nrow(rainfallFormat)){
                  memo <- paste("Removed", rowCountOrig - nrow(rainfallFormat),
                                "entries with duplicate genomic positions")
                  warning(memo)
              }
              
              # create a flag column for where these entries came from
              rainfallFormat$origin <- 'mutation'
              
              # make sure everything is of the proper type
              rainfallFormat$sample <- factor(rainfallFormat$sample, levels=unique(rainfallFormat$sample))
              rainfallFormat$chromosome <- factor(rainfallFormat$chromosome, levels=unique(rainfallFormat$chromosome))
              rainfallFormat$start <- as.integer(rainfallFormat$start)
              rainfallFormat$stop <- as.integer(rainfallFormat$stop)
              rainfallFormat$refAllele <- factor(rainfallFormat$refAllele, levels=unique(rainfallFormat$refAllele))
              rainfallFormat$variantAllele <- factor(rainfallFormat$variantAllele, levels=unique(rainfallFormat$variantAllele))
              
              # return the standardized format
              return(rainfallFormat)
              
          })

#' @rdname toRainfall-methods
#' @aliases toRainfall
#' @param object Object of class data.table
#' @param verbose Boolean specifying if status messages should be reported
#' @importFrom data.table as.data.table
#' @noRd
setMethod(f="toRainfall",
          signature="data.frame",
          definition=function(object, verbose, ...){
              
              # print status message
              if(verbose){
                  memo <- paste("converting", class(object), "to data.table")
                  message(memo)
              }
              
              # first convert to data.table
              object <- as.data.table(object)
              
              # then convert to rainfall format
              object <- toRainfall(object, verbose=verbose)
              
              return(object)
              
          })

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
              object$base_change <- paste0(object$refAllele, "2", object$variantAllele)
              
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
              
              # split the data.table into lists by chromosome to avoid
              # calculating inter mutation differences among these variables
              object <- split(object, by=c("chromosome"))
              
              # caluclate log10 of distance between mutations, data.tables must be sorted for this
              a <- function(x){
                  x <- x[order(x$start),]
                  difference <- log10(diff(x$start) + 1)
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
                  message(memo)
              }
              
              # perform quality check on the BSgenome object
              if(is.null(BSgenome)){
                  memo <- paste("BSgenome object is not specified, whole chromosomes",
                                "will not be plotted, this is not recommended!")
                  warning(memo)
                  object$chromosome <- factor(object$chromosome, levels=gtools::mixedsort(unique(as.character(object$chromosome))))
                  
                  return(object)
              } else if(is(BSgenome, "BSgenome")) {
                  if(verbose){
                      memo <- paste("BSgenome passed object validity checks")
                      message(memo)
                  }
              } else {
                  memo <- paste("class of the BSgenome object is", class(BSgenome),
                                "should either be of class BSgenome or NULL",
                                "setting this to param to NULL")
                  stop(memo)
              }
              
              # create a data table of genomic coordinates end positions
              genomeCoord_a <- data.table::as.data.table(seqlengths(BSgenome))
              colnames(genomeCoord_a) <- c("start")
              genomeCoord_a$chromosome <- names(seqlengths(BSgenome))
              genomeCoord_a$stop <- genomeCoord_a$start
              genomeCoord_a$origin <- "chrStop"
              
              # create a data table of genomic coordinate start positions
              genomeCoord_b <- data.table::data.table("start"=1, "stop"=1,
                                                      "chromosome"=genomeCoord_a$chromosome,
                                                      "origin"="chrStart")
              
              # bind the master genome coord data table together
              genomeCoord_a <- data.table::rbindlist(list(genomeCoord_a, genomeCoord_b), use.names=TRUE, fill=TRUE)
              genomeCoord_a <- genomeCoord_a[,c("chromosome", "start", "stop", "origin")]
              
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
              object$chromosome <- factor(object$chromosome, levels=gtools::mixedsort(unique(as.character(object$chromosome))))
              
              # return the new object
              return(object)
          })

#' @rdname RainfallPrimaryData-methods
#' @aliases RainfallPrimaryData
#' @param object Object of class data.table
#' @param sample Character string of sample to emphasize
#' @param verbose Boolean for status updates
#' @return data.table object with calculated mutation distances
#' @importFrom gtools mixedsort
#' @noRd
setMethod(f="formatSample",
          signature="data.table",
          definition=function(object, sample, verbose, ...){
              
              # sort samples in a smart sort
              object$sample <- factor(object$sample, levels=gtools::mixedsort(unique(object$sample)))
              
              # if sample is null do nothing further
              if(is.null(sample)){
                  return(object)
              }
              
              # print status message
              if(verbose){
                  memo <- paste("setting data for sample highlighting")
                  message(memo)
              }
              
              # perform some quality checks
              
              # check that input is a character vector
              if(!is.character(sample)){
                  memo <- paste("input to sample is not a character vector",
                                "attempting to coerce.")
                  warning(memo)
                  sample <- as.character(sample)
              }
              
              # check that the sample to highlight is in the data
              if(!sample %in% unique(object$sample)){
                  memo <- paste("Could not find the sample", toString(sample),
                                "in the data, valid entries are:", toString(unique(object$sample)))
                  stop(memo)
              }
              
              # perform the neccessary alterations to the data for sample highlighting
              index <- object$sample %in% sample
              object <- object[index,]
              
              # return the formated object
              return(object)
          })

#' @rdname buildRainfallPlot-methods
#' @aliases Rainfall
#' @param object Object of class RainfallPrimaryData
#' @return gtable object containg a rainfall plot
#' @noRd
#' @import ggplot2
setMethod(f="buildRainfallPlot",
          signature="RainfallPrimaryData",
          definition=function(object, palette, pointSize, plotALayers, verbose=verbose){
              
              # extract the data we need
              primaryData <- getData(object, name="primaryData")
              
              # print status message
              if(verbose){
                  memo <- paste("Building the Rainfall plot")
                  message(memo)
              }
              
              # perform quality checks
              if(!is.null(plotALayers)){
                  if(!is.list(plotALayers)){
                      memo <- paste("plotALayers is not a list")
                      stop(memo)
                  }
                  
                  if(any(!unlist(lapply(plotALayers, function(x) ggplot2::is.ggproto(x) | ggplot2::is.theme(x) | is(x, "labels"))))){
                      memo <- paste("plotALayers is not a list of ggproto or ",
                                    "theme objects... setting plotALayers to NULL")
                      warning(memo)
                      plotALayers <- NULL
                  }
              }
              
              if(is.null(pointSize)){
                  pointSize <- 2
              } else {
                  if(!is.numeric(pointSize)){
                      memo <- paste("pointSize should be of class numeric, attempting to coerce")
                      warning(memo)
                      pointSize <- as.numeric(pointSize)
                  }
                  if(length(pointSize) != 1){
                      memo <- paste("pointSize should be of length one, using only the first value to set size")
                      warning(memo)
                      pointSize <- pointSize[1]
                  }
              }
              
              # remove na values from plot where distance could not be calculated
              coordData <- primaryData[primaryData$origin != "mutation",]
              primaryData <- primaryData[!is.na(primaryData$log10_difference),]
              
              # add the palette for encoding transitions/transversions
              if(is.null(palette)){
                  palette <- c("#77C55D", "#A461B4", "#C1524B", "#93B5BB", "#4F433F", "#BFA753", "#ff9999")
              } else if(length(palette) != 7){
                  memo <- paste("The input to the parameter \"palette\" is not of length",
                                "7, a color should be specified for each transition/transversion in:",
                                toString(unique(primaryData$TransTranv)), "! using the default palette.")
                  warning(memo)
                  palette <- c("#77C55D", "#A461B4", "#C1524B", "#93B5BB", "#4F433F", "#BFA753", "#ff9999")
              }
              
              ############ start building the plot #############################
              
              # define a palette
              plotPalette <- scale_color_manual("Transitions/Transversions", values=palette, na.value="grey70")
              
              # base theme
              baseTheme <- theme_bw()
              
              # define the plot theme
              plotTheme <- theme(axis.text.x=element_text(angle=60, hjust=1),
                                 legend.position="bottom",
                                 strip.background=element_rect(fill="black"),
                                 strip.text=element_text(color="white", face="bold"))
              
              # define plot labels
              plotLabels <- labs(x="Genomic Position",
                                 y="log10 intra-mutation distance")
              
              # define the geom
              plotGeom <- geom_point(mapping=aes_string(color='trans_tranv'), size=pointSize)
              
              # define the faceting
              plotFacet <- facet_grid(. ~ chromosome, scales="free_x")
              
              # define the plot
              rainfallPlot <- ggplot(data=primaryData, mapping=aes_string(x='start', y='log10_difference'))
              
              # combine plot elements
              rainfallPlot <- rainfallPlot + plotGeom + plotFacet + baseTheme +
                  plotTheme + plotLabels + plotPalette + plotALayers + geom_blank(data=coordData, aes_string(x='start'))
              
              # convert to grob
              rainfallGrob <- ggplotGrob(rainfallPlot)
              return(rainfallGrob)
              
          })

#' @rdname buildDensityPlot-methods
#' @aliases Rainfall
#' @param object Object of class RainfallPrimaryData
#' @return gtable object containg a density plot based on a rainfall
#' @noRd
#' @import ggplot2
setMethod(f="buildDensityPlot",
          signature="RainfallPrimaryData",
          definition=function(object, palette, plotBLayers, verbose=verbose){
              
              # extract the data we need
              primaryData <- getData(object, name="primaryData")
              
              # print status message
              if(verbose){
                  memo <- paste("Building the Density plot")
                  message(memo)
              }

              # perform quality checks
              if(!is.null(plotBLayers)){
                  if(!is.list(plotBLayers)){
                      memo <- paste("plotBLayers is not a list")
                      stop(memo)
                  }
                  
                  if(any(!unlist(lapply(plotBLayers, function(x) ggplot2::is.ggproto(x) | ggplot2::is.theme(x) | is(x, "labels"))))){
                      memo <- paste("plotBLayers is not a list of ggproto or ",
                                    "theme objects... setting plotBLayers to NULL")
                      warning(memo)
                      plotBLayers <- NULL
                  }
              }
              
              # split up primary data into data to plot and data defining coord boundaries
              # in the same way as buildRainfallPlot
              coordData <- primaryData[primaryData$origin != "mutation",]
              primaryData <- primaryData[!is.na(primaryData$log10_difference),]
              
              # remove data by chromosome based on a # of entries threshold
              # this should prevent wierd issues with two few points when doing densit plots
              primaryData_by_chr <- split(primaryData, by="chromosome")
              primaryData_by_chr <- primaryData_by_chr[lapply(primaryData_by_chr, nrow) > 4]
              primaryData_filtered <- data.table::rbindlist(primaryData_by_chr)
              
              # add the palette for encoding transitions/transversions
              if(is.null(palette)){
                  palette <- c("#77C55D", "#A461B4", "#C1524B", "#93B5BB", "#4F433F", "#BFA753", "#ff9999")
              } else if(length(palette) != 7){
                  memo <- paste("The input to the parameter \"palette\" is not of length",
                                "7, a color should be specified for each transition/transversion in:",
                                toString(unique(primaryData$TransTranv)), "! using the default palette.")
                  warning(memo)
                  palette <- c("#77C55D", "#A461B4", "#C1524B", "#93B5BB", "#4F433F", "#BFA753", "#ff9999")
              }
              
              ############ start building the plot #############################
              
              # define a palette
              plotPalette <- scale_color_manual("Transitions/Transversions", values=palette, na.value="grey70")
              
              # base theme
              baseTheme <- theme_bw()
              
              # define the plot theme
              plotTheme <- theme(axis.text.x=element_blank(),
                                 axis.ticks.x=element_blank(),
                                 axis.title.x=element_blank(),
                                 panel.grid.major.x=element_blank(),
                                 panel.grid.minor.x=element_blank(),
                                 panel.grid.minor.y=element_blank(),
                                 legend.position="none",
                                 strip.background=element_rect(fill="black"),
                                 strip.text=element_text(color="white", face="bold"))
              
              # define plot labels
              plotLabels <- labs(y="Mutation Density")
              
              # define the geom, when doing this we need to filter the data
              # where there are two few entries to calculate an nice density
              # this should prevent wierd issues with two few points when doing densit plots
              primaryData_by_chr <- split(primaryData, by="chromosome")
              primaryData_by_chr <- primaryData_by_chr[lapply(primaryData_by_chr, nrow) > 4]
              primaryData_filtered <- data.table::rbindlist(primaryData_by_chr)
                  
              # define the geom
              plotGeom <- geom_density(data=primaryData_filtered, mapping=aes_string(x='start'), fill="lightcyan4", color="lightcyan4")
              
              # define the faceting
              plotFacet <- facet_grid(. ~ chromosome, scales="free_x")
              
              # define the plot
              densityPlot <- ggplot(data=primaryData)
              
              # combine plot elements
              densityPlot <- densityPlot + plotGeom + plotFacet + baseTheme +
                  plotTheme + plotLabels + plotPalette + plotBLayers + geom_blank(data=coordData, aes_string(x='start'))
              
              # convert to grob
              densityGrob <- ggplotGrob(densityPlot)
              return(densityGrob)
          })

#' @rdname arrangeRainfallPlot-methods
#' @aliases arrangeRainfallPlot
#' @param sectionHeights Relative heights of each plot section (should sum to one).
#' @param verbose Boolean specifying if status messages should be printed
#' @noRd
#' @importFrom gridExtra arrangeGrob
setMethod(f="arrangeRainfallPlot",
          signature="RainfallPlots",
          definition=function(object, sectionHeights, verbose, ...){
              
              # grab the data we need
              plotA <- getGrob(object, 1)
              plotB <- getGrob(object, 2)
              
              # Obtain the max width for relevant plots
              plotList <- list(plotB, plotA)
              plotList <- plotList[lapply(plotList, length) > 0]
              plotWidths <- lapply(plotList, function(x) x$widths)
              maxWidth <- do.call(grid::unit.pmax, plotWidths)
              
              # Set the widths for all plots
              for(i in 1:length(plotList)){
                  plotList[[i]]$widths <- maxWidth
              }
              
              # set section heights based upon the number of sections
              defaultPlotHeights <- c(.25, .75)
              
              if(is.null(sectionHeights) & length(plotList) == length(defaultPlotHeights)){
                  sectionHeights <- defaultPlotHeights
              } else if(length(sectionHeights) != length(plotList)){
                  memo <- paste("There are", length(sectionHeights), "section heights provided",
                                "but", length(plotList), "vertical sections...",
                                "using default values!")
                  warning(memo)
                  sectionHeights <- defaultPlotHeights
              } else if(!all(is.numeric(sectionHeights))) {
                  memo <- paste("sectionHeights must be numeric... Using",
                                "default values!")
                  warning(memo)
                  sectionHeights <- defaultPlotHeights
              }
              
              # arrange the final plot
              finalPlot <- do.call(gridExtra::arrangeGrob, c(plotList, list(ncol=1, heights=sectionHeights)))
              
              return(finalPlot)
          })