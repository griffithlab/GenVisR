################################################################################
##################### Public/Private Class Definitions #########################

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Public Class !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

#' Class Waterfall
#' 
#' An S4 class for the waterfall plot object, under development!!!
#' @name Waterfall-class
#' @rdname Waterfall-class
#' @slot PlotA gtable object for the top sub-plot.
#' @slot PlotB gtable object for the left sub-plot.
#' @slot PlotC gtable object for the main plot.
#' @slot PlotD gtable object for the bottom sub-plot.
#' @slot Grob gtable object for the arranged plot.
#' @slot primaryData data.table object storing the primary data, should have
#' column names sample, gene, mutation, label.
#' @slot simpleMutationCounts data.table object storing simplified mutation
#' counts, should have column names sample, mutation, Freq, mutationBurden
#' @slot complexMutationCounts data.table object storing mutation counts per
#' mutation type should have column names sample, mutation, Freq, mutationBurden.
#' @slot geneData data.table object storing gene counts, should have column
#' names gene, mutation, count.
#' @slot ClinicalData data.table object stroring the data used to plot the
#' clinical sub-plot.
#' @slot mutationHierarchy data.table object storing the hierarchy of mutation
#' type in order of most to least important and the mapping of mutation type to
#' color. Should have column names mutation, color, and label.
#' @exportClass Waterfall
#' @import methods
#' @importFrom gtable gtable
#' @importFrom data.table data.table
methods::setOldClass("gtable")
setClass("Waterfall",
         representation=representation(PlotA="gtable",
                                       PlotB="gtable",
                                       PlotC="gtable",
                                       PlotD="gtable",
                                       Grob="gtable",
                                       primaryData="data.table",
                                       simpleMutationCounts="data.table",
                                       complexMutationCounts="data.table",
                                       geneData="data.table",
                                       ClinicalData="data.table",
                                       mutationHierarchy="data.table"),
         validity=function(object){
         }
         )

#' Constructor for the Waterfall class.
#' 
#' @name Waterfall
#' @rdname Waterfall-class
#' @param input Object of class \code{\link{MutationAnnotationFormat}}, \code{\link{VEP}},
#' \code{\link{GMS}}, or alterantively a data frame/data table with column names "sample", "gene", "mutation".
#' @param labelColumn Character vector specifying a column name from which to
#' extract labels for cells.
#' @param samples Character vector specifying samples to plot. If not NULL
#' all samples in "input" not specified with this parameter are removed.
#' @param coverage Integer specifying the size in base pairs of the genome
#' covered by sequence data from which mutations could be called. Required for
#' the mutation burden sub-plot (see details and vignette). Optionally a named 
#' vector of integers corresponding to each sample can be supplied for more accurate
#' calculations.
#' @param mutation Character vector specifying mutations to keep, if defined
#' mutations not supplied are removed from the main plot.
#' @param mutationHierarchy Data.table object with rows specifying the order of
#' mutations from most to least deleterious and column names "mutation" and
#' "color". Used to change the default colors and/or to give priority to a
#' mutation for the same gene/sample (see details and vignette).
#' @param recurrence Numeric value between 0 and 1 specifying a
#' mutation recurrence cutoff. Genes which do not have mutations in the
#' proportion of samples defined are removed.
#' @param genes Character vector of genes to keep
#' @param geneOrder Character vector specifying the order in which to plot
#' genes.
#' @param geneMax Integer specifying the maximum number of genes to be plotted.
#' Genes kept will be choosen based on the reccurence of mutations in samples.
#' Unless geneOrder is specified.
#' @param sampleOrder Character vector specifying the order in which to plot
#' samples.
#' @param plotA String specifying the type of plot for the top sub-plot, one of
#' "burden", "frequency", or NULL for a mutation burden (requires coverage to be
#' specified), frequency of mutations, or no plot respectively.
#' @param plotALayers list of ggplot2 layers to be passed to the plot.
#' @param plotATally String specifying one of "simple" or "complex" for a
#' simplified or complex tally of mutations respectively.
#' @param plotB String specifying the type of plot for the left sub-plot, one of
#' "proportion", "frequency", or NULL for a plot of gene proportions frequencies
#' , or no plot respectively.
#' @param plotBTally String specifying one of "simple" or "complex" for a
#' simplified or complex tally of genes respectively.
#' @param plotBLayers list of ggplot2 layers to be passed to the plot.
#' @param gridOverlay Boolean specifying if a grid should be overlayed on the
#' waterfall plot.
#' @param drop Boolean specifying if mutations not in the main plot should be dropped from the
#' legend. If FALSE the legend will be based on mutations in the data before any subsets occur.
#' @param labelSize Integer specifying the size of label text
#' @param labelAngle Numeric value specifying the angle of label text
#' @param sampleNames Boolean specifying if samples should be labeled on the plot.
#' @param clinical Object of class Clinical, used for adding a clinical data subplot.
#' @param sectionHeights Numeric vector specifying relative heights of each plot section,
#' should sum to one. Expects a value for each section.
#' @param sectionWidths Numeric vector specifying relative heights of each plot section,
#' should sum to one. Expects a value for each section.
#' @param verbose Boolean specifying if status messages should be reported
#' @param plotCLayers list of ggplot2 layers to be passed to the plot.
#' @seealso \code{\link{MutationAnnotationFormat}}, \code{\link{VEP}}, \code{\link{GMS}}, \code{\link{Clinical}}
#' @details TODO
#' @examples
#' set.seed(426)
#' 
#' # create a data frame with required column names
#' mutationDF <- data.frame("sample"=sample(c("sample_1", "sample_2", "sample_3"), 10, replace=TRUE),
#'                          "gene"=sample(c("egfr", "tp53", "rb1", "apc"), 10, replace=TRUE),
#'                          "mutation"=sample(c("missense", "frame_shift", "splice_site"), 10, replace=TRUE))
#' 
#' # set the mutation hierarchy (required for DF)
#' hierarchyDF <- data.frame("mutation"=c("missense", "frame_shift", "slice_site"),
#'                           "color"=c("#3B3B98", "#BDC581", "#6A006A"))
#'                           
#' # Run the Waterfall Plot and draw the output
#' Waterfall.out <- Waterfall(mutationDF, mutationHierarchy=hierarchyDF)
#' drawPlot(Waterfall.out)
#' @export
Waterfall <- function(input, labelColumn=NULL, samples=NULL, coverage=NULL,
                      mutation=NULL, genes=NULL, mutationHierarchy=NULL,
                      recurrence=NULL, geneOrder=NULL, geneMax=NULL,
                      sampleOrder=NULL, plotA=c("frequency", "burden", NULL),
                      plotATally=c("simple", "complex"), plotALayers=NULL,
                      plotB=c("proportion", "frequency", NULL),
                      plotBTally=c("simple", "complex"), plotBLayers=NULL,
                      gridOverlay=FALSE, drop=TRUE, labelSize=5, labelAngle=0,
                      sampleNames=TRUE, clinical=NULL, sectionHeights=NULL,
                      sectionWidths=NULL, verbose=FALSE, plotCLayers=NULL){
    
    message("This function is part of the new S4 feature and is under active development, did you mean to use waterfall() with a lower case w?")
    
    # calculate all data for plots
    data <- WaterfallData(input, labelColumn=labelColumn, mutationHierarchy=mutationHierarchy,
                          samples=samples, coverage=coverage, mutation=mutation, genes=genes,
                          recurrence=recurrence, geneOrder=geneOrder, geneMax=geneMax,
                          sampleOrder=sampleOrder, verbose=verbose)
    
    # get the clinical data
    if(is.null(clinical)){
        ClinicalData <- data.table::data.table()
    } else {
        ClinicalData <- getData(clinical)
    }
    
    # construct all the plots based on the data
    plots <- WaterfallPlots(data, clinical=clinical, plotA=plotA,
                            plotATally=plotATally, plotALayers=plotALayers,
                            plotB=plotB, plotBTally=plotBTally,
                            plotBLayers=plotBLayers, plotCLayers=plotCLayers,
                            gridOverlay=gridOverlay, drop=drop,
                            labelSize=labelSize, labelAngle=labelAngle,
                            sampleNames=sampleNames, verbose=verbose)

    # align all plots together
    Grob <- arrangeWaterfallPlot(plots, sectionHeights=sectionHeights,
                                 sectionWidths=sectionWidths, verbose=verbose)
    
    new("Waterfall", PlotA=getGrob(plots, index=1), PlotB=getGrob(plots, index=2),
        PlotC=getGrob(plots, index=3), PlotD=getGrob(plots, index=4),
        Grob=Grob, primaryData=getData(data, name="primaryData"),
        simpleMutationCounts=getData(data, name="simpleMutationCounts"),
        complexMutationCounts=getData(data, name="complexMutationCounts"),
        geneData=getData(data, name="geneData"),
        ClinicalData=ClinicalData,
        mutationHierarchy=getData(data, name="mutationHierarchy"))
}

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Private Classes !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

#' Private Class WaterfallData
#' 
#' An S4 class for the Data of the Waterfall plot object
#' @name WaterfallData-class
#' @rdname WaterfallData-class
#' @slot primaryData data.table object storing the primary data, should have
#' column names sample, gene, mutation, label.
#' @slot simpleMutationCounts data.table object storing simplified mutation
#' counts, should have column names sample, mutation, Freq, mutationBurden
#' @slot complexMutationCounts data.table object storing mutation counts per
#' mutation type should have column names sample, mutation, Freq, mutationBurden.
#' @slot geneData data.table object storing gene counts, should have column
#' names gene, mutation, count.
#' @import methods
#' @importFrom data.table data.table
#' @noRd
setClass("WaterfallData",
         representation=representation(primaryData="data.table",
                                       simpleMutationCounts="data.table",
                                       complexMutationCounts="data.table",
                                       geneData="data.table",
                                       mutationHierarchy="data.table"),
         validity=function(object){
         }
)

#' Constructor for the WaterfallData class.
#' 
#' @name WaterfallData
#' @rdname WaterfallData-class
#' @param object Object of class MutationAnnotationFormat
#' @noRd
WaterfallData <- function(object, labelColumn, samples, mutationHierarchy,
                          coverage, mutation, genes, recurrence, geneOrder,
                          geneMax, sampleOrder, verbose){
    
    # assign the mapping of mutations and colors
    mutationHierarchy <- setMutationHierarchy(object, mutationHierarchy, verbose)
    
    # convert to initial data to waterfall format
    primaryData <- toWaterfall(object, mutationHierarchy, labelColumn, verbose)
    
    # subset samples if specified
    primaryData <- sampSubset(primaryData, samples, verbose)
    
    # calculate the frequency and mutation burden
    simpleMutationCounts <- calcSimpleMutationBurden(primaryData, coverage, verbose)
    complexMutationCounts <- calcComplexMutationBurden(primaryData, coverage, verbose)
    
    # remove mutations if specified
    primaryData <- rmvMutation(primaryData, mutation, verbose)
    
    # remove entries for the same gene/sample based on a hierarchy leaving one
    # (must happen before recurrenceSubset)
    primaryData <- mutHierarchySubset(primaryData, mutationHierarchy, verbose)
    
    # Get genes which should be kept based on genes param
    keepGenes_a <- geneSubset(primaryData, genes, verbose)
    
    # get genes which should be kept based on recurrence parameter
    keepGenes_b <- recurrenceSubset(primaryData, recurrence, verbose)
    
    # Filter the necessary genes from geneSubset and recurrenceSubset
    keepGenes <- unique(c(keepGenes_a, keepGenes_b))
    primaryData <- geneFilter(primaryData, keepGenes, verbose)
    
    # set the order of genes for plotting
    primaryData <- orderGenes(primaryData, geneOrder, verbose)
    
    # limit to a maximum number of genes
    primaryData <- maxGeneSubset(primaryData, geneMax, verbose)
    
    # set the order of samples for plotting
    primaryData <- orderSamples(primaryData, sampleOrder, verbose)
    
    # summarize gene level data
    geneData <- constructGeneData(primaryData, verbose)
    
    # initalize the object
    new("WaterfallData", primaryData=primaryData, simpleMutationCounts=simpleMutationCounts,
        complexMutationCounts=complexMutationCounts, geneData=geneData, mutationHierarchy=mutationHierarchy)
}

#' Private Class WaterfallPlots
#' 
#' An S4 class for the grobs of the Waterfall plot object
#' @name WaterfallPlots-class
#' @rdname WaterfallPlots-class
#' @slot PlotA gtable object for the top sub-plot.
#' @slot PlotB gtable object for the left sub-plot.
#' @slot PlotC gtable object for the main plot.
#' @slot PlotD gtable object for the bottom sub-plot.
#' @import methods
#' @importFrom gtable gtable
#' @noRd
setClass("WaterfallPlots",
         representation=representation(PlotA="gtable",
                                       PlotB="gtable",
                                       PlotC="gtable",
                                       PlotD="gtable"),
         validity=function(object){
         }
)

#' Constructor for the WaterfallPlots class.
#' 
#' @name WaterfallPlots
#' @rdname WaterfallPlots-class
#' @param object Object of class WaterfallData
#' @noRd
WaterfallPlots <- function(object, clinical, plotA, plotATally, plotALayers, 
                           plotB, plotBTally, plotBLayers, gridOverlay, drop,
                           labelSize, labelAngle, sampleNames,
                           plotCLayers, verbose){
    
    # create the top sub-plot
    PlotA <- buildMutationPlot(object, plotA, plotATally, plotALayers, verbose)
    
    # create left sub-plot
    PlotB <- buildGenePlot(object, plotB, plotBTally, plotBLayers, verbose)
    
    # add the clinical data
    if(is(clinical, "Clinical")){
        # format the clinical data object using the WaterfallData object
        ClinicalData <- formatClinicalData(object, clinical, verbose)
        
        # make sure the new formating conforms with ClinicalData class 
        ClinicalData <- ClinicalData(ClinicalData, inputFormat="long", verbose=verbose)
        
        # add layer to suppress the x-axis text in the proportion plot
        addLayer <- theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
        plotCLayers[[length(plotCLayers ) + 1]] <- addLayer
        xTitle <- FALSE
        
        # add the clinical data plot
        PlotD <- buildClinicalPlot(ClinicalData, clinicalLayers=getLayers(clinical), verbose=verbose)
    } else {
        PlotD <- gtable::gtable()
        xTitle <- TRUE
    }
    
    # create the main plot
    PlotC <- buildWaterfallPlot(object, gridOverlay, drop, labelSize,
                                labelAngle, xTitle, sampleNames, plotCLayers,
                                verbose)
    
    # initalize the object
    new("WaterfallPlots", PlotA=PlotA, PlotB=PlotB, PlotC=PlotC, PlotD=PlotD)
}


################################################################################
###################### Accessor function definitions ###########################

#' Helper function to getData from classes, under development!!!
#'
#' @rdname getData-methods
#' @aliases getData
.getData_waterfall <- function(object, name=NULL, index=NULL, ...){
    
    if(is.null(name) & is.null(index)){
        memo <- paste("Both name and index are NULL, one must be specified!")
        stop(memo)
    }
    
    if(is.null(index)){
        index <- 0
    } else {
        if(index > 5){
            memo <- paste("index out of bounds")
            stop(memo)
        }
    }
    
    if(is.null(name)){
        name <- "noMatch"
    } else {
        slotAvailableName <- c("primaryData", "simpleMutationCounts", "complexMutationCounts", "geneData", "mutationHierarchy")
        if(!(name %in% slotAvailableName)){
            memo <- paste("slot name not found, specify one of:", toString(slotAvailableName))
            stop(memo)
        }
    }
    
    if(name == "primaryData" | index == 1){
        data <- object@primaryData
    } else if(name == "simpleMutationCounts" | index == 2){
        data <- object@simpleMutationCounts
    } else if(name == "complexMutationCounts" | index == 3){
        data <- object@complexMutationCounts
    } else if(name == "geneData" | index == 4){
        data <- object@geneData
    } else if(name == "mutationHierarchy" | index == 5) {
        data <- object@mutationHierarchy
    }
    
    return(data)
}

#' @rdname getData-methods
#' @aliases getData
setMethod(f="getData",
          signature="WaterfallData",
          definition=.getData_waterfall)

#' @rdname getData-methods
#' @aliases getData
setMethod(f="getData",
          signature="Waterfall",
          definition=.getData_waterfall)

#' Helper function to extract grobs from objects, under development!!!
#'
#' @rdname getGrob-methods
#' @aliases getGrob
#' @noRd
.getGrob_Waterfall <- function(object, index=1, ...){
    if(index == 1){
        grob <- object@PlotA
    } else if(index == 2) {
        grob <- object@PlotB
    } else if(index == 3) {
        grob <- object@PlotC
    } else if(index == 4) {
        grob <- object@PlotD
    } else if(index == 5) {
        grob <- object@Grob
    } else {
        stop("Subscript out of bounds")
    }
    return(grob)
}

#' @rdname getGrob-methods
#' @aliases getGrob
setMethod(f="getGrob",
          signature="WaterfallPlots",
          definition=.getGrob_Waterfall)

#' @rdname getGrob-methods
#' @aliases getGrob
setMethod(f="getGrob",
          signature="Waterfall",
          definition=.getGrob_Waterfall)

#' @rdname drawPlot-methods
#' @aliases drawPlot
#' @importFrom grid grid.draw
#' @importFrom grid grid.newpage
#' @exportMethod drawPlot
setMethod(
    f="drawPlot",
    signature="Waterfall",
    definition=function(object, ...){
        mainPlot <- getGrob(object, index=5)
        grid::grid.newpage()
        grid::grid.draw(mainPlot)
    }
)

################################################################################
####################### Method function definitions ############################

#' @rdname setMutationHierarchy-methods
#' @aliases setMutationHierarchy
#' @noRd
#' @importFrom data.table data.table
#' @importFrom data.table setDT
#' @importFrom grDevices colors
setMethod(f="setMutationHierarchy",
          signature="data.table",
          definition=function(object, mutationHierarchy, verbose, ...){
              
              # if a mutation hierarchy is not specified attempt to create one
              if(is.null(mutationHierarchy)){
                  memo <- paste("mutationHierarchy is null, setting a mutation hierarchy randomly",
                                "this is strongly discouraged!!!!!!!")
                  warning(memo)
                  
                  if(!"mutation" %in% colnames(object)){
                      memo <- paste("column \"mutation\" was not found in input!")
                      stop(memo)
                  } else {
                      mutations <- unique(object$mutation)
                  }
                  
                  newCol <- grDevices::colors(distinct=TRUE)[!grepl("^gray", grDevices::colors(distinct=TRUE))]
                  mutationHierarchy <- data.table::data.table("mutation"=mutations, 
                                                              "color"=sample(newCol, length(mutations)))
                  
              }
                  
              # perform some quality checks on mutationHierarchy
              
              # check that mutationHiearchy is a data table
              if(!any(class(mutationHierarchy) %in% "data.table")){
                  memo <- paste("mutationHiearchy is not an object of class",
                                "data.table, attempting to coerce.")
                  warning(memo)
                  mutationHierarchy <- data.table::setDT(mutationHierarchy)
              }
              
              # check for the correct columns
              if(!all(colnames(mutationHierarchy) %in% c("mutation", "color"))){
                  missingCol <- colnames(mutationHierarchy)[!c("mutation", "color") %in% colnames(mutationHierarchy)]
                  memo <- paste("The correct columns were not found in",
                                "mutationHierarchy, missing", toString(missingCol))
                  stop(memo)
              }
              
              # check that all mutations are specified
              if(!all(object$mutation %in% mutationHierarchy$mutation)){
                  missingMutations <- unique(object$mutation[!object$mutation %in% mutationHierarchy$mutation])
                  memo <- paste("The following mutations were found in the",
                                "input however were not specified in the",
                                "mutationHierarchy!", toString(missingMutations),
                                "adding these in as least important and",
                                "assigning random colors!")
                  warning(memo)
                  newCol <- grDevices::colors(distinct=TRUE)[!grepl("^gray", grDevices::colors(distinct=TRUE))]
                  tmp <- data.table::data.table("mutation"=missingMutations,
                                                "color"=sample(newCol, length(missingMutations)))
                  mutationHierarchy <- data.table::rbindlist(list(mutationHierarchy, tmp), use.names=TRUE, fill=TRUE)
              }
              
              # add in a pretty print mutation labels
              mutationHierarchy$label <- gsub("_", " ", mutationHierarchy$mutation)
              mutationHierarchy$label <-  gsub("'", "' ", mutationHierarchy$label)
              
              # check for duplicate mutations
              if(any(duplicated(mutationHierarchy$mutation))){
                  duplicateMut <- mutationHierarchy[duplicated(mutationHierarchy$mutation),"mutation"]
                  memo <- paste("The mutation type",toString(as.character(duplicateMut)),
                                "was duplicated in the supplied mutationHierarchy!")
                  warning(memo)
                  mutationHierarchy <- mutationHierarchy[!duplicated(mutationHierarchy$mutation),]
              }
              
              # ensure columns are of the proper type
              mutationHierarchy$color <- as.character(mutationHierarchy$color)
              mutationHierarchy$mutation <- as.character(mutationHierarchy$mutation)
              
              # print status message
              if(verbose){
                  memo <- paste("Setting the hierarchy of mutations from most",
                                "to least deleterious and mapping to colors:",
                                toString(mutationHierarchy$mutation))
                  message(memo)
              }
              
              return(mutationHierarchy)
          })

#' @rdname setMutationHierarchy-methods
#' @aliases setMutationHierarchy
#' @noRd
#' @importFrom data.table as.data.table
setMethod(f="setMutationHierarchy",
          signature="data.frame",
          definition=function(object, mutationHierarchy, verbose, ...){
              
              # print status message
              if(verbose){
                  memo <- paste("Converting", class(object),
                                "to data.table in setMutationHierarchy")
                  message(memo)
              }
              
              # convert to data.table
              object <- data.table::as.data.table(object)
              
              # if mutationHierarchy is a data.frame convert that too
              if(is.data.frame(mutationHierarchy)) mutationHierarchy <- data.table::as.data.table(mutationHierarchy)
              
              # check that hierarchy is the proper class
              if(!"data.table" %in% class(mutationHierarchy)){
                  memo <- paste("mutationHierarchy must be an object of class data.table or data.frame!")
                  stop(memo)
              }
              
              # convert to waterfall format
              object <- setMutationHierarchy(object, mutationHierarchy=mutationHierarchy, verbose=verbose)
              
              return(object)
          })

#' @rdname Waterfall-methods
#' @aliases Waterfall
#' @param object Object of class data.table
#' @param verbose Boolean for status updates
#' @noRd
#' @importFrom data.table as.data.table
setMethod(f="toWaterfall",
          signature="data.table",
          definition=function(object, hierarchy, labelColumn, verbose, ...){
              
              # print status message
              if(verbose){
                  memo <- paste("Converting", class(object)[1],
                                "to expected waterfall format")
                  message(memo)
              }
              
              # check for correct columns
              correctCol <- c("sample", "gene", "mutation")
              if(!all(correctCol %in% colnames(object))){
                  missingCol <- correctCol[!correctCol %in% colnames(object)]
                  memo <- paste("Could not find correct column names, missing:",
                                toString(missingCol))
                  stop(memo)
              }
              
              # set up the label variables
              sample <- object$sample
              mutation <- object$mutation
              gene <- object$gene
              label <- NA
              labelFlag <- TRUE
              
              # if a label column exists and is proper overwrite the label variable
              # if not change the flag
              if(!is.null(labelColumn)){
                  if(length(labelColumn) != 1) {
                      memo <- paste("Parameter \"labelColumn\" must be of length 1!",
                                    "Found length to be", length(labelColumn))
                      warning(memo)
                      labelFlag <- FALSE
                  }
                  
                  if(!labelColumn %in% colnames(object)){
                      memo <- paste("could not find column:", labelColumn,
                                    " valid names are:", toString(colnames(object)))
                      warning(memo)
                      labelFlag <- FALSE
                  }
                  
                  if(labelFlag){
                      label <- object[,labelColumn, with=FALSE]
                  }
              }
              
              # combine all columns into a consistent format
              waterfallFormat <- data.table::as.data.table(cbind.data.frame(sample, gene, mutation, label))
              colnames(waterfallFormat) <- c("sample", "gene", "mutation", "label")
              
              # convert appropriate columns to factor
              waterfallFormat$sample <- factor(waterfallFormat$sample)
              
              return(waterfallFormat)
          })

#' @rdname Waterfall-methods
#' @aliases Waterfall
#' @param object Object of class data.table
#' @param verbose Boolean for status updates
#' @importFrom data.table as.data.table
#' @noRd
setMethod(f="toWaterfall",
          signature="data.frame",
          definition=function(object, hierarchy, labelColumn, verbose, ...){
              
              # print status message
              if(verbose){
                  memo <- paste("Converting", class(object),
                                "to data.table")
                  message(memo)
              }
              
              # convert to data.table
              object <- data.table::as.data.table(object)
              
              # convert to waterfall format
              object <- toWaterfall(object, hierarchy=hierarchy, labelColumn=labelColumn, verbose=verbose)
              
              return(object)
          })

#' @rdname Waterfall-methods
#' @aliases Waterfall
#' @param object Object of class data.table
#' @param samples Character vector of samples to keep
#' @param verbose Boolean for status updates
#' @return data.table object subset on "samples" if "samples" is not null.
#' @noRd
setMethod(f="sampSubset",
          signature="data.table",
          definition=function(object, samples, verbose, ...){
              
              # set the object to primary data
              primaryData <- object
              
              # Dont do anything if samples is null
              if(is.null(samples)) return(primaryData)
              
              # print status message
              if(verbose){
                  memo <- paste("Restricting input to specified samples")
                  message(memo)
              }
              
              # make sure the samples variable is as expected
              if(class(samples) != 'character')
              {
                  memo <- paste0("argument supplied to samples is not a ",
                                 "character vector, attempting to coerce")
                  warning(memo)
                  samples <- as.character(samples)
              }
              
              # add in samples not in the original data frame x
              if(!all(samples %in% primaryData$sample))
              {
                  newSamples <- samples[!samples %in% primaryData$sample]
                  memo <- paste("The following samples were specified but not",
                                "found:", toString(newSamples), "adding these",
                               "samples to the input data!")
                  warning(memo)
                  primaryData <- rbind(primaryData,
                                       data.table("sample"=newSamples), fill=TRUE)
              }
              
              # remove any samples not requested to be plotted
              primaryData <- primaryData[primaryData$sample %in% samples,]
              primaryData$sample <- factor(primaryData$sample, levels=samples)
              
              return(primaryData)
          })

#' @rdname Waterfall-methods
#' @aliases Waterfall
#' @param object Object of class waterfall
#' @param coverage Integer specifying the size in base pairs of genome from
#' which mutations could have been called (denominator of mutation burden).
#' @param verbose Boolean for status updates.
#' @return data.table object with simplified mutation counts
#' @noRd
#' @importFrom data.table as.data.table
setMethod(f="calcSimpleMutationBurden",
          signature="data.table",
          definition=function(object, coverage, verbose, ...){
              
              # set primaryData as the object
              primaryData <- object
              
              # status message
              if(verbose){
                  memo <- paste("Calculating frequency and mutation burden.")
                  message(memo)
              }
              
              # quality checks
              if(!all(is.numeric(coverage)) && !is.null(coverage)){
                  memo <- paste("coverage is not numeric, attempting to coerce.")
                  warning(memo)
                  coverage <- as.numeric(coverage)
              }
              
              # if coverage is a vector check that each sample has a corresponding coverage value
              if(length(coverage) > 1 && !is.null(coverage)){
                  
                  # check that there are no duplicate samples in coverage
                  if(any(duplicated(names(coverage)))){
                      memo <- paste("Found duplicate names in coverage, removing duplicates!")
                      warning(memo)
                      coverage <- coverage[!duplicated(names(coverage))]
                  }
                  
                  # check that there are no "extra" samples in coverage
                  if(any(!names(coverage) %in% unique(primaryData$sample))){
                      extraSamples <- names(coverage)[!names(coverage) %in% unique(primaryData$sample)]
                      memo <- paste("The following names in coverage are were not found in the data",
                                    toString(extraSamples), "removing these from coverage!")
                      warning(memo)
                      coverage <- coverage[!names(coverage) %in% extraSamples]
                  }
                  
                  # check that each sample has a value for coverage
                  if(any(!unique(primaryData$sample) %in% names(coverage))){
                      missingSamples <- unique(primaryData$sample)
                      missingSamples <- missingSamples[!missingSamples %in% names(coverage)]
                      memo <- paste("coverage has length > 1 however a coverage value could",
                                    "not be found for all samples, samples missing a coverage",
                                    "value are:", toString(missingSamples), "...setting coverage to NULL!")
                      warning(memo)
                      coverage <- NULL
                  }
              }
              
              # dont include any samples with NA values for genes
              samples <- primaryData[is.na(primaryData$gene),]$sample
              primaryData <- primaryData[!is.na(primaryData$gene),]
              
              # obtain a data table of mutation counts on the sample level
              simpleMutationCounts <- data.table::as.data.table(table(primaryData[,c('sample')], exclude=samples))
              colnames(simpleMutationCounts) <- c("sample", "Freq")
              simpleMutationCounts$mutation <- NA
              
              # add the samples with NA values back in if needed
              if(length(samples) != 0){
                  samples <- data.table::data.table("sample"=samples, "Freq"=0)
                  simpleMutationCounts <- data.table::rbindlist(list(simpleMutationCounts, samples),
                                                                use.names=TRUE, fill=TRUE)
              }
              
              # if coverage is not specified return just frequencies else use it in the calculations
              if(is.null(coverage)){
                  
                  if(verbose){
                      memo <- paste("coverage not specified, could not",
                                    "calculate the mutation burden")
                      message(memo)
                  }      

                  simpleMutationCounts$mutationBurden <- NA
                  simpleMutationCounts$coverage <- NA
                  return(simpleMutationCounts)
              } else if(length(coverage) > 1){
                  
                  # need to match coverage to sample
                  coverage <- data.table::as.data.table(coverage, keep.rownames=TRUE)
                  colnames(coverage) <- c("sample", "coverage")
                  simpleMutationCounts <- merge(simpleMutationCounts, coverage, by=c("sample"))
                  
              } else {
                  
                  # assign just the single coverage value
                  simpleMutationCounts$coverage <- coverage
              }
              
              # mutation burden calculation
              simpleMutationCounts$mutationBurden <- simpleMutationCounts$Freq/simpleMutationCounts$coverage * 1000000
              return(simpleMutationCounts)
          })

#' @rdname Waterfall-methods
#' @aliases Waterfall
#' @param object Object of class waterfall
#' @param coverage Integer specifying the size in base pairs of genome from
#' which mutations could have been called (denominator of mutation burden).
#' @param verbose Boolean for status updates.
#' @return data.table object with complex mutation counts
#' @noRd
#' @importFrom data.table setDT
#' @importFrom data.table data.table
#' @importFrom data.table as.data.table
setMethod(f="calcComplexMutationBurden",
          signature="data.table",
          definition=function(object, coverage, verbose, ...){
              
              # set primaryData as the object
              primaryData <- object
              
              # status message
              if(verbose){
                  memo <- paste("Calculating complex mutation burden.")
                  message(memo)
              }
              
              # quality checks
              if(!all(is.numeric(coverage)) && !is.null(coverage)){
                  memo <- paste("coverage is not numeric, attempting to coerce.")
                  warning(memo)
                  coverage <- as.numeric(coverage)
              }
              
              # if coverage is a vector check that each sample has a corresponding coverage value
              if(length(coverage) > 1 && !is.null(coverage)){
                  
                  # check that there are no duplicate samples in coverage
                  if(any(duplicated(names(coverage)))){
                      memo <- paste("Found duplicate names in coverage, removing duplicates!")
                      warning(memo)
                      coverage <- coverage[!duplicated(names(coverage))]
                  }
                  
                  # check that there are no "extra" samples in coverage
                  if(any(!names(coverage) %in% unique(primaryData$sample))){
                      extraSamples <- names(coverage)[!names(coverage) %in% unique(primaryData$sample)]
                      memo <- paste("The following names in coverage are were not found in the data:",
                                    toString(extraSamples), "...removing these from coverage!")
                      warning(memo)
                      coverage <- coverage[!names(coverage) %in% extraSamples]
                  }
                  
                  # check that each sample has a value for coverage
                  if(any(!unique(primaryData$sample) %in% names(coverage))){
                      missingSamples <- unique(primaryData$sample)
                      missingSamples <- missingSamples[!missingSamples %in% names(coverage)]
                      memo <- paste("coverage has length > 1 however a coverage value could",
                                    "not be found for all samples, samples missing a coverage",
                                    "value are:", toString(missingSamples), "...setting coverage to NULL!")
                      warning(memo)
                      coverage <- NULL
                  }
              }
              
              # obtain a data table of mutation counts on the sample level
              complexMutationCounts <- as.data.frame(table(primaryData[,c('sample', 'mutation')]))
              data.table::setDT(complexMutationCounts)
              
              # if coverage is not specified return just frequencies else use it in the calculations
              if(is.null(coverage)){
                  if(verbose){
                      memo <- paste("coverage not specified, could not",
                                    "calculate the mutation burden")
                      message(memo)
                  }      

                  complexMutationCounts$mutationBurden <- NA
                  complexMutationCounts$coverage <- NA
                  return(complexMutationCounts)
                  
              } else if(length(coverage) > 1){
                  
                  # need to match coverage to sample
                  coverage <- data.table::as.data.table(coverage, keep.rownames=TRUE)
                  colnames(coverage) <- c("sample", "coverage")
                  complexMutationCounts <- merge(complexMutationCounts, coverage, by=c("sample"))
                  
              } else {
                  
                  # assign just the single coverage value
                  complexMutationCounts$coverage <- coverage
              }
              
              # mutation burden calculation
              complexMutationCounts$mutationBurden <- complexMutationCounts$Freq/complexMutationCounts$coverage * 1000000
              return(complexMutationCounts)
          })

#' @rdname Waterfall-methods
#' @aliases Waterfall
#' @param object Object of class waterfall
#' @param verbose Boolean for status updates
#' @return data.table object with mutations removed from primaryData slot.
#' @noRd
setMethod(f="rmvMutation",
          signature="data.table",
          definition=function(object, mutation, verbose, ...){
              
              # set primaryData as the object
              primaryData <- object
              mutation <- unique(mutation)
              
              # do nothing if mutation is null
              if(is.null(mutation)){
                  return(primaryData)
              }
              
              # perform quality checks
              if(!is.character(mutation)){
                  memo <- paste("mutation is defined but is not a character vector,",
                                "attempting to coerce.")
                  mutation <- as.character(mutation)
                  warning(memo)
              }
              
              # store how many mutations there are originaly to figure out how
              # many mutations are removed
              totalMut <- length(na.omit(primaryData$mutation))
              
              # check if any mutations were specified but are not in the data
              if(!any(mutation %in% primaryData$mutation)){
                  missingMutation <- mutation[!mutation %in% unique(primaryData$mutation)]
                  memo <- paste("the following mutations were specified to be kept",
                                "but were not found:", toString(missingMutation))
                  warning(memo)
              }
              
              # remove mutations not specified to be kept, samples should remain
              # because they are in the factor levels
              keep <- primaryData$mutation %in% mutation
              primaryData <- primaryData[keep,]

              # figure out number of mutations that existed
              mutCount <- totalMut - length(na.omit(primaryData$mutation))
              
              # print status message
              if(verbose){
                  memo <- paste("Removed", mutCount, "which were not specified",
                                "to be kept.")
                  message(memo)
              }
              
              return(primaryData)
          })

#' @rdname Waterfall-methods
#' @aliases Waterfall
#' @param object Object of class waterfall
#' @param genes character vector giving genes to keep
#' @param verbose Boolean for status updates
#' @return Character vector of genes which should be kept.
#' @noRd
setMethod(f="geneSubset",
          signature="data.table",
          definition=function(object, genes, verbose, ...){
              # access the part of the object we want to manipulate
              primaryData <- object
              
              # Dont do anything if genes is null
              if(is.null(genes)) return(NA)
              
              # print status message
              if(verbose){
                  memo <- paste("Restricting input to specified genes")
                  message(memo)
              }
              
              # Perform quality checks
              if(class(genes) != 'character')
              {
                  memo <- paste("argument supplied to genes is not of class",
                                 "character, attempting to coerce!")
                  warning(memo)
                  genes <- as.character(genes)
              }
              
              # keep all genes specified and add NA so all samples are kept later
              genes <- c(genes, NA)
              
              return(genes)
          })

#' @rdname Waterfall-methods
#' @aliases Waterfall
#' @param object Object of class waterfall
#' @param mutationHierarchy character vector giving the order of mutations for
#' the hiearchy in order of most to least important
#' @param verbose Boolean for status updates
#' @return data.table object subset based on a mutation hierarchy keeping the
#' most important if there is more than one record for the same gene/sample.
#' @noRd
setMethod(f="mutHierarchySubset",
          signature="data.table",
          definition=function(object, mutationHierarchy, verbose, ...){
              
              # grab the data
              primaryData <- object
              rowsPrimaryData <- nrow(primaryData)
              mutationHierarchy
              
              # refactor the data frame
              primaryData$mutation<- factor(primaryData$mutation, levels=mutationHierarchy$mutation)
              
              # sort the data frame so that the duplicated call will remove the
              # proper mutation
              primaryData <- primaryData[order(primaryData$sample, primaryData$gene, primaryData$mutation),]
              
              # collapse the data on sample/gene
              primaryData <- primaryData[!duplicated(primaryData[, c("sample", "gene")]), ]
              
              # print status message
              if(verbose){
                  memo <- paste("Removed", rowsPrimaryData-nrow(primaryData),
                                "rows when setting the mutation hierarchy.")
                  message(memo)
              }
              
              return(primaryData)
          })

#' @rdname Waterfall-methods
#' @aliases Waterfall
#' @param object Object of class waterfall
#' @param recurrence Numeric value specifying a recurrence cutoff to require,
#' genes not meeting this threshold are removed.
#' @param verbose Boolean for status updates
#' @return Character vector of genes to keep based on recurrence.
#' @noRd
#' @importFrom stats na.omit
setMethod(f="recurrenceSubset",
signature="data.table",
definition=function(object, recurrence, verbose, ...){
    # access the part of the object we want to manipulate
    primaryData <- object
    
    # Dont do anything if recurrence is null
    if(is.null(recurrence)) return(NA)
    
    # Perform quality checks
    if(!is.numeric(recurrence)){
        memo <- paste("argument supplied to recurrence is not of class",
                      "numeric, attempting to coerce!")
        warning(memo)
        recurrence <- as.numeric(as.character(recurrence))
    }
    if(length(recurrence) > 1){
        memo <- paste("argument supplied to recurrence has length > 1",
                      "only the first element was used.")
        warning(memo)
        recurrence <- recurrence[1]
    }
    
    # determine the frequency of gene mutations
    gene <- NULL # appease R CMD CHECK
    mutRecur <- primaryData[, count := .N, by = list(gene)]
    mutRecur <- unique(mutRecur[,c("gene", "count")])
    mutRecur <- stats::na.omit(mutRecur)
    mutRecur$prop <- mutRecur$count/nlevels(primaryData$sample)
    
    # If recurrence cutoff specified exceeds upper limit such that no
    # useful plot would be generated, reset recurrence cutoff
    maxRecur <- max(mutRecur$prop)
    if(maxRecur < recurrence){
        memo <- paste0("The recurrence cutoff specified exceeds the recurrence",
                       " seen in the data, resetting this value to equal max ",
                       "recurrence:", maxRecur)
        warning(memo)
        recurrence <- maxRecur
    }
    
    gene_above_recur <- mutRecur[mutRecur$prop >= recurrence,]$gene
    gene_below_recur <- mutRecur[mutRecur$prop < recurrence,]$gene
    
    # print status message
    if(verbose){
        memo <- paste("Designating for removal", length(unique(gene_below_recur)),
                      "genes not meeting the recurrence cutoff threshold.")
        message(memo)
    }
    
    # add NA to the end of 'gene_above_recurrence' vector, allowing for all
    # samples having NA as a gene name to be retained in subsequent subsets
    gene_above_recur <- c(as.character(gene_above_recur), NA)
    
    return(gene_above_recur)
})

#' @rdname Waterfall-methods
#' @aliases Waterfall
#' @param object Object of class waterfall
#' @param geneOrder Character vector specifying the order in which to plot
#' genes.
#' @param verbose Boolean for status updates
#' @return data.table object genes reordered
#' @noRd
setMethod(f="orderGenes",
          signature="data.table",
          definition=function(object, geneOrder, verbose, ...){
              # access the part of the object we want to manipulate
              primaryData <- object
             
              # print status message
              if(verbose){
                  memo <- paste("Setting gene order")
                  message(memo)
              }
              
              # order based on mutation frequency
              gene_mutation_table <- table(primaryData[,c('gene', 'mutation')])
              geneOrderFreq <- names(sort(rowSums(gene_mutation_table)))
              
              # order the genes based on custom order
              if(!is.null(geneOrder)) {
                  # perform quality checks on geneOrder
                  if(!is.character(geneOrder)){
                      memo <- paste("Argument supplied to gene order is not of",
                                    "class character, attempting to coerce.")
                      geneOrder <- as.character(geneOrder)
                      warning(memo)
                  }
                  if(any(duplicated(geneOrder))){
                      memo <- paste("Detected duplicated element in geneOrder,",
                                    "removing duplicates.")
                      geneOrder <- unique(geneOrder)
                      warning(memo)
                  }
                  
                  # if there are any genes in geneOrder not in x, remove those
                  gene_to_rmv <- geneOrder[!geneOrder %in% unique(primaryData$gene)]
                  if(length(gene_to_rmv) > 0){
                      memo <- paste("The following arguments to geneOrder were",
                                    "not found in the data:", toString(gene_to_rmv))
                      warning(memo)
                      geneOrder <- geneOrder[!geneOrder %in% gene_to_rmv]
                  }
                  if(length(geneOrder) == 0) {
                      memo <- paste0("Found no genes supplied to geneOrder in the",
                                     "data, defaulting to an order based on frequency.")
                      warning(memo)
                  }
                  
                  # merge the custom gene order with the frequency order in case a user did not specify all genes
                  geneOrderFreq <- geneOrderFreq[!geneOrderFreq %in% geneOrder]
                  geneOrder <- c(geneOrderFreq, geneOrder)
                  
              } else {
                  geneOrder <- geneOrderFreq
              }
              
              primaryData$gene <- factor(primaryData$gene, levels=geneOrder)
              return(primaryData)
          })

#' @rdname Waterfall-methods
#' @aliases Waterfall
#' @param object Object of class waterfall
#' @param geneMax Integer specifying the maximum number of genes to be plotted.
#' @param verbose Boolean for status updates
#' @return data.table object subset to contain only the max number of genes.
#' @noRd
setMethod(f="maxGeneSubset",
          signature="data.table",
          definition=function(object, geneMax, verbose, ...){
              # access the part of the object we want to manipulate
              primaryData <- object

              # do nothing if null
              if(is.null(geneMax)){
                  return(primaryData)
              }
              
              # perform quality checks
              if(!is.numeric(geneMax)){
                  memo <- paste("geneMax is not numeric, attempting to convert",
                                "to an integer.")
                  geneMax <- as.integer(geneMax)
                  warning(memo)
              }
              
              if(geneMax %% 1 != 0){
                  memo <- paste("geneMax is not a whole number, rounding.")
                  geneMax <- round(geneMax)
                  warning(memo)
              }
              
              # limit to max genes
              keepGenes <- utils::tail(levels(primaryData$gene), geneMax)
              removeGenes <- unique(primaryData[!primaryData$gene %in% keepGenes,"gene"])
              primaryData <- primaryData[primaryData$gene %in% keepGenes,]
              
              # print status message
              if(verbose){
                  memo <- paste("geneMax is set to", geneMax, "removing",
                                nrow(removeGenes),"genes")
                  message(memo)
              }
              
              return(primaryData)
          })

#' @rdname Waterfall-methods
#' @aliases Waterfall
#' @param object Object of class waterfall
#' @param sampOrder Character vector specifying the order of samples
#' @param verbose Boolean for status updates
#' @return data.table object with samples reordered
#' @importFrom data.table dcast
#' @noRd
setMethod(f="orderSamples",
          signature="data.table",
          definition=function(object, sampleOrder, verbose, ...){
              # access the part of the object we want to manipulate
              primaryData <- object
              
              # print status message
              if(verbose){
                  memo <- paste("setting the order of samples.")
                  message(memo)
              }
              
              # if sampleOrder is specified reorder on that
              if(!is.null(sampleOrder)){
                  # perform quality checks
                  if(!is.character(sampleOrder)){
                      memo <- paste("sampleOrder is not a character vector,",
                                    "attempting to coerce.")
                      warning(memo)
                      if(is.list(sampleOrder)){
                          sampleOrder <- unlist(sampleOrder)
                      }
                      sampleOrder <- as.character(sampleOrder)
                  }
                  if(sum(duplicated(sampleOrder)) != 0){
                      memo <- paste("Found duplicate elements in sampleOrder, uniquing!")
                      warning(memo)
                      sampleOrder <- unique(sampleOrder)
                  }
                  
                  # check if there are any samples specified not in the data
                  newSamples <- sampleOrder[!sampleOrder %in% levels(primaryData$sample)]
                  if(length(newSamples) != 0){
                      memo <- paste("The following samples were not detected",
                                    "in the data or its subsets:",
                                    toString(newSamples),
                                    ". Binding these to the data.")
                      warning(memo)
                      newSampleDT <- data.table::data.table("sample"=newSamples, "gene"=NA, "mutation"=NA, "label"=NA)
                      primaryData <- rbind(primaryData, newSampleDT, fill=TRUE)
                  }
                  
                  # status message
                  removeSample <- unique(levels(primaryData$sample)[!levels(primaryData$sample) %in% sampleOrder])
                  if(length(removeSample) != 0){
                      memo <- paste("Removing the samples:", toString(removeSample),
                                    "which were found in the data but not sampleOrder")
                      warning(memo)
                      primaryData <- primaryData[!primaryData$sample %in% removeSample,]
                  }
                  
                  # return the data with reordered samples
                  primaryData$sample <- factor(primaryData$sample, levels=sampleOrder)
                  return(primaryData)
              }
              
              # perform a hierarchical sort if sampleOrder is null
              # recast the data going from long format to wide format,
              #values in this data are counts of a mutation call
              wide_data <- data.table::dcast(primaryData, sample ~ gene,
                                             fun.aggregate = length, value.var="mutation")
              
              # apply a boolean function to convert the data frame values to 1's and 0's
              values <- wide_data[,-1, drop=FALSE]
              sample <- wide_data[,1]
              values <- data.frame(apply(values, 2,
                                         function(x) as.numeric(as.logical(x))))
              wideBoolean <- cbind(sample, values)
              
              # reverse the columns so that genes with highest mutation's are 
              # listed first (assumes gene_sort has been run on the data frame)
              wideBoolean <- wideBoolean[,c(1, rev(2:ncol(wideBoolean))), with=FALSE]
              
              # if there are any NA values present in a sample at the gene put that
              # remove that sample and save (at this stage it should be samples)
              if(any(grepl("^NA.$", colnames(wideBoolean)))) {
                  # Find which column has the NA header
                  NA_index <- which(grepl("^NA.$", colnames(wideBoolean)))
                  
                  # Append NA column to end of data frame
                  NA_gene <- wideBoolean[,NA_index, with=FALSE]
                  wideBoolean <- wideBoolean[,-NA_index, with=FALSE]
                  
                  # Save copy and remove samples with no mutations,
                  # these will be added to the end
                  samp_no_mut <- wideBoolean[rowSums(wideBoolean[, 2:ncol(wideBoolean)]) == 0,]$sample
                  samp_no_mut <- as.character(samp_no_mut)
                  wideBoolean <- wideBoolean[!wideBoolean$sample %in% samp_no_mut,]
              } else {
                  samp_no_mut <- NULL
              }
              
              # hiearchial sort on all column's (i.e. genes) such that samples are
              # rearranged if there is a mutation in that gene
              sampleOrder <- wideBoolean[do.call(order, as.list(-wideBoolean[,2:ncol(wideBoolean), with=F])),]$sample
              
              # Put those samples not in sample order in from the original levels of the
              # data (these are samples with no mutations)
              not_in <- as.character(levels(sampleOrder)[which(!levels(sampleOrder) %in% sampleOrder)])
              not_in <- not_in[!not_in %in% samp_no_mut]
              sampleOrder <- c(as.character(sampleOrder), as.character(not_in))
              
              # Put those samples with no mutations back in
              if(!is.null(samp_no_mut)){
                  sampleOrder <- c(sampleOrder, samp_no_mut)
              }
              
              # set the new sample order
              primaryData$sample <- factor(primaryData$sample, levels=unique(sampleOrder))
              return(primaryData)
          })

#' @rdname Waterfall-methods
#' @aliases Waterfall
#' @param object Object of class waterfall
#' @param plotA String specifying the type of plot for the top sub-plot, one of
#' "burden", "frequency", or NULL for a mutation burden (requires coverage to be
#' specified), frequency of mutations, or no plot respectively.
#' @param plotATally String specifying one of "simple" or "complex" for a
#' simplified or complex tally of mutations respectively.
#' @param plotALayers list of ggplot2 layers to be passed to the plot.
#' @param verbose Boolean for status updates
#' @return gtable object containing the top sub-plot.
#' @noRd
#' @import ggplot2
#' @importFrom gtable gtable
setMethod(f="buildMutationPlot",
          signature="WaterfallData",
          definition=function(object, plotA, plotATally, plotALayers, verbose, ...){
              
              # grab only the first element for parameters
              plotA <- plotA[1]
              plotATally <- plotATally[1]
              
              # if plot type is null return an empty gtable
              if(is.null(plotA)) return(gtable::gtable())
              
              # print status message
              if(verbose) {
                  memo <- paste("Constructing top sub-plot")
                  message(memo)
              }
              
              # make sure plotATally is valid
              if(!toupper(plotATally) %in% toupper(c("simple", "complex"))){
                  memo <- paste("plotATally is not set to either \"simple\" or",
                                " \"complex\"... defaulting to \"simple\"")
                  warning(memo)
                  plotATally <- "simple"
              }
              
              # extract the data for the type of plot we need
              if(toupper(plotATally) == toupper("simple")){
                  mutationData <- getData(object, name="simpleMutationCounts")
              } else if(toupper(plotATally) == toupper("complex")) {
                  mutationData <- getData(object, name="complexMutationCounts")
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

              if(!is.null(plotA) && !toupper(plotA) %in% toupper(c("frequency", "burden"))){
                  memo <- paste("plotA is not set to \"frequency\", \"burden\"",
                                " or NULL... defaulting to \"frequency\".")
                  warning(memo)
                  plotA <- "frequency"
              }
              if(all(is.na(mutationData$mutationBurden)) && toupper(plotA) == toupper("burden")){
                  memo <- paste("plotA is set to:",toString(plotA),"but could",
                                "not find calculated mutation burden, please",
                                "specify coverage!, Resetting plotA to frequency.")
                  warning(memo)
                  plotA <- "frequency"
              }
              
              # make sure sample levels match primaryData for plotting
              mutationData$sample <- factor(mutationData$sample, levels=levels(getData(object, name="primaryData")$sample))
              if(toupper(plotATally) == toupper("complex")) {
                  mutationData$mutation <- factor(mutationData$mutation, levels=levels(getData(object, name="primaryData")$mutation))
              } else if(toupper(plotATally) == toupper("simple")) {
                  # Do Nothing
              }
              
              ############# set ggplot2 layers #################################
              
              # theme
              plotTheme <- theme(axis.ticks.x = element_blank(),
                                 axis.text.x = element_blank(),
                                 axis.title.x = element_blank(),
                                 legend.title = element_text(size=14),
                                 axis.text.y = element_text(colour = "black"),
                                 axis.title.y = element_text(colour = "black"),
                                 panel.background = element_blank(),
                                 # panel.grid.minor.y = element_line(colour = "black"),
                                 panel.grid.major.y = element_line(colour = "grey80"),
                                 panel.grid.minor.x = element_blank(),
                                 panel.grid.major.x = element_blank(),
                                 panel.border = element_rect(fill = NA)
              )
              
              # legend
              if(toupper(plotATally) == toupper("simple")){
                  plotLegend <- scale_fill_manual(name="Translational Effect",
                                                  na.value="deepskyblue4", values="deepskyblue4")
              } else if(toupper(plotATally) == toupper("complex")){
                  plotLegend <- scale_fill_manual(name="Translational Effect",
                                                  values=getData(object, name="mutationHierarchy")$color,
                                                  breaks=getData(object, name="mutationHierarchy")$mutation,
                                                  drop=FALSE)
                  plotTheme <- plotTheme + theme(legend.position="none")
              }
              
              # titles
              if(toupper(plotA) == toupper("frequency")) {
                  plotTitleY <- ylab("Mutation\nFrequency")
              } else if(toupper(plotA) == toupper("burden")) {
                  plotTitleY <- ylab("Mutations\nper MB")
              }
              
              # keep all samples
              x_scale <- scale_x_discrete(drop=FALSE)
              
              #  geom definition
              plotGeom <- geom_bar(stat='identity', alpha=1, width=1)
              
              # plot
              if(toupper(plotA) == toupper("frequency")) {
                  mutPlot <- ggplot(mutationData, aes_string(x='sample', y='Freq', fill='mutation'))
              } else if(toupper(plotA) == toupper("burden")) {
                  mutPlot <- ggplot(mutationData, aes_string(x='sample', y='mutationBurden', fill='mutation'))
              }
              mutPlot <- mutPlot + plotGeom + x_scale + plotTitleY + plotLegend + plotTheme + plotALayers
              
              # convert to gtable grob
              plotGrob <- ggplotGrob(mutPlot)
              return(plotGrob)
          })

#' @rdname Waterfall-methods
#' @aliases Waterfall
#' @param object Object of class waterfall
#' @param verbose Boolean for status updates
#' @return data.table object containing summarized gene level data
#' @noRd
setMethod(f="constructGeneData",
          signature="data.table",
          definition=function(object, verbose, ...){
              # extract the data to work with
              primaryData <- object
              
              # status message
              if(verbose){
                  memo <- paste("Constructing GeneData")
                  message(memo)
              }
              
              # construct geneData
              gene <- mutation <- NULL # appeases R CMD CHECK
              geneData <- primaryData[, count := .N, by = list(gene, mutation)]
              geneData <- unique(geneData[,c("gene", "mutation", "count")])
              
              # if a gene is present but has NA for mutation this probably means
              # you want to plot it but it's count needs to be 0
              geneData[is.na(geneData$mutation),"count"] <- 0
              
              # set the levels
              geneData$gene <- factor(geneData$gene, levels=levels(primaryData$gene))
              
              return(geneData)
          })

#' @rdname Waterfall-methods
#' @aliases Waterfall
#' @param object Object of class waterfall
#' @param plotB String specifying the type of plot for the left sub-plot, one of
#' "proportion", "frequency", or NULL for a plot of gene proportions frequencies
#' , or no plot respectively.
#' @param plotBTally String specifying one of "simple" or "complex" for a
#' simplified or complex tally of genes respectively.
#' @param plotBLayers list of ggplot2 layers to be passed to the plot.
#' @param verbose Boolean for status updates
#' @return gtable object containing the left sub-plot.
#' @noRd
#' @import ggplot2
#' @importFrom gtable gtable
setMethod(f="buildGenePlot",
          signature="WaterfallData",
          definition=function(object, plotB, plotBTally, plotBLayers, verbose, ...){
              # grab only the first element for parameters
              plotB <- plotB[1]
              plotBTally <- plotBTally[1]
              
              # if plot type is null return an empty gtable
              if(is.null(plotB)) return(gtable::gtable())
              
              # extract the data needed for this plot
              geneData <- getData(object, name="geneData")
              
              # print status message
              if(verbose) {
                  memo <- paste("Constructing left sub-plot")
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

              if(!is.null(plotB) && !toupper(plotB) %in% toupper(c("frequency", "proportion"))){
                  memo <- paste("plotB is not set to \"frequency\", \"proportion\"",
                                " or NULL... defaulting to \"proportion\".")
                  warning(memo)
                  plotB <- "proportion"
              }
              if(!toupper(plotBTally) %in% toupper(c("simple", "complex"))){
                  memo <- paste("plotBTally is not set to either \"simple\" or",
                                " \"complex\"... defaulting to \"simple\"")
                  warning(memo)
                  plotBTally <- "simple"
              }
              
              # determine proportion of gene mutations
              sampleCount <- nlevels(getData(object, name="primaryData")$sample)
              geneData$proportion <- geneData$count/sampleCount * 100
              
              ################ ggplot2 #########################################
              
              # theme
              plotTheme <- theme(axis.text.y=element_text(colour='black', face='italic'),
                                 axis.title.y=element_blank(),
                                 legend.position=('none'), panel.grid.major.y=element_blank(),
                                 panel.grid.minor.y=element_blank())
              # titles
              if(plotB == "frequency"){
                  plotYlabel <- ylab('# Mutant')
              }else if(plotB == "proportion"){
                  plotYlabel <- ylab('% Mutant')
              }
              
              # legend
              if(plotBTally == "simple"){
                  # if simple we have to reset mutations to NA
                  geneData$mutation <- "NA"
                  plotLegend <- scale_fill_manual(name="Translational Effect",
                                                  na.value="deepskyblue4", values="deepskyblue4")
              }else if(plotBTally == "complex"){
                  plotLegend <- scale_fill_manual(name="Translational Effect",
                                                  values=getData(object, name="mutationHierarchy")$color,
                                                  breaks=getData(object, name="mutationHierarchy")$mutation,
                                                  drop=FALSE)
              }
              
              # geom definition
              plotGeom <- geom_bar(position='stack', alpha=1, width=1, stat='identity')
              
              # plot
              if(plotB == "proportion"){
                  genePlot <- ggplot(geneData, aes_string(x='gene', y='proportion', fill='mutation'))
              }else if(plotB == "frequency"){
                  genePlot <- ggplot(geneData, aes_string(x='gene', y='count', fill='mutation'))
              }
              
              genePlot <- genePlot + plotGeom + theme_bw() + coord_flip() + plotTheme +
                  plotYlabel + scale_y_reverse() + plotLegend + plotBLayers
              
              plotGrob <- ggplotGrob(genePlot)
          })

#' @rdname formatClinicalData-methods
#' @aliases Waterfall
#' @noRd
#' @importFrom data.table setDT
#' @importFrom data.table is.data.table
#' @importFrom data.table melt
#' @importFrom data.table data.table
#' @importFrom data.table rbindlist
setMethod(f="formatClinicalData",
          signature="WaterfallData",
          definition=function(object, clinical, verbose, ...){
              
              # extract the data we need
              clinicalData <- getData(clinical)
              primaryData <- getData(object, name="primaryData")
              
              # print status message
              if(verbose){
                  memo <- paste("Formatting clinical data")
                  message(memo)
              }
              
              # remove clinical samples not found in the primary data
              primaryDataSamples <- levels(primaryData$sample)
              removedSamp <- unique(clinicalData$sample[!clinicalData$sample %in% primaryDataSamples])
              clinicalData <- clinicalData[clinicalData$sample %in% primaryDataSamples,]
              memo <- paste("Removed", length(removedSamp), "samples from the clinical data",
                            "not found in the primary waterfall data! These samples are:",
                            toString(removedSamp))
              if(length(removedSamp) != 0){
                  warning(memo)
              } else if(verbose){
                  message(memo)
              }
              
              # fill in missing clinical samples from primary data if necessary
              fillSamp <- unique(primaryDataSamples[!primaryDataSamples %in% clinicalData$sample])
              clinList <- list(clinicalData, data.table::data.table("sample"=fillSamp))
              clinicalData <- data.table::rbindlist(clinList, use.names=TRUE,
                                                    fill=TRUE)
              memo <- paste("Added", length(fillSamp), "samples from the primary data",
                            "to the clinical data! These samples are:", toString(fillSamp))
              if(length(fillSamp) != 0){
                  warning(memo)
              }
              
              # set levels of clinicalData to match primaryData for ordering
              clinicalData$sample <- factor(clinicalData$sample, levels=levels(primaryData$sample))
              clinicalData$value <- factor(clinicalData$value, levels=unique(clinicalData$value))
              
              # return the formated data
              return(clinicalData)
          })

#' @rdname Waterfall-methods
#' @aliases Waterfall
#' @param object Object of class waterfall
#' @param verbose Boolean for status updates
#' @param gridOverlay Boolean specifying if a grid should be overlayed on the
#' waterfall plot.
#' @param drop Boolean specifying if unused mutations should be dropped from the
#' legend.
#' @param plotCLayers list of ggplot2 layers to be passed to the plot.
#' @return gtable object containing the main plot.
#' @noRd
#' @import ggplot2
#' @importFrom gtable gtable
setMethod(f="buildWaterfallPlot",
          signature="WaterfallData",
          definition=function(object, gridOverlay, drop, labelSize, labelAngle,
                              sampleNames, xTitle, plotCLayers, verbose, ...){
              # extract the data we need
              primaryData <- getData(object, name="primaryData")
              paletteData <- getData(object, name="mutationHierarchy")
              mutationData <- getData(object, name="complexMutationCounts")
              
              # subset all mutation data to just the mutations in mutationData
              # this should contain everything in the primaryData and other plots
              # and is the staring point for the legend, using the param drop will
              # remove legend entries for those not in the primary data at all
              mutationData <- as.character(unique(mutationData[mutationData$Freq > 0,]$mutation))
              paletteData <- paletteData[paletteData$mutation %in% mutationData,]
              
              # related to adding support for genes with no mutations, we need to assign
              # any NA samples to an actual sample so the gene with no mutation is plotted
              primaryData$sample[is.na(primaryData$sample)] <- levels(primaryData$sample)[1]
              
              # this should keep the hierarchy but remove unused levels based on mutationData
              primaryData$mutation <- factor(primaryData$mutation,
                                             levels=levels(primaryData$mutation)[levels(primaryData$mutation) %in% paletteData$mutation])
              
              # there could be samples with no gene information assign these a gene but no
              # mutation so they are plotted correctly
              primaryData[is.na(primaryData$gene),"gene"] <- tail(levels(primaryData$gene))[1]
              
              # print status message
              if(verbose){
                  memo <- paste("Building the main plot")
                  message(memo)
              }
              
              # perform quality checks
              if(!is.null(plotCLayers)){
                  if(!is.list(plotCLayers)){
                      memo <- paste("plotCLayers is not a list")
                      stop(memo)
                  }
                  
                  if(any(!unlist(lapply(plotCLayers, function(x) ggplot2::is.ggproto(x) | ggplot2::is.theme(x) | is(x, "labels"))))){
                      memo <- paste("plotCLayers is not a list of ggproto or ",
                                    "theme objects... setting plotCLayers to NULL")
                      warning(memo)
                      plotCLayers <- NULL
                  }
              }
              
              ######### start building the plot ################################
              # grid overlay
              if(gridOverlay){
                  verticalPlotGrid <- geom_vline(xintercept = seq(1.5, nlevels(primaryData$sample),
                                                                  by=1),
                                                 linetype='solid', colour='grey40', size=.1)
                  
                  if(length(unique(primaryData$gene)) == 1)
                  {
                      horizontalPlotGrid <- geom_blank()
                  } else {
                      horizontalPlotGrid <- geom_hline(yintercept = seq(1.5, length(unique(primaryData$gene)),
                                                                        by=1),
                                                       linetype='solid', colour='grey40',
                                                       size=.1)
                  }
                  plotGridOverlay <- list(horizontalPlotGrid, verticalPlotGrid)
              } else {
                  plotGridOverlay <- geom_blank()
              }
              
              # construct legend
              if(drop)
              {
                  valueMap <- paletteData$color
                  names(valueMap) <- paletteData$mutation
                  plotLegend <- scale_fill_manual(name="Mutation Type",
                                                  values=valueMap,
                                                  breaks=paletteData$mutation,
                                                  labels=paletteData$label,
                                                  drop=TRUE)
              } else {
                  valueMap <- paletteData$color
                  names(valueMap) <- paletteData$mutation
                  plotLegend <- scale_fill_manual(name="Mutation Type",
                                                  values=valueMap,
                                                  breaks=paletteData$mutation,
                                                  labels=paletteData$label,
                                                  drop=FALSE)
              }
              
              # plot titles
              plotXLabel <- xlab(paste0('Sample (n=', nlevels(primaryData$sample), ')'))
              
              if(all(is.na(primaryData$label)))
              {
                  label <- geom_blank()
              } else {
                  label <- geom_text(data=primaryData,
                                     mapping=aes_string(x='sample', y='gene',
                                                        label='label'),
                                     size=labelSize, colour='white',
                                     angle=labelAngle)
              }
              
              # base theme
              plotTheme <- theme(axis.ticks=element_blank(),
                                 panel.grid.major=element_blank(),
                                 panel.grid.minor=element_blank(),
                                 axis.title.y=element_blank(),
                                 panel.background=element_rect(fill='white',
                                                               colour='white'),
                                 axis.text.y=element_blank(),
                                 plot.title=element_blank(),
                                 legend.title=element_text(size=14),
                                 panel.border = element_rect(colour = "grey20",
                                                             fill=NA))
              
              if(sampleNames){
                  sampleLabels <- theme(axis.text.x=element_text(angle=50, hjust=1))
              } else {
                  sampleLabels <- theme(axis.text.x=element_blank(),
                                        axis.title.x=element_blank())
              }
              if(xTitle){
                  plotXtitle <- theme(axis.title.x=element_text(size=20))
              } else {
                  plotXtitle <- theme(axis.title.x=element_blank())
              }
              
              # make sure all samples are plotted
              x_scale <- scale_x_discrete(drop=FALSE)
              
              # define geom
              plotGeom <- geom_tile(aes_string(fill='mutation'), position="identity")
              
              # define plot
              waterfallPlot <- ggplot(primaryData, aes_string('sample', 'gene'))
             
              # combine plot elements
              waterfallPlot <- waterfallPlot + plotGeom + plotGridOverlay + x_scale +
                  plotLegend + plotXLabel + label + plotTheme + sampleLabels + plotXtitle + plotCLayers
              
              # covert to grob
              waterfallGrob <- ggplotGrob(waterfallPlot)
              return(waterfallGrob)
          })

#' @rdname arrangeWaterfallPlot-methods
#' @aliases arrangeWaterfallPlot
#' @param sectionHeights Relative heights of each plot section (should sum to one).
#' @param sectionWidths Relative widths of each plot section (should sum to one).
#' @noRd
#' @importFrom grid nullGrob
#' @importFrom gridExtra arrangeGrob
setMethod(f="arrangeWaterfallPlot",
          signature="WaterfallPlots",
          definition=function(object, sectionHeights, sectionWidths, verbose, ...){
              
              # grab the data we need
              plotA <- getGrob(object, 1)
              plotB <- getGrob(object, 2)
              plotC <- getGrob(object, 3)
              plotD <- getGrob(object, 4)
              
              # set default widths
              
              # Future Zach, what were you thinking TODO fix this hack
              if(length(plotB) == 0){
                  sectionWidths <- c(0, 1)
              } else {
                  if(is.null(sectionWidths)){
                      sectionWidths <- c(.20, .80)
                  } else if(length(sectionWidths) != 2){
                      memo <- paste("sectionWidths should be of length 2... Using",
                                    "default values!")
                      warning(memo)
                      sectionWidths <- c(.20, .80)
                  } else if(!all(is.numeric(sectionWidths))) {
                      memo <- paste("sectionWidths must be numeric... Using",
                                    "default values!")
                      warning(memo)
                      sectionWidths <- c(.20, .80)
                  }
              }
              
              # set up a blank plot for alignment purposes
              emptyPlot <- grid::nullGrob()
              
              ## Strip out legends and plot separately
              ## https://github.com/baptiste/gridextra/wiki/arranging-ggplot#legends
              #ind_legend <- grep("guide", plotC$layout$name)
              #plotC_legend <- plotC[["grobs"]][[ind_legend]]
              #plotC_width <- sum(plotC_legend$width)
              #heatmap_grob <- ggplot2::ggplotGrob(grid.draw(plotC) + theme(legend.position="none"))
              
              # align the heights of plotB and plotC if they exist
              if(length(plotB) != 0){
                  maxHeight <- grid::unit.pmax(plotB$heights, plotC$heights)
                  plotB$heights <- as.list(maxHeight)
                  plotC$heights <- as.list(maxHeight)
              }
              
              # Obtain the max width for relevant plots
              plotList4Width <- list(plotA, plotC, plotD)
              plotList4Width <- plotList4Width[lapply(plotList4Width, length) > 0]
              plotList4Width <- lapply(plotList4Width, function(x) x$widths)
              maxWidth <- do.call(grid::unit.pmax, plotList4Width)
              
              ################# Adjust the x-label of the gene plot ############
              # Find the necessary height adjustment
              if(length(plotB) != 0){
                  axis_b_mainPlot_Height_index <- plotC$layout[which(plotC$layout$name == "axis-b"),"t"]
                  axis_b_mainPlot_Height <- plotC$heights[axis_b_mainPlot_Height_index]
                  axis_b_genePlot_Height_index <- getGrob(object, 2)$layout[which(getGrob(object, 2)$layout$name == "axis-b"),"t"]
                  axis_b_genePlot_Height <- getGrob(object, 2)$heights[axis_b_genePlot_Height_index]
                  
                  heightAdjustment <- axis_b_mainPlot_Height - axis_b_genePlot_Height
                  
                  # set the gene plot back axis-b back to it's original value
                  plotB$heights[axis_b_genePlot_Height_index] <- getGrob(object, 2)$heights[axis_b_genePlot_Height_index]
                  
                  # add a row to the gene plot below the  the x-axis title of the appropriate adjustment
                  xlab_b_genePlot_Height_index <- plotB$layout[which(plotB$layout$name == "xlab-b"),"t"]
                  plotB <- gtable::gtable_add_rows(plotB, heights=heightAdjustment, pos=xlab_b_genePlot_Height_index)
              }
              
              ############ section plots into rows and set relative widths ###################
              
              # section plots into rows, also define preliminary heights
              plotList <- list()
              plotHeight <- list()
              # set width for first row of plots
              if(length(plotA) != 0){
                  plotA$widths <- as.list(maxWidth)
                  firstRowPlots <- gridExtra::arrangeGrob(emptyPlot, plotA, ncol=2, widths=sectionWidths)
                  plotList[["firstRow"]] <- firstRowPlots
                  plotHeight[["firstRow"]] <- .2
              }
              
              # set widths for second row of plots
              plotC$widths <- as.list(maxWidth)
              if(length(plotB) != 0){
                  secondRowPlots <- gridExtra::arrangeGrob(plotB, plotC, ncol=2, widths=sectionWidths)
              } else {
                  secondRowPlots <- gridExtra::arrangeGrob(emptyPlot, plotC, ncol=2, widths=sectionWidths)
              }
              plotList[["secondRow"]] <- secondRowPlots
              plotHeight[["secondRow"]] <- .8
              
              # set widths for third row of plots
              if(length(plotD) != 0){
                  plotD$widths <- as.list(maxWidth)
                  lastRowPlots <- gridExtra::arrangeGrob(emptyPlot, plotD, ncol=2, widths=sectionWidths)
                  plotList[["lastRow"]] <- lastRowPlots
                  plotHeight[["lastRow"]] <- .2
              }
              
              # set section heights based upon the number of sections
              # TODO another hack, make this more like auto-magic
              plotHeight <- as.vector(do.call(rbind, plotHeight))
              if(length(plotHeight) == 3){
                  plotHeight <- c(.2, .6, .2)
              }
              
              if(is.null(sectionHeights)){
                  sectionHeights <- plotHeight
              } else if(length(sectionHeights) != length(plotList)){
                  memo <- paste("There are", length(sectionHeights), "provided",
                                "but", length(plotList), "vertical sections...",
                                "using default values!")
                  warning(memo)
                  sectionHeights <- plotHeight
              } else if(!all(is.numeric(sectionHeights))) {
                  memo <- paste("sectionHeights must be numeric... Using",
                                "default values!")
                  warning(memo)
                  sectionHeights <- plotHeight
              }
              
              # arrange the final plot
              # http://stackoverflow.com/questions/6681145/how-can-i-arrange-an-arbitrary-number-of-ggplots-using-grid-arrange
              finalPlot <- do.call(gridExtra::arrangeGrob, c(plotList, list(ncol=1, heights=sectionHeights)))
              
              return(finalPlot)
          })

#' @rdname geneFilter-methods
#' @aliases geneFilter
#' @param genes Character vector of genes to keep based on geneSubset and recurrenceSubset
#' @noRd
setMethod(
    f="geneFilter",
    signature="data.table",
    definition=function(object, genes, verbose, ...){
        
        # extract the data we need
        primaryData <- object
        
        # do nothing if there are no genes to remove
        if(all(is.na(genes))) return(primaryData)
        
        if(verbose){
            memo <- paste("Performing gene subsets")
            message(memo)
        }
        
        # check if a gene was specified but is not in the data
        if(!all(na.omit(genes) %in% primaryData$gene))
        {
            newGenes <- na.omit(genes[!genes %in% primaryData$gene])
            memo <- paste("The following genes were designated to be kept but were",
                          "not found in the data:", toString(newGenes), "adding these to the data!")
            warning(memo)
            newGenes <- data.table::data.table("sample"=NA, "gene"=newGenes, "mutation"=NA, "label"=NA)
            primaryData <- list(primaryData, newGenes)
            primaryData <- data.table::rbindlist(primaryData, use.names=TRUE,
                                                 fill=TRUE)
        }
        
        # subset the original data frame based on the following: keep gene if it is
        # in the gene vector in "mutation_recurrence_subset"
        primaryData <- primaryData[(primaryData$gene %in% genes), ]
        
        return(primaryData)
    }
)