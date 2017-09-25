#' Class MutSpectra
#' 
#' An S4 class for the MutSpectra plot object
#' @name MutSpectra-class
#' @rdname MutSpectra-class
#' @slot PlotA gtable object for the mutation frequencies.
#' @slot PlotB gtable object for the mutation proportions.
#' @slot PlotC gtable object for clinical data sub-plot.
#' @slot Grob gtable object for the arranged plot.
#' @slot primaryData data.table object storing the primary data, should have
#' column names sample, mutation, frequency, proportion.
#' @slot ClinicalData data.table object storing the data used to plot the
#' clinical sub-plot.
#' @exportClass MutSpectra
#' @import methods
#' @importFrom gtable gtable
#' @importFrom data.table data.table

methods::setOldClass("gtable")
setClass("MutSpectra",
         representation=representation(PlotA="gtable",
                                       PlotB="gtable",
                                       PlotC="gtable",
                                       Grob="gtable",
                                       primaryData="data.table",
                                       ClinicalData="data.table"),
         validity=function(object){
             cat("!!!!! MutSpectra~Inspector !!!!!\n")
         }
)

#' Initalizer method for the MutSpectra class
#' 
#' @name MutSpectra
#' @rdname MutSpectra-class
#' @param .Object object of class MutSpectra
#' @noRd
#' @import ggplot2
setMethod(f="initialize",
          signature="MutSpectra",
          definition=function(.Object, input, BSgenome, sorting, palette, clinical, sectionHeights,
                              sectionWidths, sampleNames, verbose, plotALayers, plotBLayers,
                              plotCLayers){
              # convert object to mutSpectra format
              .Object@primaryData <- toMutSpectra(input, BSgenome=BSgenome, verbose=verbose)

              # annotate types of transitions and transversions
              .Object@primaryData <- annoMutSpectra(.Object, verbose)
              
              # calculate rates of transitions and transversions
              .Object@primaryData <- calcMutSpectra(.Object, verbose)
              
              # sort the samples for plotting
              .Object@primaryData <- sortSamples(.Object, sorting ,verbose)
              
              # add the clinical data
              if(is.null(clinical)){
                  .Object@ClinicalData <- data.table::data.table()
              } else {
                  # fill in the clinical data slot
                  .Object@ClinicalData <- getData(clinical)
                  # format the clinical data object using the mutSpectra object
                  .Object@ClinicalData <- formatClinicalData(.Object, verbose)
                  # add layer to suppress the x-axis text in the proportion plot
                  addLayer <- theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
                  plotBLayers[[length(plotBLayers ) + 1]] <- addLayer
              }
              
              # add the clinical data plot
              .Object@PlotC <- buildClinicalPlot(.Object, clinicalLayers=clinical@clinicalLayers)
              
              # construct the frequency plot
              .Object@PlotA <- buildFrequencyPlot(.Object, plotALayers, palette, verbose)
              
              # construct the proportion plot
              .Object@PlotB <- buildProportionPlot(.Object, sampleNames, plotBLayers, palette, verbose)
              
              # align all the plots together
              .Object@Grob <- arrangeMutSpectraPlot(.Object, sectionHeights, verbose)
              
              return(.Object)
          })

#' Constructor for the MutSpectra class.
#' 
#' @name MutSpectra
#' @rdname MutSpectra-class
#' @param input Object of class MutationAnnotationFormat, GMS, VEP.
#' @param BSgenome Object of class BSgenome, used to extract reference bases if
#' not supplied by the file format.
#' @param sorting Character vector specifying how samples should be ordered in the plot, one
#' of "mutation", "sample", or a vector of length equal to the number of samples explicitly 
#' providing the order of samples.
#' @param palette Character vector specifying the colors used for encoding transitions and transversions
#' , should be of length 6. If NULL a default palette will be used.
#' @param clinical Object of class Clinical, used for adding a clinical data subplot.
#' @param sectionHeights Numeric vector specifying relative heights of each plot section,
#' should sum to one. Expects a value for each section.
#' @param sectionWidths Numeric vector specifying relative heights of each plot section,
#' should sum to one. Expects a value for each section.
#' @param sampleNames Boolean specifying if samples should be labeled on the plot.
#' @param verbose Boolean specifying if status messages should be reported
#' @param plotALayers list of ggplot2 layers to be passed to the frequency plot.
#' @param plotBLayers list of ggplot2 layers to be passed to the proportion plot.
#' @param plotCLayers list of ggplot2 layers to be passed to the clinical plot.
#' @export
MutSpectra <- function(input, BSgenome=NULL, sorting=NULL, palette=NULL, clinical=NULL, sectionHeights=NULL,
                      sectionWidths=NULL, sampleNames=TRUE, verbose=FALSE, plotALayers=NULL, plotBLayers=NULL,
                      plotCLayers=NULL){
    cat("!!!!! MutSpectra~Constructor !!!!!\n")
    new("MutSpectra", input=input, BSgenome=BSgenome, sorting=sorting, palette=palette, clinical=clinical, sectionHeights=sectionHeights,
        sectionWidths=sectionWidths, sampleNames=sampleNames, verbose=verbose, plotALayers=plotALayers, plotBLayers = plotBLayers, 
        plotCLayers=plotCLayers)
}

#' @rdname MutSpectra-methods
#' @aliases annoMutSpectra,MutSpectra
#' @param object Object of class MutSpectra
#' @param verbose Boolean for status updates
#' @return data.table object with transitions and transversions annotated
#' @noRd
setMethod(f="annoMutSpectra",
          signature="MutSpectra",
          definition=function(object, verbose, ...){
              # access the part of the object we want to manipulate
              primaryData <- object@primaryData
              
              # print status message
              if(verbose){
                  memo <- paste("annotating mutations as transitions and transversions")
                  message(memo)
              }
              
              # add an extra column with the reference to variant
              primaryData$base_change <- paste0(primaryData$refAllele, "2", primaryData$variantAllele)
              
              # annotate the grouping of the base change
              getMutChange <- function(x){
                  switch(x, A2C="A->C or T->G (TV)",
                         T2G="A->C or T->G (TV)", A2G="A->G or T->C (TI)",
                         T2C="A->G or T->C (TI)", A2T="A->T or T->A (TV)",
                         T2A="A->T or T->A (TV)", G2A="G->A or C->T (TI)",
                         C2T="G->A or C->T (TI)", G2C="G->C or C->G (TV)",
                         C2G="G->C or C->G (TV)", G2T="G->T or C->A (TV)",
                         C2A="G->T or C->A (TV)")
              }
              primaryData$trans_tranv <- sapply(primaryData$base_change, getMutChange)
              
              # make the transitions and transverions annotated a factor of all possible outcomes
              primaryData$trans_tranv <- factor(primaryData$trans_tranv,
                                                levels=c("A->C or T->G (TV)", "A->G or T->C (TI)",
                                                         "A->T or T->A (TV)", "G->A or C->T (TI)",
                                                         "G->C or C->G (TV)", "G->T or C->A (TV)"))
              
              # remove the temp base change column
              primaryData$base_change <- NULL
              return(primaryData)
          })

#' @rdname MutSpectra-methods
#' @aliases calcMutSpectra,MutSpectra
#' @param object Object of class MutSpectra
#' @param verbose Boolean for status updates
#' @return data.table object with transitions and transversions rates calculated
#' @importFrom data.table as.data.table
#' @noRd
setMethod(f="calcMutSpectra",
          signature="MutSpectra",
          definition=function(object, verbose, ...){
              # access the part of the object we want to manipulate
              primaryData <- object@primaryData
              
              # print status message
              if(verbose){
                  memo <- paste("calculating rates of transitions and transversions")
                  message(memo)
              }
              
              # calculate the frequency of transitions/transversions on a sample basis
              frequencyTransTranv <-  table(primaryData$trans_tranv, primaryData$sample)
              
              # calculate the proportion of transitions/transversions on a sample basis
              proportionTransTranv <-  prop.table(frequencyTransTranv, 2)
              
              # combine the frequency and proportion of transition/transversion
              frequencyTransTranv <- data.table::as.data.table(frequencyTransTranv)
              colnames(frequencyTransTranv) <- c("TransTranv", "Sample", "Frequency")
              
              proportionTransTranv <- data.table::as.data.table(proportionTransTranv)
              colnames(proportionTransTranv) <- c("TransTranv", "Sample", "Proportion") 
              
              primaryData <- merge(frequencyTransTranv, proportionTransTranv,
                                   by=c("TransTranv", "Sample"))
              
              # make sure columns are of appropriate type
              primaryData$Sample <- factor(primaryData$Sample)
              
              return(primaryData)
          })

#' @rdname MutSpectra-methods
#' @aliases sortSamples,MutSpectra
#' @param object Object of class MutSpectra
#' @param verbose Boolean for status updates
#' @return data.table object with samples refactored for plotting
#' @importFrom gtools mixedsort
#' @noRd
setMethod(f="sortSamples",
          signature="MutSpectra",
          definition=function(object, sorting, verbose, ...){
              # access the part of the object we want to manipulate
              primaryData <- object@primaryData
              
              # do nothing if samples is null
              if(is.null(sorting)){
                  return(primaryData)
              }
              
              # make sure the sorting parameter is a character
              if(!is.character(sorting)){
                  memo <- paste("the input to parameter \"sorting\" is not of type",
                                "character, attempting to coerce!")
                  warning(memo)
                  sorting <- as.character(sorting)
              }
              
              # Perform sorting based on the input to sorting
              if(toupper(sorting) == "SAMPLE" | toupper(sorting) == "SAMPLES"){
                  if(verbose){
                      memo <- paste("sorting samples by name")
                      message(memo)
                  }
                  sampleOrder <- gtools::mixedsort(primaryData$Sample)
                  primaryData$Sample <- factor(primaryData$Sample, levels=sampleOrder)
              }else if(toupper(sorting) == "MUTATION" | toupper(sorting) == "MUTATIONS"){
                  if(verbose){
                      memo <- paste("sorting samples by transition/transversion proportions")
                      message(memo)
                  }
                  sampleOrder <- primaryData[order(primaryData$TransTranv, -primaryData$Proportion)]
                  sampleOrder <- sampleOrder[sampleOrder$Proportion != 0,]
                  sampleOrder <- unique(sampleOrder$Sample)
                  primaryData$Sample <- factor(primaryData$Sample, levels=sampleOrder)
              }else if(length(sorting) > 1){
                  if(verbose){
                      memo <- paste("sorting samples by supplied samples in the parameter \"sorting\"")
                      message(memo)
                  }
                  # check if input to sorting is missing samples
                  if(any(!primaryData$Sample %in% sorting)){
                      expecSamples <- unique(primaryData$Sample)
                      missingSamples <-  expecSamples[!expecSamples %in% sorting]
                      memo <- paste("The following samples were missing from the parameter", 
                                    "\"sorting\": ", toString(missingSamples),
                                    "adding these to the end of the sort order")
                      warning(memo)
                      sorting <- c(sorting, missingSamples)
                  }
                  
                  # check if input to sorting has extra samples
                  if(any(!sorting %in% primaryData$Sample)){
                      extraSamples <- sorting[!sorting %in% primaryData$Sample]
                      memo <- paste("The following samples were specified in the parameter",
                                    "\"sorting\" but were not found in the primary data:",
                                    toString(extraSamples), "removing these samples!")
                      sorting <- sorting[!sorting %in% extraSamples]
                  }
                  primaryData$Sample <- factor(primaryData$Sample, levels=sorting)
              }
              
              return(primaryData)
          })

#' @rdname buildClinicalPlot-methods
#' @aliases buildClinicalPlot,MutSpectra
#' @param clinicalLayers List of ggplot2 layers to add to the clinical plot.
#' @noRd
#' @import ggplot2
#' @importFrom gtable gtable
setMethod(f="buildClinicalPlot",
          signature="MutSpectra",
          definition=function(object, clinicalLayers, ...){
              # extract necessary data
              clinicalData <- object@ClinicalData
              
              # if clinical data is empty return empty gtable
              if(nrow(clinicalData) == 0){
                  return(gtable::gtable())
              }
              
              # construct a plot
              clinicalXLabel <- xlab(paste0("Sample n=", length(unique(clinicalData$sample))))
              clinicalPlot <- ggplot(clinicalData, aes_string(x='sample',
                                                              y='variable',
                                                              fill='value')) +
                  clinicalLayers + clinicalXLabel + theme(axis.text.x=element_blank())
              
              # contruct grob
              clinicalGrob <- ggplotGrob(clinicalPlot)
          })

#' @rdname buildFrequencyPlot-methods
#' @aliases buildFrequencyPlot,MutSpectra
#' @param palette Character vector specifying the colors used for encoding transitions and transversions
#' , should be of length 6. If NULL a default palette will be used.
#' @param verbose Boolean specifying if status messages should be reported
#' @param plotALayers list of ggplot2 layers to be passed to the frequency plot.
#' @return gtable object containing the plot of frequencies for transitions and transversions
#' @noRd
#' @import ggplot2
#' @importFrom gtable gtable
setMethod(f="buildFrequencyPlot",
          signature="MutSpectra",
          definition=function(object, plotALayers, palette, verbose, ...){
              # extract the data to plot
              primaryData <- object@primaryData
              
              # print status message
              if(verbose){
                  memo <- paste("Constructing a plot for the observed",
                                "frequencies of transitions/transversions")
                  message(memo)
              }

              # construct bar geom
              barGeom <- geom_bar(data=primaryData, mapping=aes_string(x='Sample', y='Frequency', fill='TransTranv'), stat='identity', width=1)
              
              # add the palette for encoding transitions/transversions
              if(is.null(palette)){
                  palette <- c("#77C55D", "#A461B4", "#C1524B", "#93B5BB", "#4F433F", "#BFA753")
              } else if(length(palette) != 6){
                  memo <- paste("The input to the parameter \"palette\" is not of length",
                                "6, a color should be specified for each transition/transversion in:",
                                toString(unique(primaryData$TransTranv), "! using the default palette."))
                  warning(memo)
                  palette <- c("#77C55D", "#A461B4", "#C1524B", "#93B5BB", "#4F433F", "#BFA753")
              }
              plotPalette <- scale_fill_manual("Transitions/Transversions", values=palette)
              
              # add labels
              plotLabels <- labs(y="Frequency")
              
              # set theme 
              plotTheme <- theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x=element_blank(),
                                 panel.background=element_blank())
              
              # adjust limits
              yLimits <- scale_y_continuous(expand=c(0,0))
              
              # perform quality checks
              if(!is.null(plotALayers) && !is.list(plotALayers)){
                  memo <- paste("plotALayers is not a list... attempting to coerce.")
                  warning(memo)
                  plotALayers <- as.list(plotALayers)
                  if(any(lapply(plotALayers, function(x) ggplot2::is.ggproto(x) | ggplot2::is.theme(x)))){
                      memo <- paste("plotALayers is not a list of ggproto or ",
                                    "theme objects... setting plotALayers to NULL")
                      warning(memo)
                      plotALayers <- NULL
                  }
              }
              
              # initalize the plot
              frequencyPlot <- ggplot() + barGeom + plotLabels + plotTheme + yLimits +
                  plotPalette + plotALayers
              
              # contruct grob
              frequencyGrob <- ggplotGrob(frequencyPlot)
              
              return(frequencyGrob)
          })

#' @rdname buildProportionPlot-methods
#' @aliases buildProportionPlot,MutSpectra
#' @param palette Character vector specifying the colors used for encoding transitions and transversions
#' , should be of length 6. If NULL a default palette will be used.
#' @param sampleNames Boolean specifying if samples should be labeled on the plot.
#' @param verbose Boolean specifying if status messages should be reported
#' @param plotBLayers list of ggplot2 layers to be passed to the proportion plot.
#' @return gtable object containing the plot of proportions for transitions and transversions
#' @noRd
#' @import ggplot2
setMethod(f="buildProportionPlot",
          signature="MutSpectra",
          definition=function(object, sampleNames, plotBLayers, palette, verbose, ...){
              # extract the data to plot
              primaryData <- object@primaryData
              
              # print status message
              if(verbose){
                  memo <- paste("Constructing a plot for the observed",
                                "proportions of transitions/transversions")
                  message(memo)
              }
              
              # construct bar geom
              barGeom <- geom_bar(data=primaryData, mapping=aes_string(x='Sample', y='Proportion', fill='TransTranv'), stat='identity', width=1)
              
              # add the palette for encoding transitions/transversions
              if(is.null(palette)){
                  palette <- c("#77C55D", "#A461B4", "#C1524B", "#93B5BB", "#4F433F", "#BFA753")
              } else if(length(palette) != 6){
                  memo <- paste("The input to the parameter \"palette\" is not of length",
                                "6, a color should be specified for each transition/transversion in:",
                                toString(unique(primaryData$TransTranv), "! using the default palette."))
                  warning(memo)
                  palette <- c("#77C55D", "#A461B4", "#C1524B", "#93B5BB", "#4F433F", "#BFA753")
              }
              plotPalette <- scale_fill_manual("Transitions/Transversions", values=palette)
              
              # add labels
              xLabel <- paste0("Sample: n=", length(unique(primaryData$Sample)))
              plotLabels <- labs(x=xLabel, y="Proportion")
              
              # add theme
              plotTheme <- theme(axis.text.x=element_text(angle=45, hjust=1),
                                 panel.background=element_blank())
              if(!sampleNames){
                  plotTheme <- plotTheme + theme(axis.text.x=element_blank())
              }
              
              # adjust limits
              yLimits <- scale_y_continuous(expand=c(0,0))
              
              # add additional plot layers
              if(!is.null(plotBLayers) && !is.list(plotBLayers)){
                  memo <- paste("plotBLayers is not a list... attempting to coerce.")
                  warning(memo)
                  plotBLayers <- as.list(plotBLayers)
                  if(any(lapply(plotBLayers, function(x) ggplot2::is.ggproto(x) | ggplot2::is.theme(x)))){
                      memo <- paste("plotBLayers is not a list of ggproto or ",
                                    "theme objects... setting plotBLayers to NULL")
                      warning(memo)
                      plotBLayers <- NULL
                  }
              }
              
              # initalize the plot
              proportionPlot <- ggplot() + barGeom + plotLabels + plotTheme + yLimits +
                  plotPalette + plotBLayers
              
              # contruct grob
              proportionGrob <- ggplotGrob(proportionPlot)
              
              return(proportionGrob)
          })

#' @rdname arrangeMutSpectraPlot-methods
#' @aliases arrangeMutSpectraPlot,MutSpectra
#' @param sectionHeights Relative heights of each plot section (should sum to one).
#' @noRd
#' @importFrom grid nullGrob
#' @importFrom gridExtra arrangeGrob
setMethod(f="arrangeMutSpectraPlot",
          signature="MutSpectra",
          definition=function(object, sectionHeights, verbose, ...){
              
              # grab the data we need
              plotA <- object@PlotA
              plotB <- object@PlotB
              plotC <- object@PlotC
              
              # Obtain the max width for relevant plots
              plotList <- list(plotA, plotB, plotC)
              plotList <- plotList[lapply(plotList, length) > 0]
              plotWidths <- lapply(plotList, function(x) x$widths)
              maxWidth <- do.call(grid::unit.pmax, plotWidths)
              
              # Set the widths for all plots
              for(i in 1:length(plotList)){
                  plotList[[i]]$widths <- maxWidth
              }
              
              # set section heights based upon the number of sections
              defaultPlotHeights <- c(.5, .5, .25)
              
              if(is.null(sectionHeights)){
                  if(length(plotList) < 3){
                      defaultPlotHeights <- defaultPlotHeights[-length(defaultPlotHeights)]
                  }
                  sectionHeights <- defaultPlotHeights
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
              finalPlot <- do.call(gridExtra::arrangeGrob, c(plotList, list(ncol=1, heights=sectionHeights)))
              
              return(finalPlot)
          })

#' @rdname formatClinicalData-methods
#' @aliases formatClinicalData,MutSpectra
#' @noRd
#' @importFrom data.table setDT
#' @importFrom data.table is.data.table
#' @importFrom data.table melt
#' @importFrom data.table data.table
#' @importFrom data.table rbindlist
setMethod(f="formatClinicalData",
          signature="MutSpectra",
          definition=function(object, verbose, ...){
              
              # extract the data we need
              clinicalData <- object@ClinicalData
              primaryData <- object@primaryData
              
              # print status message
              if(verbose){
                  memo <- paste("Formatting clinical data")
                  message(memo)
              }
              
              # remove clinical samples not found in the primary data
              primaryDataSamples <- levels(primaryData$Sample)
              clinicalData <- clinicalData[clinicalData$sample %in% primaryDataSamples,]
              removedSamp <- unique(object@ClinicalData$sample[!object@ClinicalData$sample %in% clinicalData$sample])
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
              } else if(verbose) {
                  message(memo)
              }
              
              # set levels of clinicalData to match primaryData for ordering
              clinicalData$sample <- factor(clinicalData$sample, levels=levels(primaryData$Sample))
              clinicalData$value <- factor(clinicalData$value, levels=unique(clinicalData$value))
              
              # return the formated data
              return(clinicalData)
          })