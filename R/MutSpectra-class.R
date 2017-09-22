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
          definition=function(.Object, input, BSgenome, sorting, clinical, sectionHeights,
                              sectionWidths, verbose, plotALayers, plotBLayers,
                              plotCLayers){
              # convert object to mutSpectra format
              .Object@primaryData <- toMutSpectra(input, BSgenome=BSgenome, verbose=verbose)
              
              # annotate types of transitions and transversions
              .Object@primaryData <- annoMutSpectra(.Object, verbose)
              
              # calculate rates of transitions and transversions
              .Object@primaryData <- calcMutSpectra(.Object, verbose)
              
              # sort the samples for plotting
              .Object@primaryData <- sortSamples(.Object, sorting ,verbose)
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
#' @param clinical Object of class Clinical, used for adding a clinical data subplot.
#' @param sectionHeights Numeric vector specifying relative heights of each plot section,
#' should sum to one. Expects a value for each section.
#' @param sectionWidths Numeric vector specifying relative heights of each plot section,
#' should sum to one. Expects a value for each section.
#' @param verbose Boolean specifying if status messages should be reported
#' @param plotALayers list of ggplot2 layers to be passed to the frequency plot.
#' @param plotBLayers list of ggplot2 layers to be passed to the proportion plot.
#' @param plotCLayers list of ggplot2 layers to be passed to the clinical plot.
#' @export
MutSpectra <- function(input, BSgenome=NULL, sorting=NULL, clinical=NULL, sectionHeights=NULL,
                      sectionWidths=NULL, verbose=FALSE, plotALayers=NULL, plotBLayers=NULL,
                      plotCLayers=NULL){
    cat("!!!!! MutSpectra~Constructor !!!!!\n")
    new("MutSpectra", input=input, BSgenome=BSgenome, sorting=sorting, clinical=clinical, sectionHeights=sectionHeights,
        sectionWidths=sectionWidths, verbose=verbose, plotALayers=plotALayers, plotBLayers = plotBLayers, 
        plotCLayers=plotCLayers)
}

#' @rdname MutSpectra-methods
#' @aliases annoMutSpectra,Waterfall
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
#' @aliases calcMutSpectra,Waterfall
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
              
              return(primaryData)
          })

#' @rdname MutSpectra-methods
#' @aliases calcMutSpectra,Waterfall
#' @param object Object of class MutSpectra
#' @param verbose Boolean for status updates
#' @return data.table object with samples refactored for plotting
#' @importFrom data.table as.data.table
#' @noRd
setMethod(f="calcMutSpectra",
          signature="MutSpectra",
          definition=function(object, sorting, verbose, ...){
              # access the part of the object we want to manipulate
              primaryData <- object@primaryData
              
              # print status message
              if(verbose){
                  memo <- paste("attempting to sort samples")
                  message(memo)
              }
              
              if(toupper(sorting) == "SAMPLE" | toupper(sorting) == "SAMPLES"){
                  #START HERE
                  #sample_order <- as.vector(unique(x$sample))
                  #sample_order <- gtools::mixedsort(sample_order)
                  #x$sample <- factor(x$sample, levels=sample_order)
              }else if(toupper(sorting) == "mutation" | toupper(sorting) == "MUTATIONS"){
                  
              }

          })