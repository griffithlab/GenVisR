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
#' @import ggplot2
setMethod(f="initialize",
          signature="MutSpectra",
          definition=function(.Object, input, BSgenome, clinical, sectionHeights,
                              sectionWidths, verbose, plotCLayers){
              # convert object to mutSpectra format
              toMutSpectra(input, BSgenome=BSgenome, verbose=verbose)
          })

#' Constructor for the MutSpectra class.
#' 
#' @name MutSpectra
#' @rdname MutSpectra-class
#' @param input Object of class MutationAnnotationFormat, GMS, VEP.
#' @param BSgenome Object of class BSgenome, used to extract reference bases if
#' not supplied by the file format.
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
MutSpectra <- function(input, BSgenome=NULL, clinical=NULL, sectionHeights=NULL,
                      sectionWidths=NULL, verbose=FALSE, plotCLayers=NULL){
    cat("!!!!! MutSpectra~Constructor !!!!!\n")
    new("MutSpectra", input=input, BSgenome=BSgenome, clinical=clinical, sectionHeights=sectionHeights,
        sectionWidths=sectionWidths, verbose=verbose, plotCLayers=plotCLayers)
}