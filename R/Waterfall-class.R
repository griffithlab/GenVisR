#' Class Waterfall
#' 
#' An S4 class for the waterfall plot object
#' @name Waterfall-class
#' @rdname Waterfall-class
#' @slot plotA gtable object for the top sub-plot.
#' @slot plotB gtable object for the left sub-plot.
#' @slot plotC gtable object for the main plot.
#' @slot plotD gtable object for the bottom sub-plot.
#' @slot Grob gtable object for the arranged plot.
#' @slot Data list of data.table objects storing the plotted data.
#' @exportClass Waterfall
#' @import methods
#' @importFrom gtable gtable

methods::setOldClass("gtable")
setClass("Waterfall",
         representation=representation(plotA="gtable",
                                       plotB="gtable",
                                       plotC="gtable",
                                       plotD="gtable",
                                       Grob="gtable",
                                       Data="list")
         )

#' Initalizer method for the Waterfall class
#' 
#' @name Waterfall
#' @rdname Waterfall-class
setMethod(f="initialize",
          signature="Waterfall",
          definition=function(.Object, input){
              
              # convert to waterfall format
              primaryData <- "test"
          })

#' Constructor for the Waterfall class.
#' 
#' @name Waterfall
#' @rdname Waterfall-class
#' @param input MutationAnnotationFormat class holding genomic information.
#' @export
Waterfall <- function(input){
    cat("!!!!! Waterfall~Constructor !!!!!\n")
    new("Waterfall", input=input)
}

