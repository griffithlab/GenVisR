#' Class Clinical
#' 
#' An S4 class to store clinical information and plots
#' @name Clinical-class
#' @rdname Clinical-class
#' @slot clinicalGrob gtable object for the clinical plot.
#' @slot clinicalLayers list of ggtheme or ggproto objects used to build the plot.
#' @slot clinicalData data.table object to store clinical data
#' @exportClass Clinical
#' @import methods
#' @importFrom gtable gtable
#' @importFrom data.table data.table

methods::setOldClass("gtable")
setClass("Clinical",
         representation=representation(clinicalGrob="gtable",
                                       clinicalLayers="list",
                                       clinicalData="data.table"),
         validity=function(object){
             cat("!!!!! Clinical~Inspector !!!!!\n")
         }
)

#' Initalizer method for the Clinical class
#' 
#' @name Clinical
#' @rdname Clinical-class
#' @noRd
#' @importFrom data.table fread
#' @importFrom data.table is.data.table

setMethod(f="initialize",
          signature="Clinical",
          definition=function(.Object, path, inputData, inputFormat, legendColumns,
                              palette, clinicalLayers, verbose){
              cat("!!!!!!!!!!!!!! Clinical~Initalizer !!!!!!!!!!!!!!!!!!")
              # Obtain the clinical data
              if(!is.null(inputData)){
                  if(!is.data.table(inputData)){
                      data.table::setDT(inputData)
                  }
                  .Object@clinicalData <- inputData
              } else {
                  .Object@clinicalData <- data.table::fread(input=path,
                                                            stringsAsFactors = TRUE,
                                                            verbose=verbose)
              }
              
              # format clinical data
              .Object@clinicalData <- formatClinicalData(.Object, inputFormat=inputFormat, verbose=verbose)
              
              # construct a list of ggplot layers
              .Object@clinicalLayers <- setClinicalPlotLayers(.Object, legendColumns,
                                                              palette, clinicalLayers,
                                                              verbose)
              
              # construct a plot of clinical data
              .Object@clinicalGrob <- buildClinicalPlot(.Object, verbose)
              
              return(.Object)
          })

#' Constructor for the Clinical class.
#' 
#' @name Clinical
#' @rdname Clinical-class
#' @param path String specifying the path to clinical data, file must have the
#' column "sample".
#' @param inputData Optional data.table or data.frame object holding clinical
#' data, used only if path is not specified. Data must have the column "sample".
#' @param inputFormat String specifying the input format of the data given, one
#' of wide or long format (see details).
#' @param legendColumns Integer specifying the number of columns in the legend.
#' @param palette Named character vector supplying colors for clinical variables.
#' @param clinicalLayers list of ggplot2 layers to be passed to the plot.
#' @param verbose Boolean specifying if progress should be reported.
#' @details The Clinical() function is a constructor to create a GenVisR object
#' of class Clinical. This is used to both display clinical data in the form
#' of a heatmap and to add clinical data to various GenVisR plots.
#' Input to this function can be either the path to a file containing clinical
#' information using the parameter "path", or alternatively a data.table object
#' if this information into R. By default the input is assumed to be in a wide
#' format where each variable has it's own column, in such cases the data will
#' be coerced into a long format where there is a key->value pair mapping to
#' the data. The assumption of "wide"/"long" format can be changed with the 
#' "inputFormat" parameter, in both cases there should be a column called
#' "sample" within the data supplied which is used as an id variable.
#' @seealso \code{\link{getData}}
#' @seealso \code{\link{drawPlot}}
#' @export
Clinical <- function(path, inputData=NULL ,inputFormat=c("wide", "long"),
                     legendColumns=1, palette=NULL, clinicalLayers=NULL,
                     verbose=FALSE){
    cat("!!!!! Clinical~Constructor !!!!!\n")
    new("Clinical", path=path, inputData=inputData, inputFormat=inputFormat,
        legendColumns=legendColumns, palette=palette,
        clinicalLayers=clinicalLayers, verbose=verbose)
}

#' @rdname formatClinicalData-methods
#' @aliases formatClinicalData,MutationAnnotationFormat
#' @noRd
#' @importFrom data.table setDT
#' @importFrom data.table is.data.table
#' @importFrom data.table melt
setMethod(f="formatClinicalData",
          signature="Clinical",
          definition=function(object, inputFormat, verbose, ...){
              
              # print status message
              if(verbose){
                  memo <- paste("Formatting clinical data")
                  message(memo)
              }
              
              # extract the data we need
              clinicalData <- unique(object@clinicalData)
              inputFormat <- inputFormat[1]
              
              # quality checks
              if(!toupper(inputFormat) %in% toupper(c("wide", "long"))){
                  memo <- paste("inputFormat is:", inputFormat, "... should be",
                                "one of \"wide\", \"long\"!")
                  stop(memo)
              }
              browser()
              if(ncol(clinicalData) > 3 && toupper(inputFormat) == toupper("long")){
                  memo <- ("Number of columns found to be > 3... assuming wide format")
                  warning(memo)
                  inputFormat="wide"
              }
              if(toupper(inputFormat)==toupper("long") && all(!toupper(colnames(clinicalData)) %in% toupper(c("sample", "variable", "value")))){
                  missingColumn <- colnames(clinicalData)[!toupper(colnames(clinicalData)) %in% toupper(c("sample", "variable", "value"))]
                  memo <- paste("the columns:",missingColumn,"is missing!")
                  stop(memo)
              }
              if(!data.table::is.data.table(clinicalData)){
                  memo <- paste("input is not a data.table... attempting to coerce!")
                  data.table::setDT(clinicalData)
              }
              if(all(!toupper(colnames(clinicalData)) %in% toupper("sample"))){
                  memo <- paste("Could not find the column name \"sample\"")
                  stop(memo)
              }
              
              # if data is in wide format coerce to long
              colnames(clinicalData) <- tolower(colnames(clinicalData))
              if(toupper(inputFormat)==toupper("wide")){
                  clinicalData <- data.table::melt(clinicalData, id.vars="sample")
              }
              
              # coerce to factor for ordering in plot call
              clinicalData$value <- factor(clinicalData$value, levels=unique(clinicalData$value))
              
              # return the formated data
              return(clinicalData)
          })

#' @rdname setClinicalPlotLayers-methods
#' @aliases setClinicalPlotLayers,MutationAnnotationFormat
#' @param legendColumns Integer specifying the number of columns in the legend.
#' @param palette Named character vector supplying colors for clinical variables.
#' @param clinicalLayers list of ggplot2 layers to be passed to the plot.
#' @param verbose Boolean specifying if status messages should be printed.
#' @noRd
#' @import ggplot2
setMethod(f="setClinicalPlotLayers",
          signature="Clinical",
          definition=function(object, legendColumns, palette, clinicalLayers, verbose, ...){
              # print status message
              if(verbose){
                  memo <- paste("building clinical plot layers")
                  message(memo)
              }
              
              # quality checks
              if(legendColumns %% 1 != 0){
                  memo <- paste("legendColumns is not a whole number, attempting to coerce.")
                  legendColumns <- as.integer(legendColumns)
              }
              if(!is.null(clinicalLayers) && !is.list(clinicalLayers)){
                  memo <- paste("clinicalLayers is not a list... attempting to coerce.")
                  warning(memo)
                  clinicalLayers <- as.list(clinicalLayers)
                  if(any(lapply(clinicalLayers, function(x) ggplot2::is.ggproto(x) | ggplot2::is.theme(x)))){
                      memo <- paste("clinicalLayers is not a list of ggproto or ",
                                    "theme objects... setting clinicalLayers to NULL")
                      warning(memo)
                      clinicalLayers <- NULL
                  }
              }
              
              ################## define layers ###################
              # legend
              plotLegendGuide <- guides(fill=guide_legend(ncol=legendColumns, title="Clinical"))
              if(!is.null(palette)){
                  plotLegend <- scale_fill_manual("Clinical", values=palette,
                                                  breaks=names(palette))
              } else {
                  plotLegend <- geom_blank()
              }
              
              # theme
              plotTheme <- theme(panel.grid.major = element_blank(),
                                 panel.grid.minor=element_blank(),
                                 panel.background=element_rect(fill='white',
                                                               colour='white'),
                                 axis.ticks.x=element_blank(),
                                 axis.ticks.y=element_blank(),
                                 axis.title.x=element_text(size=16),
                                 legend.title=element_text(size=14),
                                 axis.title.y=element_blank(),
                                 axis.text.y=element_text(size=14, colour='black'),
                                 axis.text.x=element_text(angle=90),
                                 legend.position='right')
              
              # geom
              plotGeom <- geom_tile()
              
              # return list of layers
              return(list(plotGeom, plotLegendGuide, plotTheme, plotLegend, clinicalLayers))
          })

#' @rdname buildClinicalPlot-methods
#' @aliases buildClinicalPlot,MutationAnnotationFormat
#' @noRd
#' @import ggplot2
setMethod(f="buildClinicalPlot",
          signature="Clinical",
          definition=function(object, ...){
              # extract necessary data
              clinicalData <- object@clinicalData
              clinicalLayers <- object@clinicalLayers
              
              # construct a plot
              clinicalXLabel <- xlab(paste0("Sample n=", length(unique(clinicalData$sample))))
              clinicalPlot <- ggplot(clinicalData, aes_string(x='sample',
                                                              y='variable',
                                                              fill='value')) +
                  clinicalLayers + clinicalXLabel
              
              # contruct grob
              clinicalGrob <- ggplotGrob(clinicalPlot)
          })

#' @rdname getData-methods
#' @aliases getData,Clinical
setMethod(f="getData",
          signature="Clinical",
          definition=function(object, ...){
              clinData <- object@clinicalData
              return(clinData)
          })

#' @rdname drawPlot-methods
#' @aliases drawPlot,Clinical
#' @importFrom grid grid.draw
setMethod(
    f="drawPlot",
    signature="Clinical",
    definition=function(object, ...){
        mainPlot <- object@clinicalGrob
        grid::grid.draw(mainPlot)
    }
)