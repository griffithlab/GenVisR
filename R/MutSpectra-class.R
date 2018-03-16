################################################################################
##################### Public/Private Class Definitions #########################

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Public Class !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

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
         }
)

#' Constructor for the MutSpectra class.
#' 
#' @name MutSpectra
#' @rdname MutSpectra-class
#' @param object Object of class MutationAnnotationFormat, GMS, VEP.
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
#' @param sampleNames Boolean specifying if samples should be labeled on the plot.
#' @param verbose Boolean specifying if status messages should be reported
#' @param plotALayers list of ggplot2 layers to be passed to the frequency plot.
#' @param plotBLayers list of ggplot2 layers to be passed to the proportion plot.
#' @param plotCLayers list of ggplot2 layers to be passed to the clinical plot.
#' @export
MutSpectra <- function(object, BSgenome=NULL, sorting=NULL, palette=NULL, clinical=NULL, sectionHeights=NULL,
                       sampleNames=TRUE, verbose=FALSE, plotALayers=NULL, plotBLayers=NULL, plotCLayers=NULL){
    
    # initalize the MutSpectraPrimaryData object
    primaryData <- MutSpectraPrimaryData(object, BSgenome=BSgenome, sorting=sorting, verbose=verbose)
    
    # initalize the MutSpectraPlots object
    plots <- MutSpectraPlots(primaryData=primaryData, clinical=clinical,
                             plotALayers=plotALayers, plotBLayers=plotBLayers,
                             plotCLayers=plotCLayers, sampleNames=sampleNames,
                             palette=palette, verbose=verbose)
    
    # get the clinical data
    if(is.null(clinical)){
        ClinicalData <- data.table::data.table()
    } else {
        ClinicalData <- getData(clinical)
    }
    
    # align all the plots together
    Grob <- arrangeMutSpectraPlot(plots, sectionHeights, verbose)
    
    new("MutSpectra", PlotA=getGrob(plots, 1), PlotB=getGrob(plots, 2), PlotC=getGrob(plots, 3),
        Grob=Grob, primaryData=getData(primaryData, name="primaryData"), ClinicalData=ClinicalData)
}

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Private Classes !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

#' Private Class MutSpectraPrimaryData
#' 
#' An S4 class for the Primary Data of the MutSpectra plot object
#' @name MutSpectraPrimaryData-class
#' @rdname MutSpectraPrimaryData-class
#' @slot primaryData data.table object storing the primary data, should have
#' column names sample, mutation, frequency, proportion.
#' @import methods
#' @importFrom data.table data.table
#' @noRd
setClass("MutSpectraPrimaryData",
         representation=representation(primaryData="data.table"),
         validity=function(object){
         }
)

#' Constructor for the MutSpectraPrimaryData class.
#' 
#' @name MutSpectraPrimaryData
#' @rdname MutSpectraPrimaryData-class
#' @param object Object of class MutationAnnotationFormat, GMS, VEP
#' @noRd
MutSpectraPrimaryData <- function(object, BSgenome, sorting, verbose){
    # convert object to mutSpectra format
    primaryData <- toMutSpectra(object, BSgenome=BSgenome, verbose=verbose)
    
    # annotate types of transitions and transversions
    primaryData <- annoMutSpectra(primaryData, verbose=verbose)
    
    # calculate rates of transitions and transversions
    primaryData <- calcMutSpectra(primaryData, verbose=verbose)
    
    # sort the samples for plotting
    primaryData <- sortSamples(primaryData, sorting=sorting, verbose=verbose)
    
    # initalize the object
    new("MutSpectraPrimaryData", primaryData=primaryData)
}

#' Private Class MutSpectraPlots
#' 
#' An S4 class for the plots of the MutSpectra class
#' @name MutSpectraPlots-class
#' @rdname MutSpectraPlots-class
#' @slot PlotA gtable object for the mutation frequencies.
#' @slot PlotB gtable object for the mutation proportions.
#' @slot PlotC gtable object for clinical data sub-plot.
#' @import methods
#' @importFrom gtable gtable
#' @noRd
methods::setOldClass("gtable")
setClass("MutSpectraPlots",
         representation=representation(PlotA="gtable",
                                       PlotB="gtable",
                                       PlotC="gtable"),
         validity=function(object){
         }
)

#' Constructor for the MutSpectraPlots class.
#' 
#' @name MutSpectraPlots
#' @rdname MutSpectraPlots-class
#' @param object Object of class MutationAnnotationFormat
#' @importFrom gtable gtable
#' @noRd
MutSpectraPlots <- function(primaryData, clinical, plotALayers, plotBLayers, plotCLayers, sampleNames, palette, verbose){
    
    if(is(clinical, "Clinical")){
        # format the clinical data object using the mutSpectra object
        ClinicalData <- formatClinicalData(primaryData, clinical, verbose)
        
        # make sure the new formating conforms with ClinicalData class 
        ClinicalData <- ClinicalData(ClinicalData, inputFormat="long", verbose=verbose)
        
        # add layer to suppress the x-axis text in the proportion plot
        addLayer <- theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
        plotBLayers[[length(plotBLayers ) + 1]] <- addLayer
        
        # add the clinical data plot
        PlotC <- buildClinicalPlot(ClinicalData, clinicalLayers=getLayers(clinical), verbose=verbose)
    } else {
        PlotC <- gtable::gtable()
    }
    
    # construct the frequency plot
    PlotA <- buildFrequencyPlot(primaryData, plotALayers, palette, verbose)
    
    # construct the proportion plot
    PlotB <- buildProportionPlot(primaryData, sampleNames, plotBLayers, palette, verbose)
    
    # initalize the object
    new("MutSpectraPlots", PlotA=PlotA, PlotB=PlotB, PlotC=PlotC)
}

################################################################################
###################### Accessor function definitions ###########################

#' Helper function to getData from classes
#'
#' @rdname getData-methods
#' @aliases getData
.getData_MutSpectra <- function(object, name=NULL, index=NULL, ...){
    
    if(is.null(name) & is.null(index)){
        memo <- paste("Both name and index are NULL, one must be specified!")
        stop(memo)
    }
    
    if(is.null(index)){
        index <- 0
    } else {
        if(index > 2){
            memo <- paste("index out of bounds")
            stop(memo)
        }
    }
    
    if(is.null(name)){
        name <- "noMatch"
    } else {
        slotAvailableName <- c("primaryData", "ClinicalData")
        if(!(name %in% slotAvailableName)){
            memo <- paste("slot name not found, specify one of:", toString(slotAvailableName))
            stop(memo)
        }
    }
    
    if(name == "primaryData" | index == 1){
        data <- object@primaryData
    } else if(name == "ClinicalData" | index == 2){
        data <- object@ClinicalData
    }
    
    return(data)
}

#' @rdname getData-methods
#' @aliases getData
setMethod(f="getData",
          signature="MutSpectraPrimaryData",
          definition=.getData_MutSpectra)

#' @rdname getData-methods
#' @aliases getData
setMethod(f="getData",
          signature="MutSpectra",
          definition=.getData_MutSpectra)

#' @rdname drawPlot-methods
#' @aliases drawPlot
#' @importFrom grid grid.draw
#' @importFrom grid grid.newpage
#' @exportMethod drawPlot
setMethod(f="drawPlot",
          signature="MutSpectra",
          definition=function(object, ...){
              grob <- getGrob(object, index=4)
              grid::grid.newpage()
              grid::grid.draw(grob)
          })

#' Helper function to extract grobs from objects
#'
#' @rdname getGrob-methods
#' @aliases getGrob
#' @noRd
.getGrob_MutSpectra <- function(object, index=1, ...){
    if(index == 1){
        grob <- object@PlotA
    } else if(index == 2) {
        grob <- object@PlotB
    } else if(index == 3) {
        grob <- object@PlotC
    } else if(index == 4) {
        grob <- object@Grob
    } else {
        stop("Subscript out of bounds")
    }
    return(grob)
}

#' @rdname getGrob-methods
#' @aliases getGrob
setMethod(f="getGrob",
          signature="MutSpectraPlots",
          definition=.getGrob_MutSpectra)

#' @rdname getGrob-methods
#' @aliases getGrob
setMethod(f="getGrob",
          signature="MutSpectra",
          definition=.getGrob_MutSpectra)

################################################################################
####################### Method function definitions ############################

#' @rdname toMutSpectra-methods
#' @aliases toMutSpectra
#' @param object Object of class data.table
#' @param verbose Boolean specifying if status messages should be reported
#' @importFrom data.table rbindlist
#' @noRd
setMethod(f="toMutSpectra",
          signature="data.table",
          definition=function(object, verbose, ...){
              
              # print status message
              if(verbose){
                  memo <- paste("Converting", class(object)[1],
                                "to expected MutSpectra format")
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
              mutSectraFormat <- object
              
              # grab only the snvs
              snvIndex <- which(nchar(mutSectraFormat$refAllele) == 1 & nchar(mutSectraFormat$variantAllele) == 1)
              if(verbose){
                  memo <- paste("Removing", nrow(mutSectraFormat)-length(snvIndex),
                                "entries which are not snvs!")
                  message(memo)
              }
              
              # make sure snvs are of the expected base type
              expectedBases <- c("C", "G", "T", "A")
              if(!all(toupper(mutSectraFormat$refAllele) %in% expectedBases)){
                  abnormalBasesIndex_ins <-  !mutSectraFormat$refAllele %in% expectedBases
              }
              if(!all(toupper(mutSectraFormat$variantAllele) %in% expectedBases)){
                  abnormalBasesIndex_del <-  !mutSectraFormat$variantAllele %in% expectedBases
              }
              abnormalBasesIndex <- unique(c(abnormalBasesIndex_ins, abnormalBasesIndex_del))
              if(length(abnormalBasesIndex > 1)){
                  memo <- paste("Removing", length(abnormalBasesIndex), "entries with unexpected",
                                "bases in either the refAllele or variantAllele columns!")
                  mutSectraFormat <- mutSectraFormat[-abnormalBasesIndex,]
              }
              
              # unique, to make sure no duplicate variants exist to throw off the counts
              rowCountOrig <- nrow(mutSectraFormat)
              mutSectraFormat <- unique(mutSectraFormat)
              
              # print status message
              if(verbose){
                  memo <- paste("Removed", rowCountOrig - nrow(mutSectraFormat),
                                "rows from the data which harbored duplicate",
                                "genomic locations")
                  message(memo)
              }
              
              # Remove cases where there is not change between reference and variant
              mutSpectraFormat$refAllele <- as.character(mutSpectraFormat$refAllele)
              mutSpectraFormat$variantAllele <- as.character(mutSpectraFormat$variantAllele)
              alleleMatchIndex <- mutSpectraFormat$refAllele == mutSpectraFormat$variantAllele
              mutSpectraFormat <- mutSpectraFormat[!alleleMatchIndex,]
              if(verbose){
                  memo <- paste("Removed", length(alleleMatchIndex), "entries",
                                "where the reference allele matched the tumor allele")
                  message(memo)
              }
              
              # convert appropriate columns to factor
              mutSpectraFormat$sample <- factor(mutSpectraFormat$sample)
              
              return(mutSpectraFormat)
          })

#' @rdname toMutSpectra-methods
#' @aliases toMutSpectra
#' @param object Object of class data.table
#' @param verbose Boolean specifying if status messages should be reported
#' @importFrom data.table as.data.table
#' @noRd
setMethod(f="toMutSpectra",
          signature="data.frame",
          definition=function(object, verbose, ...){
              
              # print status message
              if(verbose){
                  memo <- paste("Converting", class(object),
                                "to data.table")
                  message(memo)
              }
              
              # convert object to data.table
              object <- as.data.table(object)
              
              # then convert to MutSpectra format
              object <- toMutSpectra(object, verbose=verbose)
              
              return(object)
          })

#' @rdname MutSpectra-methods
#' @aliases MutSpectra
#' @param object Object of class data.table
#' @param verbose Boolean for status updates
#' @return data.table object with transitions and transversions annotated
#' @noRd
setMethod(f="annoMutSpectra",
          signature="data.table",
          definition=function(object, verbose, ...){
              # access the part of the object we want to manipulate
              primaryData <- object
              
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
#' @aliases MutSpectra
#' @param object Object of class data.table
#' @param verbose Boolean for status updates
#' @return data.table object with transitions and transversions rates calculated
#' @importFrom data.table as.data.table
#' @noRd
setMethod(f="calcMutSpectra",
          signature="data.table",
          definition=function(object, verbose, ...){
              # access the part of the object we want to manipulate
              primaryData <- object
              
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
#' @aliases MutSpectra
#' @param object Object of class data.table
#' @param verbose Boolean for status updates
#' @return data.table object with samples refactored for plotting
#' @importFrom gtools mixedsort
#' @noRd
setMethod(f="sortSamples",
          signature="data.table",
          definition=function(object, sorting, verbose, ...){
              # access the part of the object we want to manipulate
              primaryData <- object
              
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
              
              if(length(sorting) == 1){
                  
                  # Perform sorting based on the input to sorting
                  if(toupper(sorting) == "SAMPLE" | toupper(sorting) == "SAMPLES"){
                      if(verbose){
                          memo <- paste("sorting samples by name")
                          message(memo)
                      }
                      sampleOrder <- unique(gtools::mixedsort(primaryData$Sample))
                  }else if(toupper(sorting) == "MUTATION" | toupper(sorting) == "MUTATIONS"){
                      if(verbose){
                          memo <- paste("sorting samples by transition/transversion proportions")
                          message(memo)
                      }
                      sampleOrder <- primaryData[order(primaryData$TransTranv, -primaryData$Proportion)]
                      sampleOrder <- sampleOrder[sampleOrder$Proportion != 0,]
                      sampleOrder <- unique(sampleOrder$Sample)
                  }
              } else {
                  
                  if(verbose){
                      memo <- paste("sorting samples by supplied samples in the parameter \"sorting\"")
                      message(memo)
                  }
                  
                  # check if input to sorting is missing samples
                  sampleOrder <- sorting
                  if(any(!primaryData$Sample %in% sampleOrder)){
                      expecSamples <- unique(primaryData$Sample)
                      missingSamples <-  as.character(expecSamples[!expecSamples %in% sampleOrder])
                      memo <- paste("The following samples were missing from the parameter", 
                                    "\"sorting\": ", toString(missingSamples),
                                    "adding these to the end of the sort order")
                      warning(memo)
                      sampleOrder <- c(sampleOrder, missingSamples)
                  }
                  
                  # check if input to sorting has extra samples
                  if(any(!sampleOrder %in% primaryData$Sample)){
                      extraSamples <- sorting[!sorting %in% primaryData$Sample]
                      memo <- paste("The following samples were specified in the parameter",
                                    "\"sorting\" but were not found in the primary data:",
                                    toString(extraSamples), "removing these samples!")
                      sampleOrder <- sampleOrder[!sampleOrder %in% extraSamples]
                  }
              }
             
              primaryData$Sample <- factor(primaryData$Sample, levels=sampleOrder)
                  
              return(primaryData)
          })

#' @rdname buildFrequencyPlot-methods
#' @aliases MutSpectra
#' @param object Object of class MutSpectraPrimaryData
#' @param palette Character vector specifying the colors used for encoding transitions and transversions
#' , should be of length 6. If NULL a default palette will be used.
#' @param verbose Boolean specifying if status messages should be reported
#' @param plotALayers list of ggplot2 layers to be passed to the frequency plot.
#' @return gtable object containing the plot of frequencies for transitions and transversions
#' @noRd
#' @import ggplot2
#' @importFrom gtable gtable
setMethod(f="buildFrequencyPlot",
          signature="MutSpectraPrimaryData",
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
                                toString(unique(primaryData$TransTranv)), "! using the default palette.")
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
              
              # initalize the plot
              frequencyPlot <- ggplot() + barGeom + plotLabels + plotTheme + yLimits +
                  plotPalette + plotALayers
              
              # contruct grob
              frequencyGrob <- ggplotGrob(frequencyPlot)
              
              return(frequencyGrob)
          })

#' @rdname buildProportionPlot-methods
#' @aliases MutSpectra
#' @param object Object of class MutSpectraPrimaryData
#' @param palette Character vector specifying the colors used for encoding transitions and transversions
#' , should be of length 6. If NULL a default palette will be used.
#' @param sampleNames Boolean specifying if samples should be labeled on the plot.
#' @param verbose Boolean specifying if status messages should be reported
#' @param plotBLayers list of ggplot2 layers to be passed to the proportion plot.
#' @return gtable object containing the plot of proportions for transitions and transversions
#' @noRd
#' @import ggplot2
setMethod(f="buildProportionPlot",
          signature="MutSpectraPrimaryData",
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
                                toString(unique(primaryData$TransTranv)), "! using the default palette.")
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
              
              # initalize the plot
              proportionPlot <- ggplot() + barGeom + plotLabels + plotTheme + yLimits +
                  plotPalette + plotBLayers
              
              # contruct grob
              proportionGrob <- ggplotGrob(proportionPlot)
              
              return(proportionGrob)
          })

#' @rdname arrangeMutSpectraPlot-methods
#' @aliases arrangeMutSpectraPlot
#' @param sectionHeights Relative heights of each plot section (should sum to one).
#' @noRd
#' @importFrom grid nullGrob
#' @importFrom gridExtra arrangeGrob
setMethod(f="arrangeMutSpectraPlot",
          signature="MutSpectraPlots",
          definition=function(object, sectionHeights, verbose, ...){
              
              # grab the data we need
              plotA <- getGrob(object, 1)
              plotB <- getGrob(object, 2)
              plotC <- getGrob(object, 3)
              
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
              defaultPlotHeights <- c(1, 1, .75)
              
              if(is.null(sectionHeights)){
                  if(length(plotList) < 3){
                      defaultPlotHeights <- defaultPlotHeights[-length(defaultPlotHeights)]
                  }
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

#' @rdname formatClinicalData-methods
#' @aliases formatClinicalData,MutSpectraPrimaryData
#' @noRd
#' @importFrom data.table setDT
#' @importFrom data.table is.data.table
#' @importFrom data.table melt
#' @importFrom data.table data.table
#' @importFrom data.table rbindlist
setMethod(f="formatClinicalData",
          signature="MutSpectraPrimaryData",
          definition=function(object, clinicalData, verbose, ...){

              # extract the data we need
              primaryData <- getData(object, name="primaryData")
              clinicalData <- getData(clinicalData)
              
              # print status message
              if(verbose){
                  memo <- paste("Formatting clinical data")
                  message(memo)
              }
              
              # remove clinical samples not found in the primary data
              primaryDataSamples <- levels(primaryData$Sample)
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
              } else if(verbose) {
                  message(memo)
              }
              
              # set levels of clinicalData to match primaryData for ordering
              clinicalData$sample <- factor(clinicalData$sample, levels=levels(primaryData$Sample))
              clinicalData$value <- factor(clinicalData$value, levels=unique(clinicalData$value))
              
              # return the formated data
              return(clinicalData)
          })