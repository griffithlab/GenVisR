################################################################################
##################### Public/Private Class Definitions #########################

#' Private Class cnLohPlots
#' 
#' An S4 class for the cn, somatic loh, and germline loh plots
#' @rdname cnLohPlots-class
#' @name cnLohPlots
#' @slot cnPlot gtable object for the cn plot
#' @slot somaticLohPlot gtable object for the somatic loh plot
#' @slot germlineLohPlot gtable object for the germline loh plot
#' @import methods
#' @importFrom gtable gtable
#' @noRd
setClass("cnLohPlots",
         representation=representation(cnLohPlot="list"),
         validity=function(object) {
             
         })

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Public Class !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#' Class cnLoh
#' 
#' An S4 class for the cn, somatic loh, and germline loh plots
#' @name cnLoh-class
#' @rdname cnLoh-class
#' @slot cnData data.table object for cn plot
#' @slot somaticLohData data.table object for the somatic loh plot
#' @slot germlineLohData data.table object for the germline loh plot
#' @slot cnLohPlots gtable object for the combined plots
#' @exportClass cnLoh
#' @import methods
#' @importFrom data.table data.table
#' @importFrom gtable gtable
methods::setOldClass("gtable")
setClass(Class="cnLoh",
         representation=representation(cnData="data.table",
                                       somaticLohData="data.table",
                                       germlineLohData="data.table",
                                       cnLohPlots="cnLohPlots"),
         validity=function(object){
        
    }
)

#' Constructor for the cnLoh class
#' 
#' @name cnLoh
#' @rdname cnLoh-class
#' @param cnInput Object of class cnLohDataFormat
#' @param somaticLohInput
#' @param germlineLohInput
#' @param samples Character vector specifying samples to plot. If not NULL
#' all samples in "input" not specified with this parameter are removed.
#' @param chromosomes Character vector specifying chromosomes to plot. If not NULL
#' all chromosomes in "input" not specified with this parameter are removed.
#' @param cnvType
#' @param BSgenome Object of class BSgenome to extract genome wide chromosome 
#' coordinates
#' @param getHeterozygousLohCalls
#' @param plotAColor
#' @param plotALayers
#' @param somaticLohCutoff
#' @param plotBAlpha
#' @param plotBColors
#' @param plotBLayers
#' @param plotCLimits
#' @param plotCColors
#' @param plotCLayers
#' @param sectionHeights
#' @param verbose
#' @export
cnLoh <- function(cnInput, somaticLohInput, germlineLohInput, samples, chromosomes, 
                  cnvType, BSgenome, 
                  getHeterozygousLohCalls, plotAColor, plotALayers, plotBAlpha,
                  somaticLohCutoff, plotBColors, plotBLayers,
                  plotCLimits, plotCColors, plotCLayers, 
                  sectionHeights, verbose) {
    
    ## Check each of the input parameters
    cnLohInputParameters <- checkCombinedCnLohInputParameters(cnInput=cnInput, 
                                                              somaticLohInput=somaticLohInput, 
                                                              germlineLohInput=germlineLohInput, 
                                                              samples=samples, chromosomes=chromosomes, 
                                                              cnvType=cnvType, BSgenome=BSgenome,
                                                              getHeterozygousLohCalls=getHeterozygousLohCalls,
                                                              somaticLohCutoff=somaticLohCutoff,
                                                              plotAColor=plotAColor, plotALayers=plotALayers,
                                                              plotBAlpha=plotBAlpha, plotBColors=plotBColors,
                                                              plotBLayers=plotBLayers, plotClimits=plotCLimits,
                                                              plotCColors=plotCColors, plotCLayers=plotCLayers,
                                                              sectionHeights=sectionHeights, verbose=verbose)
    
    ## Obtain cn, somatic loh, and germline loh datasets to plot
    cnLohDataset <- cnLohData(cnInput=cnInput, somaticLohInput=somaticLohInput,
                              germlineLohInput=germlineLohInput,
                              samples=cnLohInputParameters@samples, 
                              chromosomes=cnLohInputParameters@chromosomes,
                              cnvType=cnLohInputParameters@cnvType, 
                              BSgenome=cnLohInputParameters@BSgenome,
                              getHeterozygousLohCalls=cnLohInputParameters@getHeterozygousLohCalls, 
                              verbose=cnLohInputParameters@verbose)
    
    ## Generate the cn, somatic LOH, and germline LOH plots
    plots <- cnLohPlots(object=cnLohDataset,
                        cnvType=cnLohInputParameters@cnvType,
                        somaticLohCutoff=cnLohInputParameters@somaticLohCutoff, 
                        plotAColor=cnLohInputParameters@plotAColor, 
                        plotALayers=cnLohInputParameters@plotALayers, 
                        plotBAlpha=cnLohInputParameters@plotBAlpha, 
                        plotBColors=cnLohInputParameters@plotBColors,
                        plotBLayers=cnLohInputParameters@plotBLayers, 
                        plotCLimits=cnLohInputParameters@plotCLimits, 
                        plotCColors=cnLohInputParameters@plotCColors, 
                        plotCLayers=cnLohInputParameters@plotCLayers, 
                        sectionHeights=cnLohInputParameters@sectionHeights,
                        verbose=cnLohInputParameters@verbose)
    
    ## Initialize the object
    new("cnLoh", cnData=getData(cnLohDataset, index=1), 
        somaticLohData=getData(cnLohDataset, index=2), 
        germlineLohData=getData(cnLohDataset, index=3), 
        cnLohPlots=plots)
}

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Private Classes !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

#' Private Class cnLohInputParameters
#' 
#' An S4 class to check input parameters of the combinedCnLoh function
#' @name cnLohInputParameters-class
#' @noRd
setClass("cnLohInputParameters",
         representation=representation(samples="character",
                                       chromosomes="character", 
                                       cnvType="character", BSgenome="BSgenome",
                                       getHeterozygousLohCalls="logical",
                                       somaticLohCutoff="numeric",
                                       plotAColor="character", plotALayers="list",
                                       plotBAlpha="numeric", plotBColors="character",
                                       plotBLayers="list", plotCLimits="numeric",
                                       plotCColors="character", plotCLayers="list",
                                       sectionHeights="numeric", verbose="logical"), 
         validity=function(object){
             
         })

#' Constructor for the cnLohInputParameters class
#' 
#' @name cnLohInputParameters
#' @rdname cnLohInputParameters-class
#' @noRd
checkCombinedCnLohInputParameters <- function(cnInput, somaticLohInput, germlineLohInput, samples, chromosomes, 
                                              cnvType, BSgenome, getHeterozygousLohCalls,
                                              somaticLohCutoff, plotAColor, plotALayers,
                                              plotBAlpha, plotBColors, plotBLayers, plotClimits,
                                              plotCColors, plotCLayers, sectionHeights, verbose) {
    ##### Check the verbose parameter #####
    ## Check to see if verbose is a boolean
    if (!is.logical(verbose) | is.null(verbose)) {
        memo <- paste0("The verbose parameter is not a boolean (T/F). Coercing verbose to be FALSE...")
        message(memo)
        verbose <- FALSE
    }
    
    ##### Check the samples parameter ##### 
    ## Check if the samples in the somaticLohInput, germlineLohInput, and cnInput all match
    cnSamples <- as.character(cnInput@sample$sample)
    somaticLohSamples <- as.character(somaticLohInput@sample$sample)
    germlineLohSamples <- as.character(germlineLohInput@sample$sample)
    allSamples <- unique(c(cnSamples, somaticLohSamples, germlineLohSamples))
    
    ## Check is samples is NULL
    if (is.null(samples)) {
        samples <- allSamples
        samples <- factor(samples, levels=gtools::mixedsort(samples))
        memo <- paste0("Sample parameter cannot be NULL. All samples will be plotted.")
        message(memo)
    }
    ## Check if samples is a character vector
    if (!is.character(samples)) {
        memo <- paste0("samples variable not of the character class. Attempting to coerce.")
        samples <- as.character(samples)
        message(memo)
    }
    
    ## Check if the designated samples is in the sv dataset
    if (!is.null(samples)) {
        `%nin%` = Negate(`%in%`)
        discrepantSamples <- paste(samples[which(samples %nin% allSamples)], collapse=", ")
        if (length(discrepantSamples) > 0 & discrepantSamples != "") {
            memo <- paste0("The desired samples: ", discrepantSamples, " are not found ",
                           "in either of the cnv, somatic loh, or germline loh dataset. Available sample names include: ", 
                           paste(allSamples, collapse=", "), ". ",
                           "Please designate valid sample names.")
            stop(memo)
        }
    }
    
    ##### Check the chromosomes parameter #####
    if (is.null(chromosomes)) {
        chromosomes <- paste("chr", seq(1:22), sep="")
        memo <- paste0("chromosomes parameter cannot be NULL. Using all autosomes...")
        message(memo)
    }
    ## Check if chromosomes is a character vector
    if (!is.character(chromosomes)) {
        memo <- paste0("chromosomes variable not of the character class. Attempting to coerce.")
        chromosomes <- as.character(chromosomes)
        message(memo)
    }
    ## Check if it has the "chr" prefix
    # Check to see if the chromosomes variable has "chr" in front if not NULL, autosomes, or all
    if (all(chromosomes!="autosomes") & all(chromosomes!="all")) {
        if (!all(grepl("^chr", chromosomes))) {
            if (verbose) {
                memo <- paste0("Did not detect the prefix chr in the chromosomes specified ",
                               "in the `chromosomes` variable... adding prefix")
                message(memo)
                chromosomes <- paste("chr", chromosomes, sep="")
            } 
        } else if (all(grepl("^chr", chromosomes))) {
            if (verbose) {
                memo <- paste0("Detected chr in the `chromosomes` variable...",
                               "proceeding")
                message(memo)
            } 
        } else {
            memo <- paste0("Detected unknown or mixed prefixes in the `chromosomes`` variable",
                           " colum of object... should either be chr or non (i.e.) chr1 or 1")
            message(memo)
        }
    }
    
    ##### Check the cnvType parameter #####
    ## Check if cnvType is not null
    if (is.null(cnvType)){
        memo <- paste0("the cnvType parameter cannot be null. Using absolute copy number values.")
        cnvType <- "absolute"
        message(memo)
    }
    ## Check if cnvType is a character
    if (!is.character(cnvType)) {
        memo <- paste0("cnvType paramter not of class character. Attempting to coerce...")
        message(memo)
        cnvType <- as.character(cnvType)
    }
    ## Check if cnvType has a valid value
    if (!(cnvType %in% c("logratio", "absolute", "relative"))) {
        memo <- paste0("cnvType parameter is not a valid value. Valid values include: ",
                       "logration, absolute, and relative. Using absolute copy number values.")
    }
    
    ##### Check the BSgenome parameter #####
    ## Check to see if BSgenome is a BSgenome
    if (is.null(BSgenome)) {
        memo <- paste("BSgenome object is not specified. This parameter is required ",
                      "to get the lengths of the chromsomes being plotted.")
        stop(memo)
    } else if (is(BSgenome, "BSgenome")) {
        memo <- paste("BSgenome passed object validity checks")
        message(memo)
    } else {
        memo <- paste("class of the BSgenome object is", class(BSgenome),
                      ". Should be of class BSgenome. This parameter is required ",
                      "to get the lengths of the chromsomes being plotted.")
        stop(memo)
        BSgenome <- NULL
    }
    
    ##### Check the getHeterozygousLohCalls #####
    if (!is.logical(getHeterozygousLohCalls) | is.null(getHeterozygousLohCalls)) {
        memo <- paste0("The getHeterozygousLohCalls parameter is not a boolean (T/F). ",
                       "Coercing getHeterozygousLohCalls to be TRUE...")
        message(memo)
        getHeterozygousLohCalls <- FALSE
    }
    
    ##### Check the somaticLohCutoff #####
    if (is.null(somaticLohCutoff)) {
        somaticLohCutoff <- 0.1
        memo <- paste0("somaticLohCutoff parameter cannot be NULL. Setting the somaticLohCutoff value to 0.1.")
        message(memo)
    }
    ## Check if somaticLohCutoff is numeric
    if (!is.numeric(somaticLohCutoff)) {
        memo <- paste0("somaticLohCutoff variable not of the numeric class. Attempting to coerce.")
        somaticLohCutoff <- as.numeric(somaticLohCutoff)
        message(memo)
    }
    
    ##### Check the plotALayers, plotBLayers, and plotCLayers #####
    checkPlotLayer <- function(plotLayer, name) {
        if(!is.null(plotLayer)){
            if(!is.list(plotLayer)){
                memo <- paste(name, " is not a list", sep="")
                stop(memo)
            }
            
            if(any(!unlist(lapply(plotLayer, function(x) ggplot2::is.ggproto(x) | ggplot2::is.theme(x) | is(x, "labels"))))){
                memo <- paste(name, " is not a list of ggproto or ",
                              "theme objects... setting plotALayers to NULL", sep="")
                warning(memo)
                plotLayer <- NULL
            }
        }
        return(plotLayer)
    }
    plotALayers <- checkPlotLayer(plotLayer=plotALayers, "plotALayers")
    plotBLayers <- checkPlotLayer(plotLayer=plotBLayers, "plotBLayers")
    plotCLayers <- checkPlotLayer(plotLayer=plotCLayers, "plotCLayers")
    
    ##### Check the plotAColor #####
    ## Check if it is a character vector
    if (is.null(plotAColor)) {
        memo <- paste0("plotAColor was set to NULL. Using default color.")
        message(memo)
        plotAColor <- c("Blue")
    }
    if (!is.character(plotAColor)) {
        memo <- paste0("plotAColor variable not of the character class. Attempting to coerce.")
        plotAColor <- as.character(plotAColor)
        message(memo)
    }
    areColors <- function(x) {
        sapply(x, function(X) {
            tryCatch(is.matrix(col2rgb(X)), 
                     error = function(e) FALSE)
        })
    }
    if (any(areColors(plotAColor) == FALSE)) {
        ## Get the invalid color
        nonColor <- plotAColor[which(data.table(areColors(plotAColor))$V1==FALSE)]
        memo <- paste0("The ", nonColor, " designated in the plotAColor parameter is not a valid color. ",
                       "Making the plotAColor to be: Blue.")
        
        message(memo)
        plotAColor <- c("Blue")
    }
    
    ##### Check the plotBAlpha #####
    if (is.null(plotBAlpha)) {
        plotBAlpha <- 0.75
        memo <- paste0("plotBAlpha parameter cannot be NULL. Setting the plotBAlpha value to 0.75")
        message(memo)
    }
    ## Check if plotBAlpha is numeric
    if (!is.numeric(plotBAlpha)) {
        memo <- paste0("plotBAlpha variable not of the numeric class. Attempting to coerce.")
        plotBAlpha <- as.numeric(plotBAlpha)
        message(memo)
    }
    
    ##### Check the plotBColors #####
    ## Check if it is a character vector
    if (is.null(plotBColors)) {
        memo <- paste0("plotBColors was set to NULL. Using default colors. ",
                       "Tumor vaf calls will be ",
                       "dark green and normal vaf calls will be dark red.")
        message(memo)
        plotBColors <- c("Dark Green", "Dark Red")
    }
    if (!is.character(plotBColors)) {
        memo <- paste0("plotBColors variable not of the character class. Attempting to coerce.")
        plotBColors <- as.character(plotBColors)
        message(memo)
    }
    areColors <- function(x) {
        sapply(x, function(X) {
            tryCatch(is.matrix(col2rgb(X)), 
                     error = function(e) FALSE)
        })
    }
    if (any(areColors(plotBColors) == FALSE)) {
        ## Get the invalid color
        nonColor <- plotBColors[which(data.table(areColors(plotBColors))$V1==FALSE)]
        memo <- paste0("The ", nonColor, " designated in the plotBColors parameter is not a valid color. ",
                       "Making the plotBColors to be: dark green and dark red. Tumor vaf calls will be ",
                       "dark green and normal vaf calls will be dark red.")
        
        message(memo)
        plotBColors <- c("Dark Green", "Dark Red")
    }
    
    ##### Check the plotCLimits #####
    if (is.null(plotCLimits)) {
        plotCLimits <- 25
        memo <- paste0("plotCLimits parameter cannot be NULL. Setting the plotCLimits value to 25")
        message(memo)
    }
    ## Check if plotCLimits is numeric
    if (!is.numeric(plotCLimits)) {
        memo <- paste0("plotCLimits variable not of the numeric class. Attempting to coerce.")
        plotCLimits <- as.numeric(plotCLimits)
        message(memo)
    }
    
    ##### Check the plotCColors #####
    ## Check if it is a character vector
    if (is.null(plotCColors)) {
        memo <- paste0("plotCColors was set to NULL. Using default colors. ",
                       "Tumor vaf calls will be ",
                       "dark green and normal vaf calls will be dark red.")
        message(memo)
        plotCColors <- c("white", "Dark Red")
    }
    if (!is.character(plotCColors)) {
        memo <- paste0("plotCColors variable not of the character class. Attempting to coerce.")
        plotCColors <- as.character(plotCColors)
        message(memo)
    }
    areColors <- function(x) {
        sapply(x, function(X) {
            tryCatch(is.matrix(col2rgb(X)), 
                     error = function(e) FALSE)
        })
    }
    if (any(areColors(plotCColors) == FALSE)) {
        ## Get the invalid color
        nonColor <- plotCColors[which(data.table(areColors(plotCColors))$V1==FALSE)]
        memo <- paste0("The ", nonColor, " designated in the plotCColors parameter is not a valid color. ",
                       "Making the plotCColors to be: white and dark red.")
        
        message(memo)
        plotCColors <- c("white", "Dark Red")
    }
    ##### Check the sectionHeights #####
    ## Check if it not NULL
    if (is.null(sectionHeights)) {
        sectionHeights <- c(0.3, 0.4, 0.3)
        memo <- paste0("sectionHeights variable cannot be NULL. Using default values.")
        message(memo)
    }
    
    ## Check that values are numeric
    if (!is.numeric(sectionHeights)) {
        memo <- paste0("sectionHeights valures are not class numeric. Attempting to coerce...")
        message(memo)
        sectionHeights <- as.numeric(sectionHeights)
    }
    
    ## Check that the values are > 0
    if (any(sectionHeights<0)) {
        memo <- paste0("sectionHeights cannot be a negative value. Using default values.")
        message(memo)
        sectionHeights <- c(0.3, 0.4, 0.3)
    }
    
    ## Check that there are 2 values in the variable
    if (length(sectionHeights)!=3) {
        memo <- paste0("3 values must be supplied to the sectionHeights parameter, which specifies the ",
                       "relative height of the cn plot, somatic loh plot, and germline plot respecively.")
        message(memo)
        sectionHeights <- c(0.3, 0.4, 0.3)
    }
    
    ## Check that the values sum up to 1
    if (sum(sectionHeights)!=1) {
        memo <- paste0("sectionHeight values do not equal 1. Using default values.")
        message(memo)
        sectionHeights <- c(0.3, 0.4, 0.3)
    }
    
    new("cnLohInputParameters", samples=samples, chromosomes=chromosomes, BSgenome=BSgenome,
        cnvType=cnvType, getHeterozygousLohCalls=getHeterozygousLohCalls,
        somaticLohCutoff=somaticLohCutoff, plotAColor=plotAColor, plotALayers=plotALayers,
        plotBAlpha=plotBAlpha, plotBColors=plotBColors, plotBLayers=plotBLayers, 
        plotCLimits=plotCLimits, plotCColors=plotCColors, plotCLayers=plotCLayers,
        sectionHeights=sectionHeights, verbose=verbose)    
}

#' Private Class cnLohData
#' 
#' An S4 class for the data to plot cn, somatic LOH, and germline LOH plots
#' @name cnLohData
#' @rdname cnLohData-class
setClass("cnLohData",
         representation=representation(rawCnData="data.table",
                                       segCnData="data.table",
                                       rawLohData="data.table",
                                       segLohData="data.table",
                                       rawGermlineLohData="data.table",
                                       chrData="data.table"),
         validity = function(object){
             
         })

#' Constructor for the cnLohData class
#' 
#' @name cnLohData
#' @rdname cnLohData-class
#' @param object Object of class cnLohDataFormat
cnLohData <- function(cnInput, somaticLohInput, germlineLohInput, samples, 
                      getHeterozygousLohCalls, chromosomes, cnvType, BSgenome, verbose) {

    ############################################################################
    #################### Prepare copy number variant dataset ###################
    ## Obtain raw cnv data
    cnData <- getCnvData(object=cnInput, cnvType=cnvType, verbose=verbose)
    
    ## Subset copy number data by chromosome
    cnData <- chrSubset(object=cnData, chromosomes=chromosomes, verbose=verbose)
    
    ## Subset copy number data by sample
    cnData <- sampleSubset(object=cnData, samples=samples, verbose=verbose)
    
    ## Obtain chromosome boundaries from BSgenome object
    chrData <- annoGenomeCoord(object=cnData, BSgenome=BSgenome, verbose=verbose)
    
    ## Obtain copy number segmentation data
    cnSegmentation <- getCnSegmentation(object=cnData, verbose=verbose)
    
    ## Remove gaps
    cnSegmentation <- removeGapsSegmentation(object=cnSegmentation, chrData=chrData, verbose=verbose)
    
    ############################################################################
    ##################### Prepare somatic loh dataset ##########################
    ## Obtain LOH data for desired chromosomes and samples
    lohData <- getLohData(object=somaticLohInput, verbose=verbose, getHeterozygousLohCalls=TRUE)
    
    ## Subset loh data by chromosome
    lohData <- chrSubset(object=lohData, chromosomes=chromosomes, verbose=verbose)
    
    ## Subset loh data by sample
    lohData <- sampleSubset(object=lohData, samples=samples, verbose=verbose)
    
    ## Obtain loh segmentation dataset
    lohSegmentation <- getLohSegmentation(object=lohData, verbose=verbose)
    
    ## Remove gaps 
    lohSegmentation <- removeGapsSegmentation(object=lohSegmentation, chrData=chrData,
                                              verbose=verbose)
    
    ############################################################################
    ##################### Prepare germline loh dataset #########################
    ## Obtain germlineloh data by chromosome
    germlineLohData <- getLohData(object=germlineLohInput, verbose=verbose, getHeterozygousLohCalls=FALSE)
    
    ## Subset loh data by chromosome
    germlineLohData <- chrSubset(object=germlineLohData, chromosomes=chromosomes, verbose=verbose)
    
    ## Subset loh data by sample
    germlineLohData <- sampleSubset(object=germlineLohData, samples=samples, verbose=verbose)
    
    ## Initialize the object
    new("cnLohData", rawCnData=cnData, segCnData=cnSegmentation, 
        rawLohData=lohData, segLohData=lohSegmentation, rawGermlineLohData=germlineLohData,
        chrData=chrData)
} 

#' Constructor for cnLohPlots class
#' 
#' @rdname cnLohPlots-class
#' @name cnLohPlots
#' @param object Object of class data.table
#' @importFrom gtable gtable
#' @import ggplot2
#' @noRd 
cnLohPlots <- function(object, cnvType, plotAColor, plotALayers, 
                       somaticLohCutoff, plotBAlpha, plotBColors, plotBLayers, 
                       plotCLimits, plotCColors, plotCLayers, sectionHeights, verbose) {
    
    ## Build the cn plot
    cnLohPlot <- buildCnLohPlot(object=object, cnvType=cnvType, plotAColor=plotAColor, plotALayers=plotALayers, 
                             plotBAlpha=plotBAlpha, plotBColors=plotBColors, plotBLayers=plotBLayers,
                             somaticLohCutoff=somaticLohCutoff, plotCLimits=plotCLimits,
                             plotCColors=plotCColors, plotCLayers=plotCLayers, 
                             sectionHeights=sectionHeights, verbose=verbose)
    
    ## Initialize the object
    new("cnLohPlots", cnLohPlot=cnLohPlot)
}

################################################################################
###################### Accessor function definitions ###########################

#' Helper function to get data from classes
#' 
#' @rdname getData-methods
#' @aliases getData
.getData_combinedCnLoh <- function(object, name=NULL, index=NULL, ...) {
    if(is.null(name) & is.null(index)){
        memo <- paste("Both name and index are NULL, one must be specified!")
        stop(memo)
    }
    
    if(is.null(index)){
        index <- 0
    } else {
        if(index > 3){
            memo <- paste("index out of bounds")
            stop(memo)
        }
    }
    
    if(is.null(name)){
        name <- "noMatch"
    } else {
        slotAvailableName <- c("rawCnData", "rawLohData", "rawGermlineLohData")
        if(!(name %in% slotAvailableName)){
            memo <- paste("slot name not found, specify one of:", toString(slotAvailableName))
            stop(memo)
        }
    }
    
    if(name == "rawCnData" | index == 1){
        data <- object@rawCnData
    }
    if(name == "rawLohData" | index == 2){
        data <- object@rawLohData
    }
    if(name == "rawGermlineLohData" | index == 3){
        data <- object@rawGermlineLohData
    }
    
    return(data)
}

#' @rdname getData-methods
#' @aliases getData
setMethod(f="getData",
          signature="cnLohData",
          definition=.getData_combinedCnLoh)

#' @rdname getData-methods
#' @aliases getData
setMethod(f="getData",
          signature="cnLoh",
          definition=.getData_combinedCnLoh)

#' Helper function to extract grobs from objects
#'
#' @rdname getGrob-methods
#' @aliases getGrob
#' @noRd
.getGrob_combinedCnLoh <- function(object, index, ...){
    if(index == 1){
        grob <- object@cnPlot
    } else if(index == 2) {
        grob <- object@somaticLohPlot
    } else if (index == 3) {
        grob <- object@germlineLohPlot
    } else if (index == 4) {
        grob <- object@Grob
    } else {
        stop("Subscript out of bounds") 
    }
    return(grob)
}

#' @rdname getGrob-methods
#' @aliases getGrob
setMethod(f="getGrob",
          signature="cnLohPlots",
          definition=.getGrob_combinedCnLoh)

#' @rdname getGrob-methods
#' @aliases getGrob
setMethod(f="getGrob",
          signature="cnLoh",
          definition=.getGrob_combinedCnLoh)


#' @rdname drawPlot-methods
#' @aliases drawPlot
#' @importFrom grid grid.draw
#' @importFrom grid grid.newpage
#' @exportMethod drawPlot
setMethod(
    f="drawPlot",
    signature="cnLoh",
    definition=function(object, chr=NULL, sample=NULL, ...){
        ## Get the list of gtables 
        object <- object@cnLohPlots@cnLohPlot
        
        ## Get the chromosome-sample combinations
        name <- paste0(chr, "_", sample)
        
        ## See if the desired chr-sample combo can be found in the plots
        num <- which(names(object) == name)
        if (length(num) == 0) {
            memo <- paste0("The plot for the chromosome-sample combination: ",
                   name, " could not be found. Make sure to append the chr name ",
                   "with ", dQuote("chr"), " rather than just using the chromosome number (chr1 instead of 1). ",
                   "The possible combinations that could be used are: ", 
                   paste(names(object), collapse=", "))
            stop(memo)
        }
        
        mainPlot <- object[[num]]
        grid::grid.newpage()
        grid::grid.draw(mainPlot)
    }
)

################################################################################
###################### Method function definitions #############################

######################################################
##### Function to obtain chromosomes of interest #####
#' @rdname cnLoh-methods
#' @aliases cnLoh
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
              
              ## Check format of the chromosome column
              if (!all(grepl("^chr", object$chromosome))) {
                  memo <- paste0("Did not detect the prefix chr in the chromosome column ",
                                 "of x... adding prefix")
                  message (memo)
                 object$chromosome <- paste("chr", object$chromosome, sep="")
              } else if (all(grepl("^chr", object$chromosome))) {
                  memo <- paste0("Detected chr in the chromosome column of x...",
                                 "proceeding")
                  message(memo)
              } else {
                  memo <- paste0("Detected unknown or mixed prefixes in the chromosome",
                                 " colum of object... should either be chr or non (i.e.) chr1 or 1")
                  message(memo)
              }
              
              ## Determine which chromosomes to plot
              ## Only include autosomes
              if (chromosomes[1] == "autosomes") {
                  chromosomes <- as.character(c(seq(1:22)))
              }
              ## Include all chromosomes
              if (chromosomes[1] == "all") {
                  chromosomes <- unique(object$chromosome)
                  chromosomes <- chromosomes[-grep("GL", chromosomes)]
                  chromosomes <- chromosomes[-grep("MT", chromosomes)]
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

##################################################
##### Function to obtain samples of interest #####
#' @rdname sampleSubset-methods
#' @name sampleSubset
#' @aliases sampleSubset
#' @param object Object of class data.table
#' @param samples character vector of samples to retain
#' @param verbose Boolean for status updates
#' @return data.table object with calculated mutation distances
#' @noRd
setMethod(f="sampleSubset",
          signature="data.table",
          definition=function(object, samples, verbose, ...){
              
              # print status message
              if(verbose){
                  memo <- paste("Performing sample subsets")
                  message(memo)
              }
              
              ## Check for specified samples not in the original input
              missingSamp <- samples[!samples %in% unique(object$sample)]
              if (length(missingSamp) != 0) {
                  memo <- paste("The following samples were designated to be 
                                keptbut were not found:", toString(missingSamp), 
                                "\nValid csamples are", 
                                toString(unique(object$sample)))
                  warning(memo)
              } 
              
              ## Perform the subset
              object <- object[object$sample %in% samples]
              object$sample <- factor(object$sample)
              
              ## Check that the object has a size after subsets
              if(nrow(object) < 1){
                  memo <- paste("no entries left to plot after chromosome subsets")
                  stop(memo)
              }
              
              return(object)
          })

#############################################################
##### Function to generate segmentation dataset for cnv #####
#' @rdname getCnSegmentation-methods
#' @name getCnSegmentation
#' @aliases getCnSegmentation
#' @param object of class data.table
#' @noRd
setMethod(f="getCnSegmentation",
          signature="data.table",
          definition=function(object, verbose, ...) {
              
              ## Print status message
              if (verbose) {
                  message("Segmenting copy number data")
              }
              
              ## Get the sample-chr combination
              object$sample_chr_combo <- paste0(object$chromosome, "_", object$sample)
              
              ## Split object by sample
              segDfTemp <- split(object, f=object$sample_chr_combo)
              
              ## Perform segmentation
              segmentationDF <- rbindlist(lapply(segDfTemp, function(x) {
                  cnSeg <- CNA(genomdat=as.numeric(x$cn), chrom=x$chromosome,
                                maploc=x$position, data.type="logratio", sampleid = unique(x$sample_chr_combo))
                  
                  ## Run CBS
                  cnSeg <- segment(cnSeg, min.width=3, undo.splits="sdundo",
                                   undo.SD=2)
                  cnSeg <- cnSeg$output
                  return(cnSeg)
              }))
              return(segmentationDF)
          })

#############################################################
##### Function to generate segmentation dataset for cnv #####
#' @rdname removeGapsSegmentation-methods
#' @name removeGapsSegmentation
#' @aliases removeGapsSegmentation
#' @param object of class data.table
#' @noRd
setMethod(f="removeGapsSegmentation",
          signature="data.table",
          definition=function(object, chrData, verbose, ...) {
              
              ## Print status message
              if (verbose) {
                  message("Removing gaps from segmentation file")
              }
              
              ## Get the list of the chromosomes
              splitDf <- split(object, f=object$ID)
              segs <- rbindlist(lapply(splitDf, function(df) {
                  for (i in 1:(nrow(df) - 1)) {
                      ## Don't merge segments if they are far apart
                      if ((df$loc.start[i+1]-df$loc.end[i]) < 5000000) {
                          half <- floor((df$loc.end[i] + df$loc.start[i+1])/2)
                          df$loc.end[i] <- half
                          df$loc.start[i+1] <- half + 1
                      }
                  }
                  return(df)
              }))
              return(segs)
          })

#####################################################
##### Function to get the chromosome boundaries #####
#' @rdname annoGenomeCoord-methods
#' @name annoGenomeCoord
#' @aliases annoGenomeCoord
#' @param object Object of class data.table
#' @param BSgenome Object of class BSgenome, used for extracting chromosome boundaries
#' @param verbose Boolean for status updates
#' @return Data.table with chr and start/stop positions
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom data.table as.data.table
#' @importFrom data.table rbindlist
#' @importFrom gtools mixedsort
#' @noRd
setMethod(f="annoGenomeCoord", 
          signature="data.table",
          definition=function(object, BSgenome, verbose, ...){
              
              ## Print status message
              if (verbose) {
                  memo <- paste("Acquiring chromosome boundaries from BSgenome object")
                  message(memo)
              }
              
              ## Create a data table of genomic coordinates end positions
              genomeCoord <- data.table::as.data.table(seqlengths(BSgenome))
              colnames(genomeCoord) <- c("end")
              genomeCoord$chromosome <- names(seqlengths(BSgenome))
              genomeCoord$start <- 1
              
              ## Check that chromosomes between BSgenome and original input match
              chrMismatch <- as.character(unique(object[!object$chromosome %in% genomeCoord$chromosome,]$chromosome))
              if (length(chrMismatch) >= 1) {
                  memo <- paste("The following chromosomes do not match the supplied BSgenome object",
                                toString(chrMismatch))
                  warning(memo)
                  
                  ## Test if the chr mismatch is fixed by appending chr to chromosomes
                  chrMismatch_appendChr <- length(as.character(unique(object[!paste0("chr", object$chromosome) %in% genomeCoord$chromosome,]$chromosome)))
                  if(chrMismatch_appendChr < length(chrMismatch)){
                      memo <- paste("appending \"chr\" to chromosomes in attempt to fix mismatch with the BSgenome")
                      warning(memo)
                      object$chromosome <- paste0("chr", object$chromosome)
                  }
              }
              
              ## Check to see if any chromosomes in the original input dataset lack genomic coordiantes
              if (any(!unique(object$chromosome) %in% unique(genomeCoord$chromosome))) {
                  missingGenomeCoord <- unique(object$chromosome)
                  missingGenomeCoord <- missingGenomeCoord[!missingGenomeCoord %in% unique(genomeCoord_a$chromosome)]
                  memo <- paste("The following chromosomes are missing genomic coordinates", toString(missingGenomeCoord),
                                "Full genomic coordinates will not be plotted for these chromosomes")
                  warning(memo)
              }
              
              ## Filter the genomeCoord objext to only inlcude chromosomes in the input data
              genomeCoord <- genomeCoord[genomeCoord$chromosome %in% unique(object$chromosome),]
              
              return(genomeCoord)
              
          })

#############################################################
##### Function to generate segmentation dataset for loh #####
#' @rdname getLohSegmentation-methods
#' @name getLohSegmentation
#' @param object of class data.table
#' @param chrData of class data.table 
#' @aliases getLohSegmentation
setMethod(f = "getLohSegmentation", 
          signature="data.table",
          definition=function(object, ...){
              
              ## Print status message
              if (verbose) {
                  message("Determining segmeans from LOH calculations")
              }
              
              ## Get the sample-chr combination
              object$sample_chr_combo <- paste0(object$chromosome, "_", object$sample)
              
              ## Get the absolute loh difference
              object$absDiff <- abs(as.numeric(as.character(object$tumor_var_freq)) - 0.50)
              segDfTemp <- split(object, list(as.character(object$sample_chr_combo)))
              segmentationDf <- rbindlist(lapply(segDfTemp, function(x){
                  x$midpoint <- x$position
                  lohSeg <- CNA(genomdat = as.numeric(x$absDiff), chrom = x$chromosome,
                                maploc = x$midpoint, data.type = "binary", sampleid = unique(x$sample_chr_combo))
                  lohSeg <- segment(lohSeg)
                  lohSeg <- lohSeg$output
                  return(lohSeg)
              }))
              
              return(segmentationDf)
          })

########################################
##### Function to generate cn plot #####
#' @rdname buildCnLohPlot-methods
#' @name buildCnLohPlot
#' @aliases buildCnLohPlot
#' @param object of class data.table
#' @noRd
setMethod(f="buildCnLohPlot",
          signature="cnLohData",
          definition=function(object, cnvType, plotAColor, plotALayers,
                              plotBAlpha, plotBColors, plotBLayers, 
                              somaticLohCutoff, plotCLimits, plotCColors,
                              plotCLayers, sectionHeights, verbose){
              
              ## Print status message
              if (verbose) {
                  message("Building cnv plot")
              }
              object=cnLohDataset
              ## Get the cn data (segments and raw data)
              rawCnData <- object@rawCnData
              rawCnData$chr_sample_combo <- paste0(rawCnData$chromosome, "_", rawCnData$sample)
              segCnData <- object@segCnData
              
              ## Get the somatic loh data (segments and raw data)
              rawLohData <- object@rawLohData
              rawLohData$chr_sample_combo <- paste0(rawLohData$chromosome, "_", rawLohData$sample)
              segLohData <- object@segLohData
              segLohData <- segLohData[seg.mean > somaticLohCutoff]
              
              ## Get the germline loh data (segments and raw data)
              rawGermlineLohData <- object@rawGermlineLohData
              rawGermlineLohData$chr_sample_combo <- paste0(rawGermlineLohData$chromosome, "_", rawGermlineLohData$sample)
              
              ## Get the chromosome and sample data
              chrData <- object@chrData
              samples <- unique(c(as.character(unique(rawCnData$sample)), 
                                  as.character(unique(rawLohData$sample)), 
                                  as.character(unique(rawGermlineLohData$sample))))
              sample_chr_combo <- data.table(chr_sample_combo=as.vector(outer(chrData$chromosome, samples, paste, sep="_")))
              
              ## Split the dataset by chr_sample_combo
              splitDf <- split(sample_chr_combo, f=sample_chr_combo$chr_sample_combo)
              
              ## Generate combined cn/somatic loh/ germline loh plots 
              #x <- splitDf[[8]]$chr_sample_combo
              combinedPlots <- lapply(splitDf, function(x, chrData, cnvType, plotAColor, plotALayers,
                                                        plotBAlpha, plotBColors, plotBLayers, 
                                                        somaticLohCutoff, plotCLimits, plotCColors,
                                                        plotCLayers, sectionHeights, verbose){
                  ## Print status message
                  if (verbose) {
                      memo <- paste0("Generating combined cn/somatic loh/ germline loh plots")
                      message(memo)
                  }
                  
                  ## Subset out the datasets to chr and sample of interest
                  cnData <- rawCnData[chr_sample_combo==x]
                  cnSeg <- segCnData[ID==x]
                  somaticLohData <- rawLohData[chr_sample_combo==x]
                  lohSeg <- segLohData[ID==x]
                  germlineLohData  <- rawGermlineLohData[chr_sample_combo==x]
                  
                  ## Get the chromosome of interest
                  chr <- strsplit(as.character(x), split="_")[[1]][1]
                  chr <- chrData[chromosome==chr]
                  
                  ##############################################################              
                  ##### Build the copy number variation plot ###################
                  ##############################################################
                  ## Define parameters of the plot
                  plotTheme <- theme(axis.ticks.x=element_blank(),
                                     axis.text.x=element_blank(),
                                     axis.ticks.y=element_blank(),
                                     panel.grid.major=element_blank(),
                                     panel.grid.minor=element_blank(), 
                                     legend.position="none")
                  
                  ## Create hline data for cn plot
                  hline.dat <- data.table(chromosome=cnSeg$chrom,
                                          x=cnSeg$loc.start,
                                          xend=cnSeg$loc.end,
                                          y=cnSeg$seg.mean,
                                          yend=cnSeg$seg.mean)
                  
                  ## Define the hline plot
                  if (cnvType=="absolute" | cnvType == "logratio") {
                      hline <- geom_hline(yintercept = 2, lty=2)            
                  }
                  if (cnvType == "relative") {
                      hline <- geom_hline(yintercept = 0, lty=2)
                  }
                  segHLines <- geom_segment(data=hline.dat, aes(x=x, xend=xend, y=y, yend=yend), lty=1, col="red", size = 2)
                  
                  ## Define the facet
                  facet <- facet_grid(sample~chromosome, scale="free_x", space="fixed")
                  
                  ## Define the scales
                  scale_x <- scale_x_continuous(name="Position", expand=c(0,0), limits=c(chr$start, chr$end))
                  scale_y <- scale_y_continuous(name=paste0(cnvType, " Copy Number"))
                  
                  ## Build the plot
                  p1 <- ggplot(data=cnData, aes(x=position,y=cn)) + 
                      geom_point(color=plotAColor) + facet + hline + segHLines + 
                      scale_x + scale_y + plotALayers
                  
                  ##############################################################
                  ##### Build the somatic LOH plot #############################
                  ##############################################################
                  ## Print status message
                  if (verbose) {
                      message("Building somatic loh plot")
                  }
                  
                  ## Prepare somatic loh data to be plotted
                  normalDf <- somaticLohData[,c(1,2,4,5)]
                  colnames(normalDf) <- c("chromosome", "position", "VAF", "sample")
                  normalDf$Type <- "Normal"
                  tumorDf <- somaticLohData[,c(1,2,3,5)]
                  colnames(tumorDf) <- c("chromosome", "position", "VAF", "sample")
                  tumorDf$Type <- "Tumor"
                  allLohData <- rbind(normalDf, tumorDf)
                  
                  ## Define parameters of the plot
                  plotTheme <- theme(axis.ticks.x=element_blank(),
                                     axis.text.x=element_blank(),
                                     axis.ticks.y=element_blank(),
                                     panel.grid.major=element_blank(),
                                     panel.grid.minor=element_blank(), 
                                     legend.position="none")
                  
                  ## Create hline data for cn plot
                  hline.dat <- data.table(chromosome=lohSeg$chrom,
                                          x=lohSeg$loc.start,
                                          xend=lohSeg$loc.end,
                                          y=0.5+lohSeg$seg.mean,
                                          yend=0.5+lohSeg$seg.mean)
                  
                  ## Define the hline plot
                  h1 <- geom_hline(yintercept = 0.4, lty=2) 
                  h2 <- geom_hline(yintercept = 0.6, lty=2)
                  segHLines <- geom_segment(data=hline.dat, aes(x=x, xend=xend, y=y, yend=yend), lty=1, col="red", size = 2)
                  
                  ## Define the facet
                  facet <- facet_grid(sample~chromosome, scale="free_x", space="fixed")
                  
                  ## Define the scale
                  scale_x <- scale_x_continuous(name="Position", expand=c(0,0), limits=c(chr$start, chr$end))
                  scale_y <- scale_y_continuous(name="VAF", limits=c(0,1))
                  
                  ## Define the colors
                  color <- scale_color_manual(values=plotBColors)
                  
                  ## Build the plot
                  p2 <- ggplot(data=allLohData, aes(x=position,y=VAF, col=Type)) + 
                      geom_point(alpha=plotBAlpha) + facet + h1 + h2 + segHLines + 
                      scale_x + scale_y + color + plotBLayers

                  ##############################################################                  
                  ##### Build the germline LOH plot ############################
                  ##############################################################
                  ## Print status message
                  if (verbose) {
                      message("Building germline loh plot")
                  }
                  
                  ## Define parameters of the plot
                  plotTheme <- theme(panel.grid.major=element_blank(),
                                     panel.grid.minor=element_blank())
                  
                  
                  ## Define the facet
                  facet <- facet_grid(sample~chromosome, scale="free_x", space="fixed")
                  
                  ## Define the scale 
                  scale_x <- scale_x_continuous(name="Position", expand=c(0,0), limits=c(chr$start, chr$end))
                  scale_y <- scale_y_continuous(name="Normal VAF", breaks = c(0, 0.25, 0.5, 0.75, 1.0))
                  
                  ## Define the gradient
                  gradient <- scale_fill_gradient2(low=plotCColors[1], high=plotCColors[2],
                                                   limits=c(0, plotCLimits), oob=squish, trans="sqrt")
                  
                  ## Build the plot
                  p3 <- ggplot(data=germlineLohData, aes(x=position,y=normal_var_freq)) + 
                      geom_hex(binwidth=c((chr$end*.025/4),0.025)) + facet + gradient + scale_x + scale_y + 
                      plotTheme + plotCLayers
                  
                  
                  ##############################################################
                  ##### Combine all of the plots into 1 plot ###################
                  ##############################################################
                  ## Print status message
                  if (verbose) {
                      message("Combining cn, somatic loh, and germline loh plots")
                  }
                  
                  ## Obtain the max width for relevant plots
                  cnPlot <- ggplotGrob(p1)
                  somaticLohPlot <- ggplotGrob(p2)
                  germlineLohPlot <- ggplotGrob(p3)
                  plotList <- list(cnPlot, somaticLohPlot, germlineLohPlot)
                  
                  plotList <- plotList[lapply(plotList, length) > 0]
                  plotWidths <- lapply(plotList, function(x) x$widths)
                  maxWidth <- do.call(grid::unit.pmax, plotWidths)
                  
                  ## Set the widths for all plots
                  for (i in 1:length(plotList)) {
                      plotList[[i]]$widths <- maxWidth
                  }
                  
                  ## Arrange the final plot
                  finalPlot <- do.call(gridExtra::arrangeGrob, c(plotList, list(ncol=1, heights=sectionHeights)))
                  plot(finalPlot)
                  return(finalPlot)
              }, 
              chrData=chrData, cnvType=cnvType, plotAColor=plotAColor, plotALayers=plotALayers, 
              plotBAlpha=plotBAlpha, plotBColors=plotBColors, plotBLayers=plotBLayers,
              somaticLohCutoff=somaticLohCutoff, plotCLimits=plotCLimits,
              plotCColors=plotCColors, plotCLayers=plotCLayers, sectionHeights=sectionHeights, 
              verbose=verbose)
              
              return(combinedPlots)
          })
