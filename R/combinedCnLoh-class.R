################################################################################
##################### Public/Private Class Definitions #########################

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Public Class !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#' Class cnLoh
#' 
#' An S4 class for the cn, somatic loh, and germline loh plots
#' @name cnLoh-class
#' @rdname cnLoh-class
#' @slot cnData data.table object for cn plot
#' @slot cnPlot gtable object for the cn plot
#' @slot somaticLohData data.table object for the somatic loh plot
#' @slot somaticLohPlot gtable object for the somatic loh plot
#' @slot germlineLohData data.table object for the germline loh plot
#' @slot germlineLohData gtable object for the germline loh plot
#' @exportClass cnLoh
#' @import methods
#' @importFrom data.table data.table
#' @importFrom gtable gtable
methods::setOldClass("gtable")
setClass(
    Class="cnLoh",
    representation=representation(cnData="data.table",
                                  cnPlot="gtable",
                                  somaticLohData="data.table",
                                  somaticLohPlot="gtable",
                                  germlineLohData="data.table",
                                  germlineLohPlot="gtable",
                                  Grob="gtable"),
    validity=function(object){
        
    }
)

#' Constructor for the cnLoh class
#' 
#' @name cnLoh
#' @rdname cnLoh-class
#' @param input Object of class cnLohDataFormat
#' @param samples Character vector specifying samples to plot. If not NULL
#' all samples in "input" not specified with this parameter are removed.
#' @param chromosomes Character vector specifying chromosomes to plot. If not NULL
#' all chromosomes in "input" not specified with this parameter are removed.
#' @param BSgenome Object of class BSgenome to extract genome wide chromosome 
#' coordinates
#' @export
cnLoh <- function(cnInput, lohInput, samples, chromosomes, BSgenome, windowSize,
                  step, getHeterozygousLohCalls, plotAColor, plotALayers, plotBAlpha,
                  somaticLohCutoff, plotBColors, plotBLayers,
                  plotCLimits, plotCColors, plotCLayers, 
                  sectionHeights, verbose) {
    
    ## Check each of the input parameters
    cnLohInputParameters <- checkCombinedCnLohInputParameters(cnInput=cnInput, lohInput=lohInput, samples=samples,
                                                              chromosomes=chromosomes, BSgenome=BSgenome,
                                                              windowSize=windowSize, step=step, 
                                                              getHeterozygousLohCalls=getHeterozygousLohCalls,
                                                              somaticLohCutoff=somaticLohCutoff,
                                                              plotAColor=plotAColor, plotALayers=plotALayers,
                                                              plotBAlpha=plotBAlpha, plotBColors=plotBColors,
                                                              plotBLayers=plotBLayers, plotClimits=plotCLimits,
                                                              plotCColors=plotCColors, plotCLayers=plotCLayers,
                                                              sectionHeights=sectionHeights, verbose=verbose)
    
    ## Obtain cn, somatic loh, and germline loh datasets to plot
    cnLohDataset <- cnLohData(cnInput=cnInput, lohInput=lohInput, 
                              samples=checkCombinedCnLohInputParameters@samples, 
                              chromosomes=checkCombinedCnLohInputParameters@chromosomes,
                              BSgenome=checkCombinedCnLohInputParameters@BSgenome, 
                              windowSize=checkCombinedCnLohInputParameters@windowSize,
                              step=checkCombinedCnLohInputParameters@step, 
                              normal=checkCombinedCnLohInputParameters@getHeterozygousLohCalls, 
                              verbose=checkCombinedCnLohInputParameters@verbose)
    
    ## Generate the cn, somatic LOH, and germline LOH plots
    plots <- cnLohPlots(object=cnLohDataset, 
                        somaticLohCutoff=checkCombinedCnLohInputParameters@somaticLohCutoff, 
                        plotAColor=checkCombinedCnLohInputParameters@plotAColor, 
                        plotALayers=checkCombinedCnLohInputParameters@plotALayers, 
                        plotBAlpha=checkCombinedCnLohInputParameters@plotBAlpha, 
                        plotBColors=checkCombinedCnLohInputParameters@plotBColors,
                        plotBNormalColor=checkCombinedCnLohInputParameters@plotBNormalColor, 
                        plotBLayers=checkCombinedCnLohInputParameters@plotBLayers, 
                        plotCLimits=checkCombinedCnLohInputParameters@plotCLimits, 
                        plotCLowColor=checkCombinedCnLohInputParameters@plotCColors, 
                        plotCLayers=checkCombinedCnLohInputParameters@plotCLayers, 
                        verbose=checkCombinedCnLohInputParameters@verbose)
    
    ## Arrange all of the plots together
    Grob <- arrangeCnLohPlots(object=plots, 
                              sectionHeights=checkCombinedCnLohInputParameters@sectionHeights, 
                              verbose=checkCombinedCnLohInputParameters@verbose)
    
    ## Initialize the object
    new("cnLoh", cnData=getData(cnLohDataset, index=1), cnPlot=getGrob(plots, index=1),
        somaticLohData=getData(cnLohDataset, index=2), somaticLohPlot=getGrob(plots, index=2),
        germlineLohData=getData(cnLohDataset, index=3), germlineLohPlot=getGrob(plots, index=3),
        Grob=Grob)
}

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Private Classes !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

#' Private Class cnLohInputParameters
#' 
#' An S4 class to check input parameters of the combinedCnLoh function
#' @name cnLohInputParameters-class
#' @noRd
setClass("cnLohInputParameters",
         representation=representation(samples="character",
                                       chromosomes="character", BSgenome="BSgenome",
                                       windowSize="numeric", step="numeric", 
                                       getHeterozygousLohCalls="logical",
                                       somaticLohCutoff="numeric",
                                       plotAColor="character", plotALayers="list",
                                       plotBAlpha="numeric", plotBColors="character",
                                       plotBLayers="list", plotClimits="numeric",
                                       plotCColors="character", plotCLayers="list",
                                       sectionHeights="numeric", verbose="logical"), 
         validity=function(object){
             
         })

#' Constructor for the cnLohInputParameters class
#' 
#' @name cnLohInputParameters
#' @rdname cnLohInputParameters-class
#' @noRd
checkCombinedCnLohInputParameters <- function(cnInput, lohInput, samples, chromosomes, BSgenome,
                                              windowSize, step, getHeterozygousLohCalls,
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
    
    ##### TODO: Check the samples parameter ##### 
    ## Check is samples is NULL
    if (is.null(samples)) {
        samples <- unique(object@sample$sample)
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
        discrepantSamples <- paste(samples[which(samples %nin% unique(object@sample$sample))], collapse=", ")
        if (length(discrepantSamples) > 0 & discrepantSamples != "") {
            memo <- paste0("The desired samples: ", discrepantSamples, " are not found ",
                           "in the SV dataset. Available sample names include: ", 
                           paste(unique(object@sample$sample), collapse=", "), ". ",
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
    
    ##### Check the windowSize parameter #####
    if (is.null(windowSize)) {
        windowSize <- 2500000
        memo <- paste0("windowSize parameter cannot be NULL. Setting the windowSize value to 2500000.")
        message(memo)
    }
    ## Check if windowSize is numeric
    if (!is.numeric(windowSize)) {
        memo <- paste0("windowSize variable not of the numeric class. Attempting to coerce.")
        windowSize <- as.numeric(windowSize)
        message(memo)
    }
    
    ##### Check the step parameter #####
    if (is.null(step)) {
        step <- 1000000
        memo <- paste0("step parameter cannot be NULL. Setting the step value to 1000000.")
        message(memo)
    }
    ## Check if step is numeric
    if (!is.numeric(step)) {
        memo <- paste0("step variable not of the numeric class. Attempting to coerce.")
        step <- as.numeric(step)
        message(memo)
    }
    ## Check to see if step is greater than windowSize
    if (step > windowSize) {
        memo <- paste("Step value is greater than windowSize. Make sure that the step value is 
                      at most equal to the WindowSize. Using default values for both parameters.")
        warning(memo)
        step <- 1000000
        windowSize <- 2500000
        
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
        windowSize=windowSize, step=step, getHeterozygousLohCalls=getHeterozygousLohCalls,
        somaticLohCutoff=somaticLohCutoff, plotAColor=plotAColor, plotALayers=plotALayers,
        plotBAlpha=plotBAlpha, plotBColors=plotBColors, plotBLayers=plotBLayers, 
        plotClimits=plotCLimits, plotCColors=plotCColors, plotCLayers=plotCLayers,
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
cnLohData <- function(cnInput, lohInput, samples, chromosomes, BSgenome,
                      windowSize, step, normal, verbose=FALSE) {

    ############################################################################
    #################### Prepare copy number variant dataset ###################
    ## Obtain raw cnv data
    cnData <- getCnvData(object=cnInput, verbose=verbose)
    
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
    lohData <- getLohData(object=lohInput, verbose=verbose, lohSpec=TRUE, germline=FALSE)
    
    ## Subset loh data by chromosome
    lohData <- chrSubset(object=lohData, chromosomes=chromosomes, verbose=verbose)
    
    ## Subset loh data by sample
    lohData <- sampleSubset(object=lohData, samples=samples, verbose=verbose)
    
    ## Produce data.table with window position data
    windowData <- getLohSlidingWindow(object=lohData, step=step, windowSize=windowSize,
                                      verbose=verbose)
    
    ## Perform loh calculations on each chromosome and samples within each window
    lohAbsDiff <- getLohCalculation(object=lohData, windowData=windowData, normal=normal,
                                    verbose=verbose)
    
    ## Calculate avg loh for overlapping regions
    lohAbsDiffOverlap <- rbindlist(getLohStepCalculation(object=lohAbsDiff,
                                                         step=step, verbose=verbose))
    
    ## Obtain loh segmentation dataset
    lohSegmentation <- getLohSegmentation(object=lohAbsDiffOverlap, verbose=verbose)
    
    ## Remove gaps 
    lohSegmentation <- removeGapsSegmentation(object=lohSegmentation, chrData=chrData,
                                              verbose=verbose)
    
    ############################################################################
    ##################### Prepare germline loh dataset #########################
    ## Obtain germlineloh data by chromosome
    germlineLohData <- getLohData(object=lohInput, verbose=TRUE, lohSpec=FALSE, germline=TRUE)
    
    ## Subset loh data by chromosome
    germlineLohData <- chrSubset(object=germlineLohData, chromosomes=chromosomes, verbose=verbose)
    
    ## Subset loh data by sample
    germlineLohData <- sampleSubset(object=germlineLohData, samples=samples, verbose=verbose)
    
    ## Initialize the object
    new("cnLohData", rawCnData=cnData, segCnData=cnSegmentation, 
        rawLohData=lohData, segLohData=lohSegmentation, rawGermlineLohData=germlineLohData,
        chrData=chrData)
} 

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
         representation=representation(cnPlot="gtable",
                                       somaticLohPlot="gtable",
                                       germlineLohPlot="gtable"),
         validity=function(object) {
             
         })

#' Constructor for cnLohPlots class
#' 
#' @rdname cnLohPlots-class
#' @name cnLohPlots
#' @param object Object of class data.table
#' @importFrom gtable gtable
#' @import ggplot2
#' @noRd 
cnLohPlots <- function(object, plotAColor, plotALayers, 
                       somaticLohCutoff, plotBAlpha, plotBTumorColor, plotBNormalColor, plotBLayers, 
                       plotCLimits, plotCLowColor, plotCHighColor, 
                       plotCLayers, verbose) {
    
    ## Build the cn plot
    cnPlot <- buildCnPlot(object=object, plotAColor=plotAColor, plotALayers)
    
    ## Build the somatic loh plot
    somaticLohPlot <- buildSomaticLohPlot(object=object, somaticLohCutoff=somaticLohCutoff,
                                          plotBAlpha=plotBAlpha,
                                          plotBTumorColor=plotBTumorColor, 
                                          plotBNormalColor=plotBNormalColor,
                                          plotBLayers=plotBLayers, verbose=verbose)
    
    ## Build the germline loh plot
    germlineLohPlot <- buildGermlineLohPlot(object=object, 
                                            plotCLimits=plotCLimits, plotCLowColor=plotCLowColor,
                                            plotCHighColor=plotCHighColor, plotCLayers=plotCLayers,
                                            verbose=verbose)
    
    ## Initialize the object
    new("cnLohPlots", cnPlot=cnPlot, somaticLohPlot=somaticLohPlot, germlineLohPlot=germlineLohPlot)
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
    definition=function(object, ...){
        mainPlot <- getGrob(object, index=4)
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
              
              ## Split object by sample
              segDfTemp <- split(object, list(as.character(object$sample)))
              
              ## Perform segmentation
              segmentationDF <- rbindlist(lapply(segDfTemp, function(x) {
                  cnSeg <- CNA(genomdat=as.numeric(x$cn), chrom=x$chromosome,
                                maploc=x$position, data.type="logratio", sampleid = unique(x$sample))
                  
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
              chrList <- as.list(as.character(unique(object$chrom)))
              segs <- rbindlist(lapply(chrList, function(x, object, chrData) {
                  df <- object[chrom==x]
                  for (i in 1:(nrow(df) - 1)) {
                      ## Don't merge segments if they are far apart
                      if ((df$loc.start[i+1]-df$loc.end[i]) < 5000000) {
                          half <- floor((df$loc.end[i] + df$loc.start[i+1])/2)
                          df$loc.end[i] <- half
                          df$loc.start[i+1] <- half + 1
                      }
                  }
                  return(df)
              }, object=object))
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

##########################################################################
##### Function to generate window position data for loh calculations #####
#' @rdname getLohSlidingWindow-methods
#' @name getLohSlidingWindow
#' @param object of class data.table 
#' @param step integer specifying the step size between the start position of
#' each window
#' @param windowSize integer specifying the window size for loh calcuations
#' @return Data.table with window start/stop positions
#' @aliases getLohSlidingWindow
setMethod(f="getLohSlidingWindow",
          signature="data.table",
          definition=function(object, step, windowSize, ...){
              if (verbose) {
                  message("Calcuating window sizes for loh calcluations on all chromosomes in each individual sample")
              }
              
              ## Obtain lists for each sample and chromosome
              out <- split(object, list(as.character(object$chromosome),
                                        as.character(object$sample)))
              
              ## Obtain the window position values
              window <- lapply(out, function(x, step, windowSize) {
                  ## Get the min and max position on the chromosome
                  min <- integer()
                  max <- integer()
                  window_stop_1 <- integer()
                  window_num <- integer()
                  min <- as.integer(min(as.numeric(as.character(x$position))))
                  max <- as.integer(max(as.numeric(as.character(x$position))))
                  ## Get the end of the first window position
                  window_stop_1 <- min+windowSize
                  ## Calculate the number of windows necessary
                  num <- as.integer((max-min)/step)
                  num <- as.vector(1:num)
                  window_data_start <- vector()
                  window_data_stop <- vector()
                  
                  ## Calculate exact window positions
                  window_data <- lapply(num, function(x){
                      window_data_start[x] <- as.integer(min+(step*(x-1)))
                      window_data_stop[x] <- as.integer(window_stop_1+(step*(x-1)))
                      window_data <- data.table(cbind(window_data_start[x], window_data_stop[x]))
                      return(window_data)
                  })
                  window_data <- rbindlist(window_data)
                  # Get window positions whose values are below max & set max as the 
                  # final window position (end of the chromosome)
                  colnames(window_data) <- c("window_start", "window_stop")
                  window_final <- window_data[window_data$window_stop <= max,]
                  window_final[nrow(window_final), 2] <- max
                  ## Put in the chromosome 
                  window_final$chromosome <- as.character(x$chromosome[1])
                  return(window_final)
              }, 
              step = step, windowSize = windowSize)
              
              return(window)
          })

###############################################################
##### Function to perform loh calcluations in each window #####
#' @rdname getLohCalculation-methods
#' @name getLohCalculation
#' @param object of class data.table 
#' @param window_data of class data.table 
#' @param normal integer specifying normal vaf
#' @aliases getLohCalculation
setMethod(f="getLohCalculation", 
          signature="data.table",
          definition=function(object, windowData, normal, verbose, ...) {
              
              ## Print status message
              if (verbose) {
                  message("Calculating absolute mean difference between t/n VAF at each coordinate provided.")
              }
              
              ## Split object for each unqiuq sample-chr combination
              object <- split(object, list(as.character(object$chromosome),
                                           as.character(object$sample)))
              
              ## Separate out sample and window data by chromosome name
              df <- lapply(object, function(sampleData, window, 
                                            normal) {
                  chromosome <- as.character(sampleData[1,chromosome])
                  sample <- as.character(sampleData[1,sample])
                  chromosome.sample <- paste("\\b", paste(chromosome, sample, sep = "."), "\\b", sep = "")
                  window <- windowData[[grep(chromosome.sample, names(windowData))]]
                  ## For each window position, get the vaf data that falls 
                  ## within that window
                  dataset <- rbindlist(apply(window, 1, function(x, sampleData, normal){
                      ## Determine which value to use for normal
                      if (normal==FALSE) {
                          normal <- 0.5
                      }
                      if (normal == TRUE) {
                          normal <- round(sampleData[,mean(normal_var_freq)], 
                                          digits=3)
                      }
                      
                      w_start <- as.numeric(as.character(x[1]))
                      w_stop <- as.numeric(as.character(x[2]))
                      ## Filter out vaf data outside the window
                      filtered_data <- sampleData[position >= w_start &
                                                      position <= w_stop]
                      
                      ## Peroform loh calclulation to obtain avg loh in the 
                      ## window's frame
                      loh_calc_avg <- mean(abs(as.numeric(as.character(
                          filtered_data$tumor_var_freq)) - normal))
                      ## If no coordinates are found within the window,
                      ## make as NA
                      if (is.na(loh_calc_avg)) {
                          loh_calc_avg <- NA
                          w_start <- NA
                          w_stop <- NA
                      }
                      filtered_data$loh_diff_avg <- loh_calc_avg
                      filtered_data$window_start <- w_start
                      filtered_data$window_stop <- w_stop
                      return(filtered_data)
                  }, 
                  sampleData=sampleData, normal=normal))
                  dataset <- na.omit(dataset, cols = c("loh_diff_avg", 
                                                       "window_start", 
                                                       "window_stop"))
                  return(dataset)
              }, window=windowData, normal=normal)
              return(df)
          })

#######################################################################
##### Function to perform loh calcluations in overlapping windows #####
#' @rdname getLohStepCalculation-methods
#' @name getLohStepCalculation
#' @param object of class data.table
#' @param step integer 
#' @aliases getLohStepCalculation
setMethod(f = "getLohStepCalculation",
          signature="list",
          definition=function(object, step, ...) {
              
              ## Print status message
              if (verbose) {
                  message("Calculating loh in overlapping windows")
              }
              step_loh_calc <- lapply(object, function(x, step) {
                  ## Get the sample and chromosome information
                  sample <- unique(x$sample)
                  chromosome <- unique(x$chromosome)
                  
                  ## Obtain boundaries for each step-sized window
                  start <- unique(x$window_start)
                  stop <- c(start[-1], max(x$window_stop))
                  step_boundaries <- data.table(chromosome=chromosome, start=start, stop=stop)
                  step_boundaries$sample <- sample
                  
                  ## Get the average loh within each step-sized window
                  loh_df <- x
                  loh_step_avg <- apply(step_boundaries, 1, function(x, loh_df_data) {
                      start <- as.numeric(as.character(x[2]))
                      stop <- as.numeric(as.character(x[3]))
                      step_df <- loh_df_data[position >= start & 
                                                 position < stop]
                      if (nrow(step_df) == 0) {
                          loh_step_avg <- 0
                      }
                      if (nrow(step_df) > 0) {
                          loh_step_avg <- mean(step_df$loh_diff_avg)
                      }
                      return(loh_step_avg)
                  }, loh_df_data=loh_df)
                  step_boundaries$loh_step_avg <- loh_step_avg
                  return(step_boundaries)
              }, step=step)
              return(step_loh_calc)
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
              segDfTemp <- split(object, list(as.character(object$sample)))
              segmentationDf <- rbindlist(lapply(segDfTemp, function(x){
                  x$midpoint <- floor((as.numeric(x$start) + as.numeric(x$stop))/2)
                  lohSeg <- CNA(genomdat = as.numeric(x$loh_step_avg), chrom = x$chromosome,
                                maploc = x$midpoint, data.type = "binary", sampleid = unique(x$sample))
                  lohSeg <- segment(lohSeg)
                  lohSeg <- lohSeg$output
                  return(lohSeg)
              }))
              
              return(segmentationDf)
          })

########################################
##### Function to generate cn plot #####
#' @rdname buildCnPlot-methods
#' @name buildCnPlot
#' @aliases buildCnPlot
#' @param object of class data.table
#' @noRd
setMethod(f="buildCnPlot",
          signature="cnLohData",
          definition=function(object, plotAColor, plotALayers, ...){
              
              ## Print status message
              if (verbose) {
                  message("Building cnv plot")
              }
              
              ## Perform quality checks on the input variables
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
              
              ## Separate datasets
              rawCnData <- object@rawCnData
              segCnData <- object@segCnData
              
              ## Define parameters of the plot
              plotTheme <- theme(axis.ticks.x=element_blank(),
                                 axis.text.x=element_blank(),
                                 axis.ticks.y=element_blank(),
                                 panel.grid.major=element_blank(),
                                 panel.grid.minor=element_blank(), 
                                 legend.position="none")
              
              ## Create hline data for cn plot
              hline.dat <- data.table(chromosome=segCnData$chrom,
                                      x=segCnData$loc.start,
                                      xend=segCnData$loc.end,
                                      y=segCnData$seg.mean,
                                      yend=segCnData$seg.mean)
              
              ## Define the hline plot
              hline <- geom_hline(yintercept = 2, lty=2)            
              segHLines <- geom_segment(data=hline.dat, aes(x=x, xend=xend, y=y, yend=yend), lty=1, col="red", size = 2)
              
              ## Define the facet
              facet <- facet_grid(sample~chromosome, scale="free_x", space="fixed")
              
              ## Define the scales
              scale_x <- scale_x_continuous(name="Position", expand=c(0,0))
              scale_y <- scale_y_continuous(name="Absolute Copy Number")
              
              ## Build the plot
              p1 <- ggplot(data=rawCnData, aes(x=position,y=cn)) + 
                  geom_point(color=plotAColor) + facet + hline + segHLines + 
                  scale_x + scale_y + plotALayers
                 
              ## Convert to grob
              cnPlotGrob <- ggplotGrob(p1)
              plot(cnPlotGrob)
              return(cnPlotGrob)
          })

#################################################
##### Function to generate somatic loh plot #####
#' @rdname buildSomaticLohPlot-methods
#' @name buildSomaticLohPlot
#' @aliases buildSomaticLohPlot
#' @param object of class data.table
#' @noRd
setMethod(f="buildSomaticLohPlot",
          signature="cnLohData",
          definition=function(object, somaticLohCutoff, plotBAlpha, plotBTumorColor, plotBNormalColor,
                              plotBLayers, ...){
              
              ## Print status message
              if (verbose) {
                  message("Building somatic loh plot")
              }
              
              ## Separate datasets
              segLohData <- object@segLohData
              segLohData <- segLohData[seg.mean > somaticLohCutoff]
              rawLohData <- object@rawLohData
              
              ## Prepare loh data to be plotted
              normalDf <- rawLohData[,c(1,2,4,5)]
              colnames(normalDf) <- c("chromosome", "position", "VAF", "sample")
              normalDf$Type <- "Normal"
              tumorDf <- rawLohData[,c(1,2,3,5)]
              colnames(tumorDf) <- c("chromosome", "position", "VAF", "sample")
              tumorDf$Type <- "Tumor"
              rawLohData <- rbind(normalDf, tumorDf)
              
              ## Define parameters of the plot
              plotTheme <- theme(axis.ticks.x=element_blank(),
                                 axis.text.x=element_blank(),
                                 axis.ticks.y=element_blank(),
                                 panel.grid.major=element_blank(),
                                 panel.grid.minor=element_blank(), 
                                 legend.position="none")
              
              ## Create hline data for cn plot
              hline.dat <- data.table(chromosome=segLohData$chrom,
                                      x=segLohData$loc.start,
                                      xend=segLohData$loc.end,
                                      y=0.5+segLohData$seg.mean,
                                      yend=0.5+segLohData$seg.mean)
              
              ## Define the hline plot
              h1 <- geom_hline(yintercept = 0.4, lty=2) 
              h2 <- geom_hline(yintercept = 0.6, lty=2)
              segHLines <- geom_segment(data=hline.dat, aes(x=x, xend=xend, y=y, yend=yend), lty=1, col="red", size = 2)
              
              ## Define the facet
              facet <- facet_grid(sample~chromosome, scale="free_x", space="fixed")
              
              ## Define the scale
              scale_x <- scale_x_continuous(name="Position", expand=c(0,0))
              scale_y <- scale_y_continuous(name="VAF", limits=c(0,1))
              
              ## Define the colors
              color <- scale_color_manual(values=c(plotBNormalColor, plotBTumorColor))

              ## Build the plot
              p1 <- ggplot(data=rawLohData, aes(x=position,y=VAF, col=Type)) + 
                  geom_point(alpha=plotBAlpha) + facet + h1 + h2 + segHLines + 
                  scale_x + scale_y + color + plotBLayers
              
              ## Convert to grob
              somaticLohPlotGrob <- ggplotGrob(p1)
              plot(somaticLohPlotGrob)
              return(somaticLohPlotGrob)
          })

##################################################
##### Function to generate germline loh plot #####
#' @rdname buildGermlineLohPlot-methods
#' @name buildGermlineLohPlot
#' @aliases buildGermlineLohPlot
#' @param object of class data.table
#' @noRd
setMethod(f="buildGermlineLohPlot",
          signature="cnLohData",
          definition=function(object, plotCLimits, plotCLowColor,
                              plotCHighColor, plotCLayers, verbose=verbose){
              
              ## Print status message
              if (verbose) {
                  message("Building germline loh plot")
              }
              
              ## Separate datasets
              germlineLohData <- object@rawGermlineLohData
              chrData <- object@chrData
              
              ## Define parameters of the plot
              plotTheme <- theme(panel.grid.major=element_blank(),
                                 panel.grid.minor=element_blank())
              
            
              ## Define the facet
              facet <- facet_grid(sample~chromosome, scale="free_x", space="fixed")
              
              ## Define the scale 
              scale_x <- scale_x_continuous(name="Position", expand=c(0,0))
              scale_y <- scale_y_continuous(name="Normal VAF", breaks = c(0, 0.25, 0.5, 0.75, 1.0))
              
              ## Define the gradient
              gradient <- scale_fill_gradient2(low=plotCLowColor, high=plotCHighColor,
                                               limits=plotCLimits, oob=squish, trans="sqrt")

              ## Build the plot
              p1 <- ggplot(data=germlineLohData, aes(x=position,y=normal_var_freq)) + 
                  geom_hex(binwidth=c((chrData$end[1]*.025/4),0.025)) + facet + gradient + scale_x + scale_y + 
                  plotTheme + plotCLayers
              
              ## Convert to grob
              germlineLohPlot <- ggplotGrob(p1)
              plot(germlineLohPlot)
              return(germlineLohPlot)
              
          }) 

#########################################################
##### Function to arrange lohSpec and lohFreq plots #####
#' @rdname arrangeCnLohPlots-methods
#' @name arrangeCnLohPlots
#' @param object of class cnLohData
#' @aliases arrangeCnLohPlots
#' @noRd
setMethod(f="arrangeCnLohPlots",
          signature="cnLohPlots",
          definition=function(object, sectionHeights, verbose, ...) {
              
              ## Print status message
              if (verbose) {
                  message("Combining cn, somatic loh, and germline loh plots")
              }
              
              ## Grab the data we need
              plotA <- object@cnPlot
              plotB <- object@somaticLohPlot
              plotC <- object@germlineLohPlot
              
              ## obtain the meax width for relevant plots
              plotList <- list(plotA, plotB, plotC)
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
          })