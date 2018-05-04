################################################################################
##################### Public/Private Class Definitions #########################

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Public Class !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#' Class lohSpec
#' 
#' An S4 class for the lohSpec plot object
#' @name lohSpec
#' @rdname lohSpec-class
#' @slot lohFreq_plot gtable object for the lohFreq plot
#' @slot lohSpec_plot gtable object for the lohSpec plot
#' @slot lohData data.table object soring loh data with column names: sample, 
#' chromosome, position, t_vaf, n_vaf. 
#' @exportClass lohSpec
#' @importFrom data.table data.table
#' @importFrom gtable gtable
methods::setOldClass("gtable")
setClass(
    Class="lohSpec",
    representation=representation(lohFreq_plot="gtable",
                                  lohSpec_plot="gtable",
                                  Grob="gtable",
                                  lohData="data.table"),
    validity = function(object) {
        
    }
)

#' Constructor for the lohSpec class
#' 
#' @name lohSpec
#' @rdname lohSpec-class
#' @param input Object of class VarScan.
#' @param samples Character vector specifying samples to plot. If not NULL
#' all samples in "input" not specified with this parameter are removed.
#' @param chromosomes Character vector specifying chromosomes to plot. If not NULL
#' all chromosomes in "input" not specified with this parameter are removed. If 
#' autosomes, it will plot the non-sex chromosomes. Defaults to autosomes. 
#' @param BSgenome Object of class BSgenome to extract genome wide chromosome 
#' coordinates
#' @param step Integer value specifying the step size (i.e. the number of base
#' pairs to move the window). required when method is set to slide
#' (see details). Defaults to 1000000
#' @param windowSize Integer value specifying the size of the window in base
#' pairs in which to calculate the mean Loss of Heterozygosity (see details). 
#' Defaults to 2500000. 
#' @param normalVAF Boolean specifiying what value to use for normal VAF when 
#' calcualting average LOH difference. Defaults to .50\% if FALSE. 
#' If TRUE, will use average normal VAF in each individual sample as value 
#' to calculate LOH.
#' @param gradientMidpoint Integer value specifying the midpoint of the loh 
#' color gradient, which ranges from 0-0.5. Defaules to 0.2. 
#' @param gradientColors Character vector with valid colors specifying the 
#' gradient to visualize the loh spectrum. Defaults to "#ffffff", "#b2b2ff", and
#'  "#000000"
#' @param plotAType Character vector (either "proportion" or "frequency") that 
#' directrs the function to plot the proportion or the frequency of samples with 
#' LOH in a specific region across the entire cohort. Defaults to proportion.  
#' @param plotALohCutoff Integer specifying the VAF cutoff to consider LOH events
#' @param plotAColor Character vector specifying the color for loh frequency bars
#' @param plotALayers List of ggplot2 layers to be passed to loh frequency plot
#' @param plotBLayers List of ggplot2 layers to be passed to loh heatmap plot
#' @param sectionHeights Integer vector specifying the relative heights for each of the plots
#' @param verbose Boolean specifying if status messages should be reported
#' @export
LohSpec <- function(input, getHeterozygousCalls=TRUE, chromosomes="autosomes", samples=NULL, 
                    BSgenome=BSgenome, step=1000000, windowSize=2500000, 
                    normalVAF=FALSE, gradientMidpoint=.2, gradientColors=c("#ffffff", "#b2b2ff", "#000000"),
                    plotAType="proportion", plotALohCutoff=0.1, plotAColor="#98F5FF",
                    plotALayers=NULL, plotBLayers=NULL, sectionHeights=c(0.25, 0.75), verbose=FALSE){
    
    ## Constructor to check input parameters for the lohSpec function
    lohInputParameters <- checkLohInputParameters(object=object, getHeterozygousCalls=getHeterozygousCalls,
                                                  chromosomes=chromosomes, sample=sample,
                                                  BSgenome=BSgenome, step=step, 
                                                  windowSize=windowSize, normalVAF=normalVAF, 
                                                  gradientMidpoint=gradientMidpoint,
                                                  gradientColors=gradientColors,
                                                  plotAType=plotAType, plotALohCutoff=plotALohCutoff,
                                                  plotAColor=plotAColor, plotALayers=plotALayers,
                                                  plotBLayers=plotBLayers, sectionHeights=sectionHeights, 
                                                  verbose=verbose)
    
    ## Calculate all data for plots
    lohDataset <- lohData(object=input, getHeterozygousCalls=lohInputParameters@getHeterozygousCalls, 
                          chromosomes=lohInputParameters@chromosomes, 
                          samples=lohInputParameters@sample, 
                          BSgenome=lohInputParameters@BSgenome, 
                          step=lohInputParameters@step, 
                          plotALohCutoff=lohInputParameters@plotALohCutoff,
                          windowSize=lohInputParameters@windowSize, 
                          normal=lohInputParameters@normalVAF, 
                          verbose=lohInputParameters@verbose)
    
    ## Initialize the lohSpecPlots object
    plots <- lohSpecPlots(object=lohDataset, 
                          plotALohCutoff=lohInputParameters@plotALohCutoff, 
                          plotAType=lohInputParameters@plotAType, 
                          plotAColor=lohInputParameters@plotAColor, 
                          plotALayers=lohInputParameters@plotALayers, 
                          plotBLayers=lohInputParameters@plotBLayers, 
                          gradientMidpoint=lohInputParameters@gradientMidpoint, 
                          gradientColors=lohInputParameters@gradientColors,
                          verbose=lohInputParameters@verbose)
    
    ## Arrange freq and spectrum plots 
    Grob <- arrangeLohPlots(object=plots, sectionHeights=sectionHeights, 
                            verbose=verbose)
    
    ## Initialize the object
    new("lohSpec", lohFreq_plot=getGrob(plots, index=1), lohSpec_plot=getGrob(plots, index=2), 
        lohData=getData(lohDataset, name="primaryData"), Grob=Grob)
}
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Private Classes !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

#' Private Class lohInputParameters
#' 
#' An S4 class to check input parameters of the StructuralVariant function
#' @name lohInputParameters-class
#' @noRd
setClass("lohInputParameters",
         representation=representation(getHeterozygousCalls="logical", chromosomes="character", 
                                       sample="character", BSgenome="BSgenome", 
                                       step="numeric", windowSize="numeric", 
                                       normalVAF="logical", gradientMidpoint="numeric", 
                                       gradientColors="character",
                                       plotAType="character", plotALohCutoff="numeric", 
                                       plotAColor="character", plotALayers="list", 
                                       plotBLayers="list", sectionHeights="numeric", 
                                       verbose="logical"), 
         validity=function(object){
             
         })

#' Constructor for the lohInputParameters class
#' 
#' @name lohInputParameters
#' @rdname lohInputParameters-class
#' @noRd
checkLohInputParameters <- function(object, getHeterozygousCalls, chromosomes, sample,
                                    BSgenome, step, windowSize, normalVAF, 
                                    gradientMidpoint, gradientColors,
                                    plotAType, plotALohCutoff, plotAColor, plotALayers,
                                    plotBLayers, sectionHeights, verbose) {
    ##### Check the verbose parameter #####
    ## Check to see if verbose is a boolean
    if (!is.logical(verbose) | is.null(verbose)) {
        memo <- paste0("The verbose parameter is not a boolean (T/F). Coercing verbose to be FALSE...")
        message(memo)
        verbose <- FALSE
    }
    
    ##### Check the getHeterozygousCalls parameter #####
    if (!is.logical(getHeterozygousCalls) | is.null(getHeterozygousCalls)) {
        memo <- paste0("The getHeterozygousCalls parameter is not a boolean (T/F). Coercing getHeterozygousCalls to be TRUE...")
        message(memo)
        getHeterozygousCalls <- FALSE
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
    
    ##### Check the samples parameter #####
    ## Check is sample is NULL
    if (is.null(sample)) {
        sample <- unique(object@sample$sample)
        sample <- factor(sample, levels=gtools::mixedsort(sample))
        memo <- paste0("Sample parameter cannot be NULL. All samples will be plotted.")
        message(memo)
    }
    ## Check if sample is a character vector
    if (!is.character(sample)) {
        memo <- paste0("sample variable not of the character class. Attempting to coerce.")
        sample <- as.character(sample)
        message(memo)
    }
    
    ## Check if the designated samples is in the sv dataset
    if (!is.null(sample)) {
        `%nin%` = Negate(`%in%`)
        discrepantSamples <- paste(sample[which(sample %nin% unique(object@sample$sample))], collapse=", ")
        if (length(discrepantSamples) > 0 & discrepantSamples != "") {
            memo <- paste0("The desired samples: ", discrepantSamples, " are not found ",
                           "in the SV dataset. Available sample names include: ", 
                           paste(unique(object@sample$sample), collapse=", "), ". ",
                           "Please designate valid sample names.")
            stop(memo)
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
    
    ##### Check the normalVAF parameter #####
    if (!is.logical(normalVAF) | is.null(normalVAF)) {
        memo <- paste0("normalVAF parameter should be a boolean (T/F). True if ", 
                       "user wants to use normal VAF from varscan to identify tumor/normal LOH difference. ", 
                       "False if user wants to use 0.5 to identify tumor/normal LOH difference. Setting normalVAF ",
                       "to TRUE.")
        message(memo)
        normalVAF <- TRUE
    }
    
    ##### Check the gradientMidpoint parameter #####
    if (is.null(gradientMidpoint)) {
        gradientMidpoint <- 0.2
        memo <- paste0("gradientMidpoint parameter cannot be NULL. Setting the gradientMidpoint value to 0.")
        message(memo)
    }
    ## Check if gradientMidpoint is numeric
    if (!is.numeric(gradientMidpoint)) {
        memo <- paste0("gradientMidpoint variable not of the numeric class. Attempting to coerce.")
        gradientMidpoint <- as.numeric(gradientMidpoint)
        message(memo)
    }
    ## Check if the value for gradinetMidpoint is between 0 and 0.6
    if (gradientMidpoint > 0.6 | gradientMidpoint < 0.0) {
        memo <- paste0("gradientMidpoint parameter must be between 0 and 0.6. This value defines the midpoint of colors ",
                       "that are used to visualize regions of LOH. Attempting to divide by 100 (i.e. convert 20 to 0.2).")
        message(memo)
        gradientMidpoint <- gradientMidpoint/100
    }
    
    ##### Check the gradientColors parameter #####
    ## Check if it is a character vector
    if (is.null(gradientColors)) {
        memo <- paste0("gradientColors was set to NULL. Using default colors.")
        message(memo)
        gradientColors <- c("#ffffff", "#b2b2ff", "#000000")
    }
    if (!is.character(gradientColors)) {
        memo <- paste0("gradientColors variable not of the character class. Attempting to coerce.")
        gradientColors <- as.character(gradientColors)
        message(memo)
    }
    areColors <- function(x) {
        sapply(x, function(X) {
            tryCatch(is.matrix(col2rgb(X)), 
                     error = function(e) FALSE)
        })
    }
    if (any(areColors(gradientColors) == FALSE)) {
        ## Get the invalid color
        nonColor <- gradientColors[which(data.table(areColors(gradientColors))$V1==FALSE)]
        memo <- paste0("The ", nonColor, " designated in the gradientColors parameter is not a valid color. ",
                       "Making the gradient colors to depict regions of LOH to be: #ffffff, #b2b2ff, and #000000.")
        
        message(memo)
    }
    
    ##### Check the plotAType parameter #####
    if (is.null(plotAType) | !(plotAType %in% c("proportion", "frequency"))){
        memo <- paste0("The plotAType parameter cannot be NULL. Must be either 'proportion' or 'frequency' ", 
                       "Will plot the proportion of samples with LOH.")
        message(memo)
        plotAType <- "proportion"
    }
    
    ##### Check the plotALohCutoff parameter ##### 
    if (is.null(plotALohCutoff)) {
        plotALohCutoff <- 0.1
        memo <- paste0("plotALohCutoff parameter cannot be NULL. Setting the plotALohCutoff value to 0.1.")
        message(memo)
    }
    ## Check if plotALohCutoff is numeric
    if (!is.numeric(plotALohCutoff)) {
        memo <- paste0("plotALohCutoff variable not of the numeric class. Attempting to coerce.")
        plotALohCutoff <- as.numeric(plotALohCutoff)
        message(memo)
    }
    ##### Check the plotAColor parameter #####
    ## Check if it is a character vector
    if (is.null(plotAColor)) {
        memo <- paste0("plotAColor was set to NULL. Using default color.")
        message(memo)
        plotAColor <- c("#98F5FF")
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
        memo <- paste0("The ", nonColor, " designated in the gradientColors parameter is not a valid color. ",
                       "Making the gradient colors to depict regions of LOH to be: #98F5FF.")
        
        message(memo)
        plotAColor <- c("#98F5FF")
    }
    
    ##### Check the plotALayers and plotBLayers parameter #####
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
    
    ##### Check the sectionHeights parameter #####
    ## Check if it not NULL
    if (is.null(sectionHeights)) {
        sectionHeights <- c(0.25, 0.75)
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
        sectionHeights <- c(0.25, 0.75)
    }
    
    ## Check that there are 2 values in the variable
    if (length(sectionHeights)!=2) {
        memo <- paste0("2 values must be supplied to the sectionHeights parameter, which specifies the ",
                       "relative height of the loh frequency plot and the loh heatmap.")
        message(memo)
        sectionHeights <- c(0.25, 0.75)
    }
    
    ## Check that the values sum up to 1
    if (sum(sectionHeights)!=1) {
        memo <- paste0("sectionHeight values do not equal 1. Using default values.")
        message(memo)
        sectionHeights <- c(0.25, 0.75)
    }
    
    
    new("lohInputParameters", getHeterozygousCalls=getHeterozygousCalls, chromosomes=chromosomes, 
        sample=sample, BSgenome=BSgenome, step=step, windowSize=windowSize, 
        normalVAF=normalVAF, gradientMidpoint=gradientMidpoint, gradientColors=gradientColors,
        plotAType=plotAType, plotALohCutoff=plotALohCutoff, plotAColor=plotAColor, 
        plotALayers=plotALayers, plotBLayers=plotBLayers, sectionHeights=sectionHeights, 
        verbose=verbose)
}

#' Private Class lohData
#' 
#' An S4 class for the Data of the loh plot object
#' @name lohData-class
setClass("lohData",
         representation=representation(primaryData="data.table",
                                       windowData="data.table",
                                       windowCalcData="data.table",
                                       chrData="data.table", 
                                       lohFreqData="data.table"),
         validity = function(object){
             
         }
)

#' Constructor for the lohData class.
#' 
#' @name lohData
#' @rdname lohData-class
#' @param object Object of class VarScan 
lohData <- function(object, getHeterozygousCalls, chromosomes, samples, BSgenome, step, windowSize, 
                    normal,  plotALohCutoff, verbose) {
    
    ## Obtain LOH data for desired chromosomes and samples
    primaryData <- getLohData(object=object, verbose=verbose, getHeterozygousCalls=TRUE, 
                              germline=FALSE)
    
    ## Subset data to only the desired chromosomes to be plotted
    primaryData <- chrSubset(object=primaryData, chromosomes=chromosomes, 
                             verbose=verbose)
    
    ## Subset data to only the desired samples to be plotted
    primaryData <- sampleSubset(object=primaryData, samples=samples, 
                                verbose=verbose)
    
    ## Obtain chromosome boundaries from BSgenome object
    chrData <- annoGenomeCoord(object=primaryData, BSgenome=BSgenome, 
                               verbose=verbose)
    
    ## Produce data.table with window position data
    windowData <- getLohSlidingWindow(object=primaryData, step=step, 
                                      windowSize=windowSize, verbose=verbose)
    
    ## Perform loh calculations on each chromosome and sample within each window
    lohAbsDiff <- getLohCalculation(object=primaryData, 
                                    windowData=windowData, normal=normal, 
                                    verbose=verbose)
    
    ## Calculate avg loh for overlapping regions
    lohAbsDiffOverlap <- rbindlist(getLohStepCalculation(object=lohAbsDiff, 
                                                         step=step, verbose=verbose))
    
    ## Obtain LOH segmentation dataset
    lohSegmentation <- getLohSegmentation(object=lohAbsDiffOverlap, 
                                          verbose=verbose)
    
    ## Obtain LOH frequency/proportion dataset
    lohFreq <- getLohFreq(object=lohSegmentation, plotALohCutoff=plotALohCutoff,
                          chrData=chrData, verbose=verbose)
    
    ## Initialize the object
    new("lohData", primaryData=primaryData, windowData=rbindlist(windowData), 
        windowCalcData=lohAbsDiffOverlap, chrData=chrData, 
        lohFreqData=lohFreq)
}

#' Private Class lohSpecPlots
#' 
#' An S4 class for the plots of the lohSpec class
#' @name lohSpecPlots-class
#' @rdname lohSpecPlots-class
#' @slot PlotA gtable object for the loh spectrum
#' @slot PlotB gtable object for the loh frequency/proportion
#' @import methods
#' @importFrom gtable gtable
#' @noRd
setClass("lohSpecPlots", 
         representation=representation(PlotA="gtable",
                                       PlotB="gtable"),
         validity = function(object) {
             
         })

#' Constructor for the lohSpecPlots class
#' 
#' @name lohSpecPlots
#' @rdname lohSpecPlots-class
#' @param object Object of class lohData
#' @importFrom gtable gtable
#' @noRd
lohSpecPlots <- function(object, plotALohCutoff, plotAType, plotAColor, 
                         plotALayers, plotBLayers, gradientMidpoint, gradientColors, verbose) {
    ## Use the loh segmentation data to generate lohFreq plots
    lohFreqPlot <- buildLohFreq(object=object, plotALohCutoff=plotALohCutoff, 
                                plotAType=plotAType, plotAColor=plotAColor,
                                plotALayers=plotALayers, verbose=verbose)
    
    ## Use the lohData to generate lohSpec plots
    lohSpecPlot <- lohSpec_buildMainPlot(object=object, gradientMidpoint=gradientMidpoint,
                                         gradientColors=gradientColors, plotBLayers=plotBLayers, verbose=verbose)
    
    new("lohSpecPlots", PlotA=lohFreqPlot, PlotB=lohSpecPlot)
}

################################################################################
###################### Accessor function definitions ###########################

#' Helper function to get data from classes
#' 
#' @rdname getData-methods
#' @aliases getData
.getData_lohSpec <- function(object, name=NULL, index=NULL, ...) {
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
          signature="lohData",
          definition=.getData_lohSpec)

#' @rdname getData-methods
#' @aliases getData
setMethod(f="getData",
          signature="lohSpec",
          definition=.getData_lohSpec)

#' Helper function to extract grobs from objects
#'
#' @rdname getGrob-methods
#' @aliases getGrob
#' @noRd
.getGrob_lohSpec <- function(object, index=1, ...){
    if(index == 1){
        grob <- object@PlotA
    } else if(index == 2) {
        grob <- object@PlotB
    } else if (index == 3) {
        grob <- object@Grob
    } else {
        stop("Subscript out of bounds") 
    }
    return(grob)
}

#' @rdname getGrob-methods
#' @aliases getGrob
setMethod(f="getGrob",
          signature="lohSpecPlots",
          definition=.getGrob_lohSpec)

#' @rdname getGrob-methods
#' @aliases getGrob
setMethod(f="getGrob",
          signature="lohSpec",
          definition=.getGrob_lohSpec)


#' @rdname drawPlot-methods
#' @aliases drawPlot
#' @importFrom grid grid.draw
#' @importFrom grid grid.newpage
#' @exportMethod drawPlot
setMethod(
    f="drawPlot",
    signature="lohSpec",
    definition=function(object, ...){
        mainPlot <- getGrob(object, index=3)
        grid::grid.newpage()
        grid::grid.draw(mainPlot)
    }
)

################################################################################
###################### Method function definitions #############################

######################################################
##### Function to obtain chromosomes of interest #####
#' @rdname getLohData-methods
#' @aliases getLohData
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
              
              ## Check format of the chromosome column
              if (!all(grepl("^chr", object$chromosome))) {
                  memo <- paste0("Did not detect the prefix chr in the chromosome column",
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
                  chromosomes <- paste("chr", as.character(c(seq(1:22))), sep="")
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
#' @rdname getLohData-methods
#' @aliases getLohData
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

#####################################################
##### Function to get the chromosome boundaries #####
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
              
              ## Perform quality check on BSgenome object
              if (is.null(BSgenome)) {
                  memo <- paste("BSgenome object is not specified, whole chromosomes",
                                "will not be plotted, this is not recommended!")
                  warning(memo)
                  object$chromosome <- factor(object$chromosome, levels=gtools::mixedsort(unique(as.character(object$chromosome))))
                  return(object)
              } else if (is(BSgenome, "BSgenome")) {
                  if(verbose){
                      memo <- paste("BSgenome passed object validity checks")
                  }
              } else {
                  memo <- paste("class of the BSgenome object is", class(BSgenome),
                                "should either be of class BSgenome or NULL",
                                "setting this to param to NULL")
                  warning(memo)
                  BSgenome <- NULL
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
#' @param object of class lohData 
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
              
              ## Perform quality check on input variables
              
              ## Check that step and windowSize are numeric vectors with length of 1
              if (!is.numeric(windowSize)) {
                  memo <- paste("WindowSize input value is not a numeric vector, attempting to coerce...")
                  warning(memo)
              }
              if (!is.numeric(step)) {
                  memo <- paste("Step input value is not a numeric vector, attempting to coerce...")
                  warning(memo)
              }
              if (length(windowSize) > 1) {
                  memo <- paste("Use only 1 numeric value to specify window size.")
                  warning(memo)
                  stop()
              }
              if (length(step) > 1) {
                  memo <- paste("Use only 1 numeric value to specify step size.")
                  warning(memo)
                  stop()
              }
              if (step > windowSize) {
                  memo <- paste("Step value is greater than windowSize. Make sure that the step value is 
                                at most equal to the WindowSize. Changing step value to match the windowSize value.")
                  warning(memo)
                  step <- windowSize
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
#' @param object of class lohData 
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
#' @param object of class lohData
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

############################################################################
##### Function to create dataset for the loh frequency/proportion plot #####
#' @rdname getLohFreq-methods
#' @param object of class lohData
#' @param chrData of class data.table 
#' @aliases getLohFreq
setMethod(f="getLohFreq", 
          signature="data.table",
          definition=function(object, plotALohCutoff, chrData, verbose, ...){
              
              ## Print status message
              if (verbose) {
                  message("Determining proportion/frequency of samples with LOH in each region.")
              }
              
              x <- object[,c("chrom", "loc.start", 
                             "loc.end", "seg.mean", "ID")]
              colnames(x) <- c("chromosome", "start", "end", "segmean", "sample")
              
              ## Remove any NA values in the data
              if (any(is.na(x))) {
                  na_rows_removed <- nrow(x) - nrow(na.omit(x))
                  memo <- paste0("Removing", na_rows_removed, " rows containing NA values")
                  message(memo)
                  x <- na.omit(x)
              }
              
              ## Make sure windows are consistent, if not disjoin them
              tmp <- split(x, x$sample)
              tmp_vec <- tmp[[1]]$end
              if (any(!unlist(sapply(tmp, function(x) x[,"end"] %in% tmp_vec), use.names=F))) {
                  memo <- paste0("Did not detect identical genomic segments for all samples",
                                 " ...Performing disjoin operation")
                  message(memo) 
                  # here we split the DF up in an attempt to avoid complaints that lists are to large 
                  split_df <- split(x, f=x$chromosome)
                  x <- rbindlist(lapply(split_df, function(x) {
                      # Create the Granges object for the data
                      granges <- GenomicRanges::GRanges(seqnames=x$chromosome,
                                                        ranges=IRanges::IRanges(start=x$start, end=x$end),
                                                        "sample"=x$sample, "segmean"=x$segmean)
                      
                      # disjoin with grange, get a mapping of meta columns and expand it
                      disJoint <- GenomicRanges::disjoin(granges, with.revmap=TRUE)
                      revmap <- GenomicRanges::mcols(disJoint)$revmap
                      disJoint <- rep(disJoint, lengths(revmap))
                      
                      # exract the meta columns and map them back to the disJoint GRanges object
                      sample <- unlist(IRanges::extractList(GenomicRanges::mcols(granges)$sample, revmap))
                      segmean <- unlist(IRanges::extractList(GenomicRanges::mcols(granges)$segmean, revmap))
                      GenomicRanges::mcols(disJoint)$sample <- sample
                      GenomicRanges::mcols(disJoint)$segmean <- segmean
                      
                      # convert the GRanges Object back to a data table
                      disJoint <- as.data.table(disJoint)[,c("seqnames", "start", "end", "width",
                                                             "sample", "segmean")]
                      colnames(disJoint) <- c("chromosome", "start", "end", "width", "sample", "segmean")
                      return(disJoint)
                  }))
              }
              
              ## Calculate columns of observed LOH and observed samples in the 
              ## cohort for each segment
              gainFreq <- function(x){length(x[x>=plotALohCutoff])}
              gainFrequency <- aggregate(segmean~chromosome + start + end, 
                                         data=x, gainFreq)$segmean
              x <- aggregate(segmean~chromosome + start + end, data=x, length)
              colnames(x)[which(colnames(x) %in% "segmean")] <- "sampleFrequency"
              x$gainFrequency <- gainFrequency
              
              ## Calculate the proportion
              x$gainProportion <- as.numeric(x$gainFrequency)/length(samples)
              x <- data.table(x)
              return(x)
          })

#################################################
##### Function to create loh frequency plot #####
#' @rdname buildLohFreq-methods
#' @param object of class lohData
#' @aliases buildLohFreq
setMethod(f = "buildLohFreq",
          signature="lohData",
          definition=function(object, plotALohCutoff, plotAType, plotAColor, 
                              plotALayers, verbose, ...) {
              
              ## Print status message
              if (verbose) {
                  message("Building LOH frequency or proportion plot")
              }
              
              ## Perform quality checks on the input variables
              if (!is.numeric(plotALohCutoff)) {
                  memo <- paste("LOH cutoff value is not numeric, attempting to coerce...")
                  message(memo)
              }
              if (!grepl("^#(\\d|[a-f]){6,8}$", plotAColor, ignore.case=TRUE)){
                  memo <- paste("LOH frequency/proportion color is not a valid hexadecimal code.")
                  message(memo)
              }
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
              finalDf <- object@lohFreqData
              
              ## Sort the chromosomes
              chr <- gtools::mixedsort(as.character((unique(finalDf$chromosome))))
              sample <- gtools::mixedsort(as.character((unique(finalDf$sample))))
              finalDf$chromosome <- factor(finalDf$chromosome, levels=chr, labels=chr)
              finalDf$sample <- factor(finalDf$sample, levels=sample, labels=sample)
              
              ## Build the frequency/proportion plot
              ## Define parameters of the plot
              plotTheme <- theme(axis.ticks.x=element_blank(),
                                 axis.text.x=element_blank(),
                                 axis.ticks.y=element_blank(),
                                 panel.grid.major=element_blank(),
                                 panel.grid.minor=element_blank(), 
                                 legend.position="none")
              
              ## Define the facet
              facet <- facet_grid(.~chromosome, scales="free_x", space="fixed")
              
              ## Assign the x axis label
              xlabel <- xlab("Chromosome")
              
              ## Choose whether to plot aesthetics for proportion or frequency
              if(grepl("^PROP", plotAType, ignore.case=TRUE)){
                  ylabel <- ylab("Proportion of Loss of Heterozygosity")
                  ymax <- 1
                  finalDf$gain <- finalDf$gainProportion
              } else if(grepl("^FREQ", plotAType, ignore.case=TRUE)){
                  ylabel <- ylab("Frequency of Loss of Heterozygosity")
                  ymax <- max(as.numeric(as.character(x$sampleFrequency)), na.rm=TRUE)
                  finalDf$gain <- finalDf$gainFrequency
              } else {
                  memo <- paste0("did not recognize plotAType ", plotAType,
                                 ", please specify one of \"proportion\" or \"frequency\"")
                  stop(memo)
              }
              
              ## Initiate the plot
              finalDf$gain <- as.numeric(finalDf$gain)
              finalDf$start <- as.numeric(finalDf$start)
              finalDf$end <- as.numeric(finalDf$end)
              p1 <- ggplot(data=finalDf, mapping=aes_string(xmin='start', 
                                                            xmax='end',
                                                            ymin=0,
                                                            ymax='gain')) +
                  geom_rect(fill=plotAColor) + scale_x_continuous(expand=c(0,0)) + 
                  scale_y_continuous(expand=c(0,0), limits=c(0,1))
              
              p1 <- p1 + geom_hline(aes(yintercept=0), linetype="dotted")
              
              # build the plot
              p1 <- p1 + ylabel + xlabel + facet + theme_bw() + plotTheme + plotALayers
              print(p1)
              
              ## Convert to grob
              lohFreqGrob <- ggplotGrob(p1)
              return(lohFreqGrob)
              
          })

################################################
##### Function to generate lohSpec heatmap #####
#' @rdname lohSpec_buildMainPlot-methods
#' @param object of class lohData
#' @aliases lohSpec_buildMainPlot
setMethod(f = "lohSpec_buildMainPlot",
          signature="lohData",
          definition=function(object, gradientMidpoint, gradientColors, 
                              plotBLayers, verbose, ...) {
              
              ## Print status message
              if (verbose) {
                  message("Building main LOH spectrum plot")
              }
              
              ## Perform quality checks on the input variables
              if (!is.numeric(gradientMidpoint)) {
                  memo <- paste("Gradient midpoint value is not numeric, attempting to coerce...")
                  message(memo)
              }
              sapply(gradientColors, function(x) {
                  if (!is.character(x)) {
                      memo <- paste("Gradient colors for LOH spectrum figure is not a character vector, 
                                    attempting to coerce...")
                      message(memo)
                  }
                  hexColor <- (grepl("^#(\\d|[a-f]){6,8}$", x))
                  if (hexColor == FALSE) {
                      memo <- paste("Specified colors in the gradient are not hexadecimal.")
                      message(memo)
                  }
              })
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
              
              x <- object@windowCalcData
              x <- x[loh_step_avg >= 0]
              x$start <- as.numeric(x$start)
              x$stop <- as.numeric(x$stop)
              x$loh_step_avg <- as.numeric(x$loh_step_avg)
              
              ## Set the order of the chromosomes
              chr <- gtools::mixedsort(as.character((unique(x$chromosome))))
              sample <- gtools::mixedsort(as.character((unique(x$sample))))
              x$chromosome_f <- factor(x$chromosome, levels=chr)
              x$sample <- factor(x$sample, levels=sample, labels=sample)
              
              # Define the main plot
              data <- geom_rect(data=x, aes_string(xmin='start',
                                                   xmax='stop',
                                                   ymin=-1,
                                                   ymax=1, fill='loh_step_avg'))
              
              # Define additional plot parameters
              facet <- facet_grid(sample ~ chromosome_f, scales="free_x", space="fixed")
              
              x_scale <- scale_x_continuous(expand = c(0, 0))
              y_scale <- scale_y_continuous(expand = c(0,0))
              
              lab_x <- xlab("Chromosome")
              lab_y <- ylab("Sample")
              
              # Define plot aesthetics
              BWscheme <- theme_bw()
              plotTheme <- theme(axis.ticks.x=element_blank(),
                                 axis.text.x=element_blank(),
                                 axis.ticks.y=element_blank(),
                                 axis.text.y=element_blank(),
                                 panel.grid.major=element_blank(),
                                 panel.grid.minor=element_blank())
              
              # plot an additional layer if specified
              if(!is.null(plotBLayers))
              {
                  plotLayer <- plotBLayers
              } else {
                  plotLayer <- geom_blank()
              }
              
              LOHgradient <- scale_fill_gradient2(midpoint = gradientMidpoint,
                                                  guide="colourbar",
                                                  high=gradientColors[3],
                                                  mid=gradientColors[2],
                                                  low=gradientColors[1],
                                                  space='Lab')
              
              # Build the plot
              tmp <- data.frame(x=0, y=0)
              p1 <- ggplot(data=tmp, aes(y=0)) +
                  data + facet + x_scale + y_scale + 
                  lab_x + lab_y + BWscheme + LOHgradient + plotTheme + plotLayer
              print(p1)
              
              ## Convert to grob
              lohSpecGrob <- ggplotGrob(p1)
              return(lohSpecGrob)
              })

#########################################################
##### Function to arrange lohSpec and lohFreq plots #####
#' @rdname arrangeLohPlots-methods
#' @param object of class lohData
#' @aliases arrangeLohPlots
setMethod(f="arrangeLohPlots",
          signature="lohSpecPlots",
          definition=function(object, sectionHeights, verbose, ...) {
              
              ## Print status message
              if (verbose) {
                  message("Combining LOH frequency/proportion and LOH spectrum plot")
              }
              
              
              ## Perform quality checkes on input parameters
              if (!is.numeric(sectionHeights)) {
                  memo <- paste("Values specified for the section heights are 
                                not numeric, attempting to coerce...")
                  message(memo)
              }
              if (length(sectionHeights) != 2) {
                  memo <- paste("Heights for both LOH figures are not specified. The sectionHegihts
                                variable should be a numeric vector of legth 2 specifying the heights of each of the 
                                2 LOH figures.")
                  message(memo)
                  stop()
              } 
              
              ## Grab the data we need
              plotA <- object@PlotA
              plotB <- object@PlotB
              
              ## obtain the meax width for relevant plots
              plotList <- list(plotA, plotB)
              plotList <- plotList[lapply(plotList, length) > 0]
              plotWidths <- lapply(plotList, function(x) x$widths)
              maxWidth <- do.call(grid::unit.pmax, plotWidths)
              
              ## Set the widths for all plots
              for (i in 1:length(plotList)) {
                  plotList[[i]]$widths <- maxWidth
              }
              
              ## Set section heights based upon the number of sections
              defaultPlotHeights <- c(0.25, 0.75)
              
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
              
              ## Arrange the final plot
              finalPlot <- do.call(gridExtra::arrangeGrob, c(plotList, list(ncol=1, heights=sectionHeights)))
              plot(finalPlot)
              return(finalPlot)
          })