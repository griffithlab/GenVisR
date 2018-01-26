################################################################################
##################### Public/Private Class Definitions #########################

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Public Class !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#' Class cnLoh
#' 
#' An S4 class for the cn, somatic loh, and germline loh plots
#' @name cnLoh
#' @rdname cnLoh-class
#' @slot cnData data.table object for cn plot
#' @slot cnPlot gtable object for the cn plot
#' @slot somaticLohData data.table object for the somatic loh plot
#' @slot somaticLohPlot gtable object for the somatic loh plot
#' @slot germlineLohData data.table object for the germline loh plot
#' @slot germlineLohData gtable object for the germline loh plot
#' @exportClass cnLoh
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

cnLoh <- function(cnInput, lohInput, samples, chromosomes, BSgenome, windowSize,
                  step, normal, plotAColor, plotALayers, plotBAlpha,
                  somaticLohCutoff, plotBTumorColor, plotBNormalColor, plotBLayers,
                  plotCLimits, plotCLowColor, plotCHighColor, plotCLayers, 
                  sectionHeights, verbose) {
    
    ## Obtain cn, somatic loh, and germline loh datasets to plot
    cnLohDataset <- cnLohData(cnInput=cnInput, lohInput=lohInput, samples=samples, chromosomes=chromosomes,
                              BSgenome=BSgenome, windowSize=windowSize,
                              step=step, normal=normal, verbose=verbose)
    
    ## Generate the cn, somatic LOH, and germline LOH plots
    plots <- cnLohPlots(object=cnLohDataset, plotAColor=plotAColor, plotALayers=plotALayers, 
                        plotBAlpha=plotBAlpha, somaticLohCutoff=somaticLohCutoff, plotBTumorColor=plotBTumorColor, 
                        plotBNormalColor=plotBNormalColor, plotBLayers=plotBLayers, 
                        plotCLimits=plotCLimits, 
                        plotCLowColor=plotCLowColor, plotCHighColor=plotCHighColor, 
                        plotCLayers=plotCLayers, verbose=verbose)
    
    ## Arrange all of the plots together
    Grob <- arrangeCnLohPlots(object=plots, sectionHeights=sectionHeights, verbose=verbose)
    
    ## Initialize the object
    new("cnLoh", cnData=getData(cnLohDataset, index=1), cnPlot=getGrob(plots, index=1),
        somaticLohData=getData(cnLohDataset, index=2), somaticLohPlot=getGrob(plots, index=2),
        germlineLohData=getData(cnLohDataset, index=3), germlineLohPlot=getGrob(plots, index=3),
        Grob=Grob)
}

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Private Classes !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

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
    
    ## Obtain copy number segmentation data
    cnSegmentation <- getCnSegmentation(object=cnData, verbose=verbose)  
    
    ## Obtain chromosome boundaries from BSgenome object
    chrData <- annoGenomeCoord(object=cnData, BSgenome=BSgenome, verbose=verbose)
    
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
.getGrob_combinedCnLoh <- function(object, index=1, ...){
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
#' @rdname cnLohData-methods
#' @aliases cnLohData
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
              
              # perform quality checks on the chromosome parameter arguments
              
              # check for character vector
              if(!is.character(chromosomes)){
                  memo <- paste("Input to chromosomes should be a character vector,
                                specifying which chromosomes to plot, 
                                attempting to coerce...")
                  warning(memo)
                  chromosomes <- as.character(chromosomes)
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
#' @rdname cnLohData-methods
#' @aliases cnLohData
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
              
              ## If samples is null, we don't want to do anything and just return the object
              if (is.null(samples)) {
                  return(object)
              }
              
              ## Perform quality checkes on the sample parameter arguments
              if (!is.character(samples)) {
                  memo <- paste("Input to samples should be a character vector, 
                                attempting to coerce...")
                  warning(memo)
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
#' @aliases getCnSegmentation
#' @param object of class data.table
setMethod(f="getCnSegmentation",
          signature="data.table",
          definition=function(object, ...) {
              
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

#####################################################
##### Function to get the chromosome boundaries #####
#' @rdname annoGenomeCoord-methods
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
              
              ## Perform quality checkes on the input parameters
              if (!is.logical(normal)) {
                  memo <- ("Input to specify normal VAF should be a boolean (T/F). True if 
                           user wants to use normal VAF from varscan to identify tumor/normal LOH difference. 
                           Flase if user wants to use 0.5 to identify tumor/normal LOH difference.")
                  message(memo)
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
#' @aliases buildCnPlot
#' @param object of class data.table
setMethod(f="buildCnPlot",
          signature="cnLohData",
          definition=function(object, plotAColor, plotALayers, ...) {
              
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
#' @aliases buildSomaticLohPlot
#' @param object of class data.table
setMethod(f="buildSomaticLohPlot",
          signature="cnLohData",
          definition=function(object, somaticLohCutoff, plotBAlpha, plotBTumorColor, plotBNormalColor,
                              plotBLayers, ...) {
              
              ## Print status message
              if (verbose) {
                  message("Building somatic loh plot")
              }
              
              ## Perform quality checks on the input variables
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
#' @aliases buildGermlineLohPlot
#' @param object of class data.table
setMethod(f="buildGermlineLohPlot",
          signature="cnLohData",
          definition=function(object, plotCLimits, plotCLowColor,
                              plotCHighColor, plotCLayers, verbose=verbose){
              
              ## Print status message
              if (verbose) {
                  message("Building germline loh plot")
              }
              
              ## Perform quality checks on the input variables
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
#' @param object of class cnLohData
#' @aliases arrangeCnLohPlots
setMethod(f="arrangeCnLohPlots",
          signature="cnLohPlots",
          definition=function(object, sectionHeights, verbose, ...) {
              
              ## Print status message
              if (verbose) {
                  message("Combining cn, somatic loh, and germline loh plots")
              }
              
              ## Perform quality checkes on input parameters
              if (!is.numeric(sectionHeights)) {
                  memo <- paste("Values specified for the section heights are 
                                not numeric, attempting to coerce...")
                  message(memo)
              }
              if (length(sectionHeights) != 3) {
                  memo <- paste("Heights for cn/loh figures are not specified. The sectionHegihts
                                variable should be a numeric vector of legth 3 specifying the heights of each of the 
                                3 figures (cn, somatic loh, and germline loh).")
                  message(memo)
                  stop()
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
              
              ## Set section heights based upon the number of sections
              defaultPlotHeights <- c(0.3, 0.4, 0.3)
              
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