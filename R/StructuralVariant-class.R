################################################################################
##################### Public/Private Class Definitions #########################

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Public Class !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

#' Class StructuralVariant
#' 
#' An S4 class for the Structural Variant plot object
#' @rdname StructuralVariant-class
#' @name StructuralVariant
#' @slot primaryData data.table object storing primarydata used for plotting
#' @slot geneData data.table object storing annotated gene files
#' @slot Grob gtable object for the structural variant plot
#' @import methods
#' @importFrom gtable gtable
#' @importFrom data.table data.table
#' @exportClass StructuralVariant
setClass("StructuralVariant",
         representation=representation(svData="data.table",
                                       geneData="data.table",
                                       svPlots="svPlots"),
         validity=function(object) {
             
         })

#' Constuctor for the Structural Variant class
#' 
#' @name svData
#' @rdname svData-class
#' @param object Object of class VCF
#' @rdname StructuralVariant-class
#' @name StructuralVariant
#' @param object OBject of class VCF
#' @param BSgenome Object of class BSgenome to extract genome wide chromosome
#' coordinates
#' @param filter Boolean specifying if SV calls that did not pass should be removed
#' @param svType Character vector specifying which structural variant types to annotate/visualize
#' @param svOrder Character vector specifying the deleterious order of sv types (most to least deleterious)
#' @param maxSvSize Numeric specifying the maximum size of SV events (DEL/DUP/INV only)
#' @param sample Character vector specifying which samples to annotate/visualize
#' @param chromosomes Character vector specifying chromosomes to annotate/visualize
#' @param ensembl Object of class Mart to use in biomaRt query
#' @param attributes Character vector specifying which attributes to retrieve from biomaRt query
#' @param filters Character vector specifying which filters to use in biomaRt query
#' @param annotate Boolean specifying if the user wants to obtain mutated gene counts and annotate SV events
#' @param geneAnnotationFlank Integer specifying the size of the flanks of each SV event
#' to include in the annotation step
#' @param plotSpecificGene Character vector specifying which genes to plot
#' @param plotTraGenes Boolean specifying if TRA genes should be plotted
#' @param plotOtherGenes Boolean specifying if non-TRA genes should be plotted
#' @param chrGap Integer specifying the size of the gap between the 1st and 2nd chromosome
#' @param genome Character vector specifying which genome to use to obtain chromosome bands. 
#' Serves as input into the getCytobands function of karyoploteR. 
#' @param cytobandColor Character vector specifying what to color the chromosome bands
#' @param sampleColor Character vector specifying colors to plot for each sample
#' @param plotALAyers List of ggplot2 layers to be passed to translocation plot
#' @param plotBLayers List of ggplot2 layers to be passed to chromosome plot
#' @param plotCLayers List of ggplot2 layers to be passed to non-translocation plot
#' @param outputDir Character value for directory to output SV visualizations
#' @param plotWidth Integer for width of SV visualizations
#' @param plotHeight Integer for height of SV visualizations
#' @param verbose Boolean specifying if status messages should be reported
#' @export
StructuralVariant <- function(input, BSgenome=NULL, filterSvCalls=TRUE, svType=NULL,
                              svOrder=c("TRA", "BND", "DEL", "DUP", "INV", "INS"),
                              maxSvSize=NULL, sample=NULL, chromosomes=NULL, 
                              ensembl=ensembl, attributes=attributes, filters=filters,
                              annotate=TRUE, geneAnnotationFlank=10000, 
                              plotSV=plotSV, plotSpecificGene=FALSE, plotTraGenes=FALSE, 
                              plotOtherGenes=FALSE, chrGap=5000000,
                              genome="hg19", cytobandColor=c("White", "Grey"), 
                              sampleColor=NULL, verbose=FALSE, plotALayers=NULL, 
                              plotBLayers=NULL, plotCLayers=NULL, sectionHeights=c(0.4, 0.1, 0.5)) {

    ## Check the input parameters
    inputParameters <- checkSvInputParameters(object=object, BSgenome=BSgenome, filterSvCalls=filterSvCalls, 
                                              svType=svType, svOrder=svOrder,
                                              maxSvSize=maxSvSize, sample=sample, chromosomes=chromosomes, 
                                              ensembl=ensembl, attributes=attributes, 
                                              filters=filters, chrGap=chrGap, annotate=annotate, 
                                              geneAnnotationFlank=geneAnnotationFlank, genome=genome, 
                                              plotSV=plotSV, plotSpecificGene=plotSpecificGene, 
                                              plotTraGenes=plotTraGenes, plotOtherGenes=plotOtherGenes, 
                                              cytobandColor=cytobandColor, 
                                              plotALayers=plotALayers, plotBLayers=plotBLayers,
                                              plotCLayers=plotCLayers, sectionHeights=sectionHeights, 
                                              sampleColor=sampleColor, verbose=verbose)
    
    ## Calculate all data for the plots
    svDataset <- svData(object=input, BSgenome=inputParameters@BSgenome, 
                        filterSvCalls=inputParameters@filterSvCalls, 
                        svType=inputParameters@svType, svOrder=inputParameters@svOrder,
                        maxSvSize=inputParameters@maxSvSize, sample=inputParameters@sample, 
                        chromosomes=inputParameters@chromosomes, 
                        ensembl=inputParameters@ensembl, attributes=inputParameters@attributes, 
                        filters=inputParameters@filters, chrGap=inputParameters@chrGap, 
                        annotate=inputParameters@annotate, 
                        geneAnnotationFlank=inputParameters@geneAnnotationFlank, 
                        genome=inputParameters@genome, verbose=inputParameters@verbose)

    ## Create the plots from svData
    structuralVariantPlots <- svPlots(object=svDataset, plotSV=inputParameters@plotSV, 
                                      plotSpecificGene=inputParameters@plotSpecificGene, 
                                      plotTraGenes=inputParameters@plotTraGenes, 
                                      plotOtherGenes=inputParameters@plotOtherGenes, 
                                      cytobandColor=inputParameters@cytobandColor, 
                                      plotALayers=inputParameters@plotALayers, 
                                      plotBLayers=inputParameters@plotBLayers,
                                      plotCLayers=inputParameters@plotCLayers, 
                                      sectionHeights=inputParameters@sectionHeights, 
                                      sample=inputParameters@sample, sampleColor=inputParameters@sampleColor, 
                                      verbose=inputParameters@verbose)
    
    ## Intialize the object
    new("StructuralVariant", svData=getData(object=svDataset, name="primaryData"), 
        geneData=getData(object=svDataset, name="geneData"),
        svPlots=structuralVariantPlots)
}

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Private Classes !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

#' Private Class svInputParameters
#' 
#' An S4 class to check input parameters of the StructuralVariant function
#' @name svInputParameters-class
#' @noRd
setClass("svInputParameters",
         representation=representation(BSgenome="BSgenome", filterSvCalls="logical", svType="character", 
                                       svOrder="character", maxSvSize="numeric", 
                                       sample="character", chromosomes="character", 
                                       ensembl="Mart", attributes="character", 
                                       filters="character", chrGap="numeric", 
                                       annotate="logical", geneAnnotationFlank="numeric", 
                                       genome="character", plotSV="logical", 
                                       plotSpecificGene="character", plotTraGenes="logical", 
                                       plotOtherGenes="logical", cytobandColor="character", 
                                       plotALayers="list", plotBLayers="list",
                                       plotCLayers="list", sectionHeights="numeric", 
                                       sampleColor="character", verbose="logical"), 
         validity=function(object){
             
         })

#' Constructor for the svInputParameters class
#' 
#' @name svInputParameters
#' @rdname svInputParameters-class
#' @importFrom karyoploteR getCytobands
#' @noRd
checkSvInputParameters <- function(object, BSgenome, filterSvCalls, svType, svOrder, maxSvSize, 
                                   sample, chromosomes, ensembl, attributes, 
                                   filters, chrGap, annotate, geneAnnotationFlank, 
                                   genome, plotSV, plotSpecificGene, plotTraGenes, 
                                   plotOtherGenes, cytobandColor, plotALayers, plotBLayers,
                                   plotCLayers, sectionHeights, sampleColor, verbose) {
    
    ##### Check verbose parameter #####
    ## Check to see if verbose is a booelean
    if (!is.logical(verbose) | is.null(verbose)) {
        memo <- paste0("The verbose parameter is not a boolean (T/F). Coercing verbose to be FALSE...")
        message(memo)
        verbose <- FALSE
    }
    
    ##### Check BSgenome parameter #####
    ## Check to see if BSgenome is a BSgenome
    if (is.null(BSgenome)) {
        memo <- paste("BSgenome object is not specified, whole chromosomes",
                      "will not be plotted, this is not recommended!")
        warning(memo)
    } else if (is(BSgenome, "BSgenome")) {
        memo <- paste("BSgenome passed object validity checks")
        message(memo)
    } else {
        memo <- paste("class of the BSgenome object is", class(BSgenome),
                      "should either be of class BSgenome or NULL",
                      "setting this to param to NULL")
        warning(memo)
        BSgenome <- NULL
    }
    
    ##### Check filterSvCalls parameter ##### 
    ## Check to see if filterSvCalls is a boolean
    if (is.null(filterSvCalls) | !is.logical(filterSvCalls)) {
        memo <- paste0("The filterSvCalls parameter is not a boolean (T/F). Coercing filterSvCalls to be FALSE...")
        message(memo)
        filterSvCalls <- FALSE
    }
    
    ##### Check svType parameter #####
    ## Check if it null
    if (is.null(svType)) {
        svType <- as.character(object@vcfObject@svType$svtype)
        memo <- paste0("svType variable cannot be NULL. Using all of the sv types ",
                       "available in the sv dataset")
        message(memo)
    }
    ## check if svType is a character vector
    if (!is.character(svType)) {
        memo <- paste0("svType variable not of the character class. Attempting to coerce.")
        svType <- as.character(svType)
        message(memo)
    }
    
    ##### Check svOrder parameter #####
    ## Check if svOrder is NULL
    if (is.null(svOrder)) {
        svOrder <- svType
        memo <- paste0("svOrder variable cannot be NULL. Setting the deleterious order ",
                       "of sv events to that designated in the svType: ", 
                       paste(svType, collapse=", "))
        message(memo)
    }
    if (!is.character(svOrder)) {
        memo <- paste0("svOrder variable not of the character class. Attempting to coerce.")
        svOrder <- as.character(svOrder)
        message(memo)
    }
    
    ##### Check maxSvSize parameter #####
    ## Check if maxSvSize is NULL
    if (is.null(maxSvSize)) {
        maxSvSize <- 0
        memo <- paste0("maxSvSize parameter cannot be NULL. Setting the maxSvSize value to 0.")
        message(memo)
    }
    ## Check if maxSvSize is numeric
    if (!is.numeric(maxSvSize)) {
        memo <- paste0("maxSvSize variable not of the numeric class. Attempting to coerce.")
        maxSvSize <- as.numeric(maxSvSize)
        message(memo)
    }
    
    ##### Check sample parameter #####
    ## Check is sample is NULL
    if (is.null(sample)) {
        sample <- unique(object@vcfObject@sample$sample)
        memo <- paste0("Sample parameter cannot be NULL. All samples will be plotted.")
    }
    ## Check if sample is a character vector
    if (!is.character(sample)) {
        memo <- paste0("sample variable not of the character class. Attempting to coerce.")
        sample <- as.character(sample)
        message(memo)
    }
    
    ##### Check chromosomes parameter #####
    ## Check is chromosomes is NULL
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
    ##### Check annotate parameter #####
    ## Check if annotate is a boolean
    if (!is.logical(annotate) | is.null(annotate)) {
        memo <- paste0("The annotate parameter is not a boolean (T/F). Coercing annotate to be FALSE...")
        message(memo)
        annotate <- FALSE
    }
    ##### Check ensembl parameter #####
    ## Check if ensembl is of class "mart"
    if (is.null(ensembl)) {
        memo <- paste("ensembl object cannot be NULL if SV annotations are desired, ",
                      "in which case, the ensembl object must be of class Mart...")
        if (!annotate) {
            warning(memo)
        }
        if (annotate) {
            stop(memo)
        }
    } 
    else if (is(ensembl, "Mart")) {
        memo <- paste("ensembl passed object validity checks")
        message(memo)
    } 
    else {
        memo <- paste("class of the ensembl object is", class(ensembl),
                      "should either be of class ensembl or NULL",
                      "setting this to param to NULL and will not perform ",
                      "sv annotations.")
        warning(memo)
        ensembl <- NULL
        annotate <- FALSE
        
    }
    
    
    ##### Check if attributes and filters are valid #####
    if (annotate) {
        ## If ensembl is not NULL, check if these are character vectors
        if (!is.null(ensembl) & (!is.character(attributes))) {
            if (is.null(attributes)) {
                memo <- paste0("If annotations are desired, the attributes parameter ",
                               "cannot be NULL. These values are used to specify the ouput ",
                               "for biomaRt annotations.")
                stop(memo)
            }
            memo <- paste0("attributes variable not of the character class. Attempting to coerce.")
            attributes <- as.character(attributes)
            message(memo)
        }
        if (!is.null(ensembl) & (!is.character(filters))) {
            if (is.null(filters)) {
                memo <- paste0("If annotations are desired, the filters parameter ",
                               "cannot be NULL. These values are used to specify the input ",
                               "for biomaRt annotations.")
                stop(memo)
            }
            memo <- paste0("filters variable not of the character class. Attempting to coerce.")
            filters <- as.character(filters)
            message(memo)
        }
        
        ## If ensembl is not NULL, check that these are valid inputs for getBM
        if (!is.null(ensembl)) {
            temp <- data.table(listAttributes(mart=ensembl))
            if (!all(attributes %in% temp$name)) {
                `%nin%` = Negate(`%in%`)
                discrepantAttributes <- attributes[which(attributes %nin% temp$name)]
                memo <- paste0("The following attributes: ", paste(discrepantAttributes, collapse="|"),
                               " are not valid inputs for the designated ensembl database. Please run ",
                               "biomaRt::listAttributes(ensembl) to get valid attributes.")
                stop(memo)
            }
            temp <- data.table(listFilters(mart=ensembl))
            if (!all(filters %in% temp$name)) {
                `%nin%` = Negate(`%in%`)
                discrepantFilters <- filters[which(filters %nin% temp$name)]
                memo <- paste0("The following filters: ", paste(discrepantFilters, collapse="|"),
                               " are not valid inputs for the designated ensembl database. Please run ",
                               "biomaRt::listFilters(ensembl) to get valid filters")
                stop(memo)
            }
        }
        
        ## If ensembl is NULL, set these variables to NULL
        if (is.null(ensembl)) {
            attributes <- NULL
            filters <- NULL
        }   
    }
    
    ##### Check chrGap parameter #####
    ## Check if chrGap is NULL
    if (is.null(chrGap)) {
        chrGap <- 5000000
        memo <- paste0("chrGap variable cannot be NULL. Using the default value of 5,000,000.")
    }
    ## Check if chrGap is numeric
    if (!is.numeric(chrGap)) {
        memo <- paste0("chrGap variable not of the numeric class. Attempting to coerce.")
        chrGap <- as.numeric(chrGap)
        message(memo)
    }
    
    ##### Check geneAnnotationFlank parameter #####
    ## Check if geneAnnotationFlank is NULL
    if (is.null(geneAnnotationFlank)) {
        geneAnnotationFlank <- 10000
        memo <- paste0("geneAnnotationFlank variable cannot be NULL. Setting the variable's ",
                       "value to 10,000 base pairs.")
        message(memo)
    }
    ## Check if geneAnnotationFlank is numeric and is greater than 0
    if (!is.numeric(geneAnnotationFlank)) {
        memo <- paste0("geneAnnotationFlank variable not of the numeric class. Attempting to coerce.")
        geneAnnotationFlank <- as.numeric(geneAnnotationFlank)
        message(memo)
        if (geneAnnotationFlank < 0) {
            memo <- paste0("geneAnnotationFlank cannot be a negative number. Changing ",
                           "geneAnnotationFlank to be 0.")
            message(memo)
            geneAnnotationFlank <- 0
        }
    }
    
    ##### Check genome parameter #####
    ## Check if genome is not NULL
    if (is.null(genome)) {
        memo <- paste0("The genome variable cannot be NULL. Valid options are those used by ",
                       "the karyoploteR package (e.g. hg19, mm10, etc...)")
        stop(memo)
        
    }
    ## Check if genome is of length 1
    if (length(genome) > 1) {
        memo <- paste0("The genome variable must be of length 1. Using the first value.")
        genome <- genome[1]
        message(memo)
    }
    ## Check if genome is a character
    if (!is.character(genome) & !is.null(genome)) {
        memo <- paste0("genome variable not of the character class. Attempting to coerce.")
        genome <- as.character(genome)
        message(memo)
    }
    ## Check if the genome exists in KaryoploteR
    temp <- suppressMessages(getCytobands(genome))
    if (nrow(as.data.table(temp@seqnames))==0) {
        memo <- paste0("The inputted genome is not available in the karyoploteR package, ",
                       "which is used to generate the cytoband positions... Please submit a request to the ",
                       "karyoploteR github page: https://github.com/bernatgel/karyoploteR")
        stop(memo)
    }
    
    ##### Check plotSV parameter #####
    ## Check to see if filter is a boolean
    if (!is.logical(plotSV)) {
        memo <- paste0("plotSV parameter is not a boolean (T/F). Coercing plotSpecificGene to be FALSE...")
        plotSV <- FALSE
        message(memo)
    }
    
    ##### Check plotSpecificGene parameter #####
    ## Check if plotSpecificGene is NULL
    if (is.null(plotSpecificGene)) {
        plotSpecificGene <- ""
    }
    ## Check to see if plotSpecificGene is a booelean
    if (!is.character(plotSpecificGene)) {
        memo <- paste0("The plotSpecificGene variable not of the character class. Attempting to coerce...")
        message(memo)
        plotSpecificGene <- as.character(plotSpecificGene)
    }
    ##### Check plotTraGenes parameter #####
    ## Check to see if plotTraGenes is a booelean
    if (!is.logical(plotTraGenes)) {
        memo <- paste0("The plotTraGenes parameter is not a boolean (T/F). Coercing plotTraGenes to be FALSE...")
        message(memo)
        plotTraGenes <- FALSE
    }
    ##### Check plotOtherGenes parameter #####
    ## Check to see if plotOtherGenes is a booelean
    if (!is.logical(plotOtherGenes)) {
        memo <- paste0("The plotOtherGenes parameter is not a boolean (T/F). Coercing plotOtherGenes to be FALSE...")
        message(memo)
        plotOtherGenes <- FALSE
    }
    
    ##### Check cytobandColor parameter #####
    ## Check if it is a character vector
    if (is.null(cytobandColor)) {
        memo <- paste0("cytobandColor was set to NULL. Setting the colors to Dark grey and light grey.")
        message(memo)
        cytobandColor <- c("Dark Grey", "Light Grey")
    }
    if (!is.character(cytobandColor)) {
        memo <- paste0("cytobandColor variable not of the character class. Attempting to coerce.")
        cytobandColor <- as.character(cytobandColor)
        message(memo)
    }
    ## Check if desired colors are valid 
    areColors <- function(x) {
        sapply(x, function(X) {
            tryCatch(is.matrix(col2rgb(X)), 
                     error = function(e) FALSE)
        })
    }
    if (any(areColors(cytobandColor) == FALSE)) {
        ## Get the invalid color
        nonColor <- cytobandColor[which(data.table(areColors(cytobandColor))$V1==FALSE)]
        memo <- paste0("The ", nonColor, " designated in the cytobandColor parameter is not a valid color. ",
                       "Making the cytoband colors dark grey and light grey.")
    }
    
    ##### Check sampleColor parameter #####
    ## Check if sampleColor is NULL
    if (is.null(sampleColor)) {
        ## Set the sampleColor to be of the same length as the number of samples
        sampleColor <- rainbow(length(sample))
        memo <- paste0("sampleColor parameter cannot be NULL...attempting to generate distinctive colors.")
        message(memo)
    }
    if (!is.null(sampleColor)) {
        ## Check if it is a character vector
        if (!is.character(sampleColor)) {
            memo <- paste0("sampleColor variable not of the character class. Attempting to coerce.")
            sampleColor <- as.character(sampleColor)
            message(memo)
        }
        
        ## Check if desired colors are valid 
        areColors <- function(x) {
            sapply(x, function(X) {
                tryCatch(is.matrix(col2rgb(X)), 
                         error = function(e) FALSE)
            })
        }
        if (any(areColors(sampleColor) == FALSE)) {
            ## Get the invalid color
            nonColor <- sampleColor[which(data.table(areColors(sampleColor))$V1==FALSE)]
            memo <- paste0("The ", nonColor, " designated in the sampleColor parameter is not a valid color. ",
                           "Making the cytoband colors dark grey and light grey.")
        }
        
        ## If sampleColor is not NULL, check if it's length is the same as 
        ## the desired number of samples 
        if (length(sampleColor) != length(sample)) {
            memo <- paste0("The number of colors for each sample (designated with the sampleColor variable) ",
                           "does not equal the number of samples in the sv dataset (n=", length(sample), ") .")
            stop(memo)
        }
    }
    
    ##### Check plotALayers, plotBLayers, and plotC Layers parameter #####
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
    ##### Check sectionHeights parameter #####
    ## Check if it not NULL
    if (is.null(sectionHeights)) {
        sectionHeights <- c(0.4, 0.1, 0.5)
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
        sectionHeights <- c(0.4, 0.1, 0.5)
    }
    
    ## Check that there are 3 values in the variable
    if (length(sectionHeights)!=3) {
        memo <- paste0("3 values must be supplied to the sectionHeights parameter, which specifies the ",
                       "relative height of the plot for translocations, the chromosomes, and non-translocations ",
                       "respectively.")
        message(memo)
        sectionHeights <- c(0.4, 0.1, 0.5)
    }
    
    ## Check that the values sum up to 1
    if (sum(sectionHeights)!=1) {
        memo <- paste0("sectionHeight values do not equal 1. Using default values.")
        message(memo)
        sectionHeights <- c(0.4, 0.1, 0.5)
    }
    
    new("svInputParameters", BSgenome=BSgenome, filterSvCalls=filterSvCalls, 
        svType=svType, svOrder=svOrder,
        maxSvSize=maxSvSize, sample=sample, chromosomes=chromosomes, 
        ensembl=ensembl, attributes=attributes, filters=filters, chrGap=chrGap, annotate=annotate, 
        geneAnnotationFlank=geneAnnotationFlank, genome=genome, plotSV=plotSV, plotSpecificGene=plotSpecificGene, 
        plotTraGenes=plotTraGenes, plotOtherGenes=plotOtherGenes, cytobandColor=cytobandColor, 
        plotALayers=plotALayers, plotBLayers=plotBLayers,
        plotCLayers=plotCLayers, sectionHeights=sectionHeights, 
        sampleColor=sampleColor, verbose=verbose)
}

#' Private Class svData
#' 
#' An S4 class for the data of the sv plot object
#' @name svData-class
#' @noRd
setClass("svData", 
         representation=representation(primaryData="data.table",
                                       geneData="data.table",
                                       chrData="data.table",
                                       svWindow="data.table",
                                       cytobands="data.table"),
         validity=function(object){
             
         })

#' Constructor for the svData class
#' 
#' @name svData
#' @rdname svData-class
#' @name StructuralVariant
#' @importFrom data.table data.table
#' @noRd
svData <- function(object, BSgenome, filterSvCalls, svType, svOrder, maxSvSize, sample, 
                   chromosomes, ensembl, attributes, filters, annotate, geneAnnotationFlank, chrGap, genome,
                   verbose) {
    ## Subset data to only passed sv calls
    primaryData <- getVcfData(object=object, filterSvCalls=filterSvCalls, maxSvSize=maxSvSize, 
                              svType=svType, verbose=verbose)
    
    ## Subset data to only the chromosomes desired to be plotted
    primaryData <- chrSubsetSv(object=primaryData, chromosomes=chromosomes, verbose=verbose)
    
    ## Subset data to only the samples desired to be plotted
    primaryData <- sampleSubset(object=primaryData, samples=sample, verbose=verbose)
    
    ## Obtain chromosome boundaries from BSgenome object
    chrData <- annoGenomeCoordSv(object=primaryData, BSgenome=BSgenome, 
                               verbose=verbose)
    
    ## Annotate the sv calls
    primaryData <- annotateSV(object=primaryData, ensembl=ensembl, attributes=attributes, filters=filters,
                              annotate=annotate, geneAnnotationFlank=geneAnnotationFlank, chromosomes=chromosomes, 
                              verbose=verbose)
    
    ## Get the proportion of samples that have each mutated gene
    geneData <- countGenes(object=primaryData, annotate=annotate, svOrder=svOrder, verbose=verbose)
    
    ## Get the cytoband data
    chrCytobands <- svCytobands(object=primaryData, genome=genome, chrData=chrData, verbose=verbose)
    
    ## Adjust the primaryData to account for centromeres
    adjustedPrimaryData <- adjustCentromeres(object=primaryData, chrCytobands=chrCytobands, verbose=verbose)
    
    ## Get the new positions for SV calls and cytobands
    svWindow <- getStructuralVariantWindow(object=adjustedPrimaryData, chromosomes=chromosomes, chrCytobands=chrCytobands, chrData=chrData, 
                                           chrGap=chrGap, verbose=verbose)

    ## Initialize the object
    new("svData", primaryData=primaryData, geneData=geneData, chrData=chrData, svWindow=svWindow, cytobands=chrCytobands)
}

#' Private Class svPlots
#' 
#' An S4 class for the of the svData class
#' @name svPlots-class
#' @rdname svPlots-class
#' @slot Plots list of gtables for each chr combo
#' @import methods
#' @importFrom gtable gtable
#' @noRd
setClass("svPlots", 
         representation=representation(plots="list"),
         validity = function(object) {
             
         })

#' Constructor for the svPlots class
#' 
#' @name svPlots
#' @rdname svPlots-class
#' @param object Object of class svData
#' @importFrom gtable gtable
#' @noRd
svPlots <- function(object, plotSV, plotSpecificGene, plotTraGenes, plotOtherGenes, cytobandColor, 
                    plotALayers, plotBLayers, plotCLayers, sectionHeights, 
                    sample, sampleColor, verbose, ...) {
    
    ## Create the gtable for the plots 
    svGtables <- buildSvPlot(object=object, plotSV=plotSV, plotSpecificGene=plotSpecificGene, 
                             plotTraGenes=plotTraGenes, plotOtherGenes=plotOtherGenes, cytobandColor=cytobandColor, 
                             plotALayers=plotALayers, plotBLayers=plotBLayers,
                             plotCLayers=plotCLayers, sectionHeights=sectionHeights, sample=sample,
                             sampleColor=sampleColor, verbose=verbose)
    
    ## Initialize the object
    new("svPlots", plots=svGtables)
   
}

################################################################################
###################### Accessor function definitions ###########################

#' Helper function to get data from classes
#' 
#' @rdname getData-methods
#' @aliases getData
.getData_structuralVariants <- function(object, name=NULL, index=NULL, ...) {
    if(is.null(name) & is.null(index)){
        memo <- paste("Both name and index are NULL, one must be specified!")
        stop(memo)
    }
    
    if(is.null(index)){
        index <- 0
    } else {
        if(index > 5){
            memo <- paste("index out of bounds")
            stop(memo)
        }
    }
    
    if(is.null(name)){
        name <- "noMatch"
    } else {
        slotAvailableName <- c("primaryData", "geneData", "chrData", "svWindow", "cytobands")
        if(!(name %in% slotAvailableName)){
            memo <- paste("slot name not found, specify one of:", toString(slotAvailableName))
            stop(memo)
        }
    }
    
    if(name == "primaryData" | index == 1){
        data <- object@primaryData
    }
    if(name == "geneData" | index == 2){
        data <- object@geneData
    }
    
    return(data)
}

#' @rdname getData-methods
#' @aliases getData
setMethod(f="getData",
          signature="svData",
          definition=.getData_structuralVariants)

#' @rdname drawPlot-methods
#' @aliases drawPlot
#' @importFrom grid grid.draw
#' @importFrom grid grid.newpage
#' @export
setMethod(f="drawPlot",
          signature="StructuralVariant",
          definition=function(object, chr1=NULL, chr2=NULL, ...) {
              ## Get the list of gtabls
              object <- object@svPlots@plots
              
              ## Get the chr combo
              order <- paste(gtools::mixedsort(c(chr1, chr2)), collapse="_")
              
              ## See if the desired chr combo can be found in the plots
              num <- which(names(object) == order)
              if (length(num) == 0) {
                  memo <- paste0("The plot for the chromosome combination: ",
                                 chr1, "_", chr2, " could not be found. Make sure to append the chr name with ",
                                 dQuote("chr"), " rather than just using the chromosome number (chr1 instead of 1).",
                                 "The possible combinations that could be used are: ", 
                                 paste(names(object), collapse=", "))
                  stop(memo)
              }
              finalPlot <- object[[num]]
              grid::grid.newpage()
              grid::grid.draw(finalPlot)
          })

################################################################################
#################### Method function definitions ###############################

######################################################
##### Function to obtain chromosomes of interest #####
#' @rdname svData-methods
#' @aliases svData
#' @param object Object of class data.table
#' @param chromosomes character vector of chromosomes to retain
#' @param verbose Boolean for status updates
#' @return data.table object with calculated mutation distances
#' @importFrom data.table data.table
#' @importFrom data.table rbindlist
#' @noRd
setMethod(f="chrSubsetSv",
          signature="data.table",
          definition=function(object, chromosomes, verbose, ...){
              # print status message
              if(verbose){
                  memo <- paste("Performing chromosome subsets")
                  message(memo)
              }
              
              # if chromosomes is null we dont want to do anything just return all autosomes
              if(is.null(chromosomes)){
                  chromosomes <- "autosomes"
              }              
             
              ## Check format of the chromosome1 column
              if (!all(grepl("^chr", object$chromosome))) {
                  if (verbose) {
                      memo <- paste0("Did not detect the prefix chr in the chromosome1 column ",
                                     "of x... adding prefix")
                      message (memo)
                      object$chromosome <- paste("chr", object$chromosome, sep="")
                  }
              } else if (all(grepl("^chr", object$chromosome))) {
                  if (verbose) {
                      memo <- paste0("Detected chr in the chromosome1 column of x...",
                                 "proceeding")
                      message(memo)
                  }
              } else {
                      memo <- paste0("Detected unknown or mixed prefixes in the chromosome1",
                                 " colum of object... should either be chr or non (i.e.) chr1 or 1")
                      message(memo)
              }
              
              ## Check format of the chromosome2 column
              if (!all(grepl("^chr", object$chromosome2))) {
                  if (verbose) {
                      memo <- paste0("Did not detect the prefix chr in the chromosome2 column",
                                     "of x... adding prefix")
                      message (memo)
                      object$chromosome2 <- paste("chr", object$chromosome2, sep="")
                  }
              } else if (all(grepl("^chr", object$chromosome2))) {
                  if (verbose) { 
                      memo <- paste0("Detected chr in the chromosome2 column of x...",
                                     "proceeding")
                      message(memo)    
                  }
              } else {
                  memo <- paste0("Detected unknown or mixed prefixes in the chromosome2",
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
                  chromosomes <- unique(c(object$chromosome, object$chromosome2))
                  chromosomes <- chromosomes[-grep("GL", chromosomes)]
                  chromosomes <- chromosomes[-grep("MT", chromosomes)]
              }
              
              # check for specified chromosomes not in the original input
              missingChr <- chromosomes[!(chromosomes %in% unique(object$chromosome) & 
                                              chromosomes %in% unique(object$chromosome2))]
              if(length(missingChr) != 0){
                  memo <- paste("The following chromosomes were designated to be kept but were not found:",
                                toString(missingChr), "\nValid chromosomes are", toString(unique(object$chromosome)))
                  warning(memo)
              }
              
              # perform the subset - remove GL and MT chromosomes
              object <- object[!grepl("GL", object$chromosome)]
              object <- object[!grepl("MT", object$chromosome)]
              object <- object[!grepl("GL", object$chromosome2)]
              object <- object[!grepl("MT", object$chromosome2)]
              
              ## Remove rows that have nothing to do with the desired chromosomes
              ## Keep DEL/DUP/INV/INS events that occur on other chromosomes that 
              ## have translocations with the desired chromosomes
              allStructuralVariants <- object 
              ## Get the chromosome combination
              chr_combo <- data.table(paste(allStructuralVariants$chromosome, 
                                 allStructuralVariants$chromosome2,
                                 sep="_"))
              otherChromosomes <- apply(chr_combo, 1, function(x, chromosomes){
                  chr <- strsplit(x, split="_")[[1]]
                  otherChr <- data.table(chr[-which(chr %in% chromosomes)])
                  if (length(otherChr)!=0) {
                      return(otherChr)
                  }
              }, chromosomes=chromosomes)
              otherChromosomes <- unique(rbindlist(otherChromosomes))
              otherChromosomes <- paste(otherChromosomes$V1, otherChromosomes$V1, sep="_")
              
              ## Write message to say that only non-TRA events are shown
              if (length(otherChromosomes)==0) {
                  memo <- paste0("No translocation events detected in the dataset.")
                  if (verbose){
                      message(memo)
                  }
              }
               
              object$chr_combo <- paste(object$chromosome, object$chromosome2, sep="_")
              object <- object[object$chromosome %in% chromosomes | object$chromosome2 %in% chromosomes | 
                                   object$chr_combo %in% otherChromosomes,]
              object$chromosome <- factor(object$chromosome)
              object$chromosome2 <- factor(object$chromosome2)
              
              ## Remove rows that are duplciated in the ID column
              object <- object[!duplicated(object$ID)]
              object <- object[,-c("chr_combo")]
              
              # check that the object has a size after subsets
              if(nrow(object) < 1){
                  memo <- paste("no entries left to plot after chromosome subsets")
                  stop(memo)
              }
              
              return(object)
          })

#####################################################
##### Function to get the chromosome boundaries #####
#' @rdname annoGenomeCoordSv-methods
#' @aliases annoGenomeCoordSv
#' @param object Object of class data.table
#' @param BSgenome Object of class BSgenome, used for extracting chromosome boundaries
#' @param verbose Boolean for status updates
#' @return Data.table with chr and start/stop positions
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom data.table data.table
#' @importFrom data.table rbindlist
#' @importFrom gtools mixedsort
#' @noRd
setMethod(f="annoGenomeCoordSv", 
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
              
              ## Get all of the chromosomes that are in the SV dataset
              all_chr <- unique(c(as.character(object$chromosome), as.character(object$chromosome2)))
              
              ## Check that chromosomes between BSgenome and original input match
              chrMismatch <- all_chr[!all_chr %in% genomeCoord$chromosome]
              if (length(chrMismatch) >= 1) {
                  memo <- paste("The following chromosomes do not match the supplied BSgenome object",
                                toString(chrMismatch))
                  warning(memo)
                  
                  ## Test if the chr mismatch is fixed by appending chr to chromosomes
                  all_chr <- paste("chr", all_chr, sep="")
                  chrMismatch_appendChr <- all_chr[!all_chr %in% genomeCoord$chromosome]
                  if(chrMismatch_appendChr < length(chrMismatch)){
                      memo <- paste("appending \"chr\" to chromosomes in attempt to fix mismatch with the BSgenome")
                      warning(memo)
                      object$chromosome <- paste0("chr", object$chromosome)
                  }
              }
              
              ## Check to see if any chromosomes in the original input dataset lack genomic coordiantes
              if (any(!all_chr %in% unique(genomeCoord$chromosome))) {
                  missingGenomeCoord <- unique(object$chromosome)
                  missingGenomeCoord <- missingGenomeCoord[!missingGenomeCoord %in% unique(genomeCoord_a$chromosome)]
                  memo <- paste("The following chromosomes are missing genomic coordinates", toString(missingGenomeCoord),
                                "Full genomic coordinates will not be plotted for these chromosomes")
                  warning(memo)
              }
              
              ## Filter the genomeCoord objext to only inlcude chromosomes in the input data
              genomeCoord <- genomeCoord[genomeCoord$chromosome %in% all_chr,]
              
              return(genomeCoord)
              
          })

##################################################
##### Function to obtain samples of interest #####
#' @rdname svData-methods
#' @aliases svData
#' @param object Object of class data.table
#' @param samples character vector of samples to retain
#' @param verbose Boolean for status updates
#' @return data.table object with calculated mutation distances
#' @importFrom data.table data.table
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
                                kept but were not found:", toString(missingSamp), 
                                "\nValid csamples are", 
                                toString(unique(object$sample)))
                  warning(memo)
              } 
              
              ## Perform the subset
              object <- object[object$sample %in% samples]
              object$sample <- factor(object$sample)
              
              ## Remove rows that are duplciated in the ID column
              object <- object[!duplicated(object$ID)]
              
              ## Check that the object has a size after subsets
              if(nrow(object) < 1){
                  memo <- paste("no entries left to plot after chromosome subsets")
                  stop(memo)
              }
              
              return(object)
          })

##########################################
##### Function to annotate SV events #####
#' @rdname svData-methods
#' @aliases svData
#' @param object Object of class data.table
#' @param samples character vector of samples to retain
#' @param verbose Boolean for status updates
#' @return data.table object with calculated mutation distances
#' @importFrom data.table data.table
#' @noRd
setMethod(f="annotateSV",
          signature="data.table",
          definition=function(object, annotate, ensembl, attributes, filters, geneAnnotationFlank, 
                              chromosomes, verbose, ...){
    
              if (annotate == TRUE) {
                  # print status message
                  if(verbose){
                      memo <- paste("Annotating sv positions")
                      message(memo)
                  }
                  
                  ## Define the chromosomes to annotate
                  all_chr <- unique(c(as.character(object$chromosome), as.character(object$chromosome2)))
                  
                  ## Go through each row of the primaryData dataset and run through biomaRt
                  annotatedDf <- data.table::rbindlist(apply(object, 1, function(x, object, ensembl, attributes, filters, 
                                                                                 geneAnnotationFlank, verbose){
                      x <- data.table(t(x))
                      if (verbose) {
                          num <- which(object$chromosome==x$chromosome & 
                                           as.numeric(object$position)==as.numeric(as.character(x$position)) &
                                           object$chromosome2==x$chromosome2 & 
                                           as.numeric(object$position2)==as.numeric(as.character(x$position2)))
                          print(paste("Annotating call: ", num, "/", nrow(object), sep=""))
                      }
                      ## Get the breakpoint information
                      chr1 <- as.character(x$chromosome[1])
                      chr2 <- as.character(x$chromosome2[1])
                      chr1 <- gsub(pattern="chr", replacement="", x=chr1)
                      chr2 <- gsub(pattern="chr", replacement="", x=chr2)
                      pos1 <- as.numeric(x$position[1])
                      leftPos1 <- pos1
                      rightPos1 <- pos1
                      pos2 <- as.numeric(x$position2[1])
                      leftPos2 <- pos2
                      rightPos2 <- pos2
                      if (geneAnnotationFlank > 0) {
                          leftPos1 <- pos1 - geneAnnotationFlank
                          rightPos1 <- pos1 + geneAnnotationFlank
                          leftPos2 <- pos2 - geneAnnotationFlank
                          rightPos2 <- pos2 + geneAnnotationFlank 
                      }
                      
                      ## Annotate the first breakpoint (TRA and BND)
                      if (x$svtype == "BND" | x$svtype == "TRA") {
                          gene1 <- getBM(attributes=attributes, filters=filters, 
                                         values=list("chr"=chr1, "start"=leftPos1, "end"=rightPos1), mart=ensembl)$hgnc_symbol
                          gene2 <- getBM(attributes=attributes, filters=filters, 
                                         values=list("chr"=chr2, "start"=leftPos2, "end"=rightPos2), mart=ensembl)$hgnc_symbol
                          
                          genes <- paste(c(gene1, gene2), collapse="|")
                      }
                      if (x$svtype == "DEL" | x$svtype == "DUP" | x$svtype == "INV" | x$svtype == "INS") {
                          genes <- as.character(getBM(attributes=attributes, filters=filters, 
                                                      values=list("chr"=chr1, "start"=leftPos1, "end"=rightPos2), 
                                                      mart=ensembl)$hgnc_symbol)
                          genes <- paste(genes, collapse="|")
                      }
                      
                      ## Substitute "pseudogene" and "" for "No Gene"
                      genes <- gsub(pattern="Pseudogene", replacement="No Gene", x=genes)
                      if (genes=="") {
                          genes <- "No Gene"
                      }

                      ## Append genes to the dataset
                      x$genes <- genes
                      
                      return(x)
                  },
                  object=object, ensembl=ensembl, attributes=attributes, filters=filters, 
                  geneAnnotationFlank=geneAnnotationFlank, verbose=verbose))
                
                  ## Get the columns of interest
                  cols <- c("chromosome", "position", "chromosome2", "position2", "direction", 
                            "svtype", "total_read_support", "sample", "ID", 
                            "tumorSample", "genes")
                  annotatedDf <- annotatedDf[,which(colnames(annotatedDf) %in% cols),with=FALSE]
              }
              if (annotate == FALSE) {
                  annotatedDf <- object
                  annotatedDf$genes <- ""
              }
              return(annotatedDf)
          })

#################################################################################
##### Function to count the number/proportion of samples with mutated genes #####
#' @rdname svData-methods
#' @aliases svData
#' @param object Object of class data.table
#' @param samples character vector of samples to retain
#' @param verbose Boolean for status updates
#' @return data.table object with calculated mutation distances
#' @importFrom data.table as.data.table
#' @noRd
setMethod(f="countGenes", 
          signature="data.table",
          definition=function(object, annotate, svOrder, verbose, ...) {
              
              if (annotate) {
                  ## Print status message
                  if (verbose) {
                      message("Calculating proportion of samples with mutations in each gene.")
                  }
                  
                  ## Get the list of mutated genes
                  object <- object[-which(object$genes == "No Gene")]
                  genes <- data.table(unique(unlist(strsplit(object$genes, split="|", fixed=TRUE))))
                  
                  ## Get the total number of samples 
                  sample_num <- length(unique(object$sample))
                  
                  ## Check if the svOrders are in the sv dataset
                  if (!all(svOrder %in% object@vcfObject@svType$svtype)) {
                      svOrder <- svOrder[which(svOrder %in% object@vcfObject@svType$svtype)]
                      if (length(svOrder) == 0) {
                          memo <- paste0("Structural variant types in the svOrder variable are not ",
                                         "found in the SV dataset. Assigning the following order based on the ", 
                                         "mutations found in the SV dataset: ",
                                         paste(object@vcfObject@svType$svtype, collapse=" > "))
                          message(memo)
                          svOrder <- object@vcfObject@svType$svtype
                      }
                  }
                  
                  ## Go through the list of genes and see how many times/samples it is mutated
                  final_df <- data.table::rbindlist(apply(genes, 1, function(x, object, svOrder, sample_num) {
                      ## Get the rows with the genes
                      mt <- object[unique(grep(x, object$genes)),
                                   c("svtype", "sample","total_read_support")]
                      mt$gene <- x
                      mt$total_sample_num <- sample_num

                      ## Split the rows by sample and get top sv type
                      samples <- split(mt, mt$sample)
                      mt <- data.table::rbindlist(lapply(samples, function(y, svOrder) {
                          ## Set the order
                          setDT(y)[,y := factor(svtype, levels=svOrder)]
                          y <- y[order(svtype, -total_read_support),]
                          
                          ## Get the top row
                          final <- y[1,c(1:5)]
                          
                          ## Reorder the columns
                          final <- final[,c("gene", "sample", "svtype", 
                                            "total_read_support")]
                          return(final)
                          
                      }, 
                      svOrder=svOrder))
                      
                      ## Count the number of samples with each svtype
                      mutated_samples <- paste(mt$sample, collapse="|")
                      svtype <- paste(mt$svtype, collapse="|")
                      trs <- paste(mt$total_read_support, collapse="|")
                      proportion <- length(unique(mt$sample))/sample_num
                      
                      ## Make the final dataset
                      final <- data.table::data.table(gene=x, sample=mutated_samples,
                                          svtype=svtype, total_read_support=trs,
                                          proportion=proportion, total_sample_num=sample_num)

                      return(final)
                      
                  }, 
                  object=object, svOrder=svOrder, sample_num=sample_num))
                  
                  ## Order the final_df dataset by proportion
                  final_df <- final_df[order(-proportion, sample, svtype),]
              }
              if (annotate==FALSE) {
                  final_df <- data.table::data.table(gene="", sample="", svtype="",
                                         total_read_support="", proportion="",
                                         total_sample_num="")
              }
              return(final_df)
          })

################################################
##### Function to create cytobands dataset #####
#' @rdname svData-methods
#' @aliases svData
#' @param object Object of class data.table
#' @param samples character vector of samples to retain
#' @param verbose Boolean for status updates
#' @return data.table object with calculated mutation distances
#' @importFrom data.table data.table
#' @importFrom data.table as.data.table
#' @importFrom karyoploteR getCytobands
#' @importFrom karyoploteR getCytobandColors
#' @noRd
setMethod(f="svCytobands",
          signature="data.table",
          definition=function(object, genome, chrData, verbose) {
              
              ## Print status message
              if (verbose) {
                  message("Subsetting cytoband dataset.")
              }
              
              cytoband <- data.table::as.data.table(suppressMessages(getCytobands(genome=genome, use.cache = TRUE)))
              colnames(cytoband) <- c("chromosome", "start", "end", "width", "strand", "band", "stain")
              
              ## Subset the dataset by the chromosome
              cytoband <- cytoband[which(chromosome %in% chrData$chromosome)]
              
              ## Rename the cytoband 
              cytoband$band <- paste(cytoband$chromosome, cytoband$band, sep="")
              cytoband$band <- gsub("chr", "", cytoband$band)
              
              ## Get the centromere positions
              temp <- split(cytoband, f=as.character(cytoband$chromosome))
              finalCytoband <- data.table::rbindlist(lapply(temp, function(x, chrData){
                  ## Get the chromosome length
                  chrLength <- as.numeric(chrData[chromosome==x$chromosome[1],"end", with=FALSE])
                  ## Get the size of what the centromere should be
                  centromere <- chrLength*0.05
                  pBands <- x[grep("p", x$band),]
                  ## Get the start/stop position of the centromere
                  centromereStart <- pBands$end[nrow(pBands)]
                  centromereEnd <- centromereStart + centromere
                  centromereData <- data.table(chromosome=as.character(x$chromosome[1]), start=centromereStart, end=centromereEnd,
                                               width=0, strand="*", band="centromere", stain="black")
                  ## Get the new positions for the q bands
                  qBands <- x[grep("q", x$band),]
                  qBands$start <- qBands$start + centromere
                  qBands$end <- qBands$end + centromere
                  
                  final <- rbind(pBands, centromereData, qBands)
                  return(final)
              }, chrData=chrData))
              return(finalCytoband)
          })

###########################################################
##### Function to adjust sv positions for centromeres #####
#' @rdname svData-methods
#' @aliases svData
#' @param object Obejct of class data.table
#' @param chrCytobands Data.table with cytoband and centromere data
#' @param verbose Boolean for status updates
#' @return data.table object with adjusted sv positions
#' @importFrom data.table data.table
#' @noRd
setMethod(f="adjustCentromeres", 
          signature="data.table",
          definition=function(object, chrCytobands, verbose) {
              
              ## Print status message
              if (verbose){
                  message("Adjusting positions of structural variants to account for centromere to aid in visualization.")
              }
              
              ## Convert column values
              object$chromosome <- as.character(object$chromosome)
              object$position <- as.numeric(object$position)
              object$chromosome2 <- as.character(object$chromosome2)
              object$position2 <- as.numeric(object$position2)

              ## Adjust the sv positions for the centromeres 
              all_chr <- unique(c(object$chromosome, object$chromosome2))
              for (i in 1:length(all_chr)) {
                  chr <- as.character(all_chr[i])
                  
                  ## Get the centromere positions
                  centromereStart <- chrCytobands[chromosome==chr & band=="centromere", start]
                  centromereEnd <- chrCytobands[chromosome==chr & band=="centromere", end]
                  centromereLength <- centromereEnd-centromereStart
                  
                  ## For all genomic coordinates that are downstream of the centromere, add the centromere length
                  num <- which(object$chromosome==chr & object$position >= centromereStart)
                  if (length(num) > 0) {
                      object$position[num] <- object$position[num] + centromereLength
                  }
                  num <- which(object$chromosome2==chr & object$position2 >= centromereStart)
                  if (length(num) > 0) {
                      object$position2[num] <- object$position2[num] + centromereLength
                  }
              }
              return(object)
          })

#######################################################
##### Function to change positions of chromosomes #####
#' @rdname svData-methods
#' @aliases svData
#' @param object Object of class data.table
#' @param samples character vector of samples to retain
#' @param verbose Boolean for status updates
#' @return data.table object with calculated mutation distances
#' @importFrom data.table data.table
#' @importFrom ggforce geom_bezier
#' @noRd
setMethod(f="getStructuralVariantWindow",
          signature="data.table",
          definition=function(object, chromosomes, chrData, chrCytobands, chrGap, verbose){
              
              ## Print status message
              if (verbose) {
                  message("Adjusting chromosome boundaries for visualization of structural variants.")
              }
              
              ## Split the data by the sample
              sampleData <- object
              sampleData$chr_combo <- paste(sampleData$chromosome, sampleData$chromosome2, sep="_")

              ## Get the chromosomes of interest
              if (!is.null(chromosomes)) {
                  coi <- data.table(chromosomes)
              }
              if (is.null(chromosomes)) {
                  coi <- data.table(chromosomes=unique(chrData$chromosome))
              }
              
              ## For each of the COI-chr combination, generate a dataset to plot
              #chr <- coi$chromosomes[1]
              finalDf <- data.table::rbindlist(apply(coi, 1, function(chr, sampleData, chrCytobands, chrGap){
                  ## Get rows in cohort that have sv events involving COI
                  dataset <- sampleData[chromosome==chr | chromosome2 == chr,]
                  
                  ## Get the rows in the dataset that have sv events 2ndary to COI
                  otherChr1 <- unique(dataset[!which(chromosome %in% chr)]$chromosome)
                  otherChr2 <- unique(dataset[!which(chromosome2 %in% chr)]$chromosome2)
                  otherChr <- unique(c(otherChr1, otherChr2))
                  otherChrCombo <- paste(otherChr, otherChr, sep="_")
                  otherDataset <- sampleData[chr_combo %in% otherChrCombo]
                  
                  ## Combine COI TRA with 2ndary SV events
                  dataset <- rbind(dataset, otherDataset)
                  
                  ## Get the mate chromosome for TRA events
                  otherChr <- unique(c(as.character(dataset$chromosome), as.character(dataset$chromosome2)))
                  otherChr <- as.data.table(otherChr[-which(otherChr %in% chr)])
                  
                  ## For the individual COI, go through each of the chr combination if there are other chromosomes to plot
                  if (nrow(otherChr) > 0) {
                      #mateChr <- otherChr$V1[12]
                      temp <- data.table::rbindlist(apply(otherChr, 1, function(mateChr, chr, dataset, sampleData, chrCytobands, chrGap)
                          {
                          ## Get the chromosome combination to see rows which rows of the COI dataset have the correct
                          ## matching of the chr and mate chr
                          combo <- c(paste(chr, chr, sep="_"), paste(chr, mateChr, sep="_"), 
                                     paste(mateChr, chr, sep="_"), paste(mateChr, mateChr, sep="_"))
                          final <- dataset[which(chr_combo %in% combo), c("chromosome", "position", "chromosome2", "position2", "direction",
                                                                          "svtype", "total_read_support", "sample", "genes")]
                          
                          ## Get the chr boundaries for the two chromosomes
                          chrDataTemp <- chrCytobands[which(chromosome %in% c(mateChr, chr))]
                          
                          ## Figure out the order to plot the chromosomes
                          chrOrder <- gtools::mixedsort(unique(as.character(chrDataTemp$chromosome)))
                          
                          ## Order the chrDataTemp dataset by chrOrder
                          chrDataTemp$chromosome <- factor(chrDataTemp$chromosome, levels=chrOrder)
                          chrDataTemp <- chrDataTemp[order(chrDataTemp$chromosome)]
                          
                          ## Get the "end" for each of the chromosomes
                          chrLength <- data.table::rbindlist(lapply(split(chrDataTemp, f=chrDataTemp$chromosome), function(x) {
                              chr <- x$chromosome[1]
                              end <- max(x$end)
                              final <- data.table(chromosome=chr, chrLength=end)
                              return(final)
                          }))
                          
                          ## Add the cytoband/centromere data to the "final" dataset
                          temp <- data.table(chromosome=chrDataTemp$chromosome, position=chrDataTemp$start,
                                             chromosome2=chrDataTemp$chromosome, position2=chrDataTemp$end,
                                             direction="cytoband", svtype=chrDataTemp$band, total_read_support="",
                                             sample="", genes="")
                          
                          final <- rbind(final, temp)
                          
                          ## Add the length of the first chromosome to all of the sv calls on the second chromosome
                          num <- which(final$chromosome == chrLength$chromosome[2])
                          final$position[num] <- final$position[num] + chrGap + chrLength$chrLength[1]
                          num <- which(final$chromosome2 == chrLength$chromosome[2])
                          final$position2[num] <- final$position2[num] + chrGap + chrLength$chrLength[1]
                          
                          ## Get the midpoints for the dataset
                          final$midpoint <- (as.numeric(as.character(final$position))+
                                                 as.numeric(as.character(final$position2)))/2
                          
                          ## Get the chr-combo ID 
                          final$chr_combo <- paste(chrLength$chromosome[1], chrLength$chromosome[2], sep="_")
                          
                          return(final)
                          
                      }, 
                      chr=chr, dataset=dataset, sampleData=sampleData, chrCytobands=chrCytobands, chrGap=chrGap))   
                  }
                  if (nrow(otherChr) == 0) {
                      ## Get the positions for the sv events 
                      dataset$midpoint <- (dataset$position + dataset$position2)/2 
                      temp <- dataset[,c("chromosome", "position", "chromosome2", "position2",
                                         "direction", "svtype", "total_read_support",
                                         "sample", "genes", "midpoint", "chr_combo")]
                      
                      ## Add in the cytoband information
                      cytoTemp <- chrCytobands[chromosome==chr]
                      cyto <- data.table(chromosome=cytoTemp$chromosome, 
                                         position=cytoTemp$start,
                                         chromosome2=cytoTemp$chromosome,
                                         position2=cytoTemp$end,
                                         direction="cytoband", svtype=cytoTemp$band, 
                                         total_read_support="", sample="", genes="",
                                         midpoint=(cytoTemp$start+cytoTemp$end)/2,
                                         chr_combo=paste(cytoTemp$chromosome, cytoTemp$chromosome, sep="_"))
                      
                      ## Combine all of the data
                      temp <- rbind(temp, cyto)
                  }
                  
                  return(temp)
              }, 
              sampleData=sampleData, chrCytobands=chrCytobands, chrGap=chrGap))
              
              ## Get the genes for each SV call
              finalDf$gene <- as.character(finalDf$genes)
              finalDf$gene <- gsub(pattern="No Gene\\|", replacement="", finalDf$gene)
              finalDf$gene <- gsub(pattern="\\|No Gene", replacement="", finalDf$gene)
              finalDf$gene <- gsub(pattern="No Gene", replacement="", finalDf$gene)
              
              return(finalDf)
          })

################################################################
##### Function to generate SV plots for SVs on the same chr#####
#' @rdname svPlots-methods
#' @aliases svPlots
#' @param object object of class svData
#' @return list object of gtables for each chr_combo
#' @importFrom data.table data.table
#' @importFrom gtable gtable
#' @noRd
setMethod(f="buildSvPlot", 
          signature="svData", 
          definition=function(object, plotSV, plotSpecificGene, plotTraGenes, plotOtherGenes, 
                              cytobandColor, sample, sampleColor, plotALayers, 
                              plotBLayers, plotCLayers, sectionHeights,  verbose) {
              
              if (plotSV == FALSE) {
                  ## Print status message
                  if (verbose) {
                      message("plotSV was set to FALSE, so no SV plots will be generated.")
                  }
                  svPlots <- list()
                  return(svPlots)
              }
              
              if (plotSV == TRUE) {
                  ## Print status message
                  if (verbose) {
                      message("Generating SV plots")
                  }
                  
                  ## Get the svWindow
                  svWindow <- object@svWindow
                  
                  ## Convert BND notation to readable format
                  svWindow$direction <- gsub("N\\[P\\[", "3' to 5'", svWindow$direction)
                  svWindow$direction <- gsub("N]P]", "3' to 3'", svWindow$direction)
                  svWindow$direction <- gsub("]P]N", "5' to 3'", svWindow$direction)
                  svWindow$direction <- gsub("\\[P\\[N", "5' to 5'", svWindow$direction)
                  svWindow$svtype <- gsub("BND", "TRA", svWindow$svtype)
                  
                  ## Assign colors for samples 
                  names(sampleColor) <- sample
                  
                  ## Split the sv window by chr_combo
                  window <- split(svWindow, svWindow$chr_combo)
                  
                  ## Go through each window dataset and generate a plot
                  svPlots <- suppressWarnings(lapply(window, function(dataset, plotSpecificGene, 
                                                                      cytobandColor, sectionHeights, 
                                                                      sampleColor, 
                                                                      plotTraGenes, plotOtherGenes,
                                                                      plotALayers, plotBLayers, 
                                                                      plotCLayers) {

                      ## Split the dataset by sample to assign color names
                      df <- split(dataset, f=dataset$sample)
                      dataset <- data.table::rbindlist(lapply(df, function(x, sampleColor){
                          if (nrow(x) > 0) {
                              sampleName <- as.character(x$sample[1])
                              if (!is.null(sampleColor)) {
                                  x$sampleColor <- sampleColor[which(names(sampleColor) == sampleName)]
                              }
                              if (is.null(sampleColor)) {
                                  x$sampleColor <- sampleName
                              }
                          }
                          return(x)
                      }, sampleColor=sampleColor))
                      
                      colnames(dataset) <- c("Chromosome", "Position", "Chromosome2", "Position2", "Direction",
                                             "SV_Type", "Total_Read_Support", "Sample", "Genes", 
                                             "Midpoint", "chr_combo", "gene", "sampleColor")
                      ## Sort dataset
                      dataset$Sample <- factor(dataset$Sample, levels=gtools::mixedsort(unique(dataset$Sample)))
                      
                      ## Create bins for the chr positions (remove position transformation)
                      chrOrder <- gtools::mixedsort(unique(c(as.character(dataset$Chromosome), as.character(dataset$Chromosome2))))
                      chr1Length <- max(dataset[Direction=="cytoband" & Chromosome == chrOrder[1], Position2])
                      ## Get chr1 data
                      chr1OldBreaks <- round(seq(0, chr1Length, by=chr1Length/5), digits=0) 
                      chr1NewBreaks <- round(seq(0, chr1Length, by=chr1Length/5), digits=0) 
                      temp <- data.table(chr=chrOrder[1], newBreaks=chr1NewBreaks, oldBreaks=chr1OldBreaks)
                      ## Get chr2 data
                      if (length(chrOrder) > 1) {
                          chr2Start <- min(dataset[Direction=="cytoband" & Chromosome == chrOrder[2], Position])
                          chr2Length <- max(dataset[Direction=="cytoband" & Chromosome == chrOrder[2], Position2])
                          chr2OldBreaks <- round(seq(chr2Start, chr2Length, by=(chr2Length-chr2Start)/5), digits=0)
                          chr2NewBreaks <- round(seq(0, chr2Length-chr1Length, by=(chr2Length-chr1Length)/5), digits=0)
                          chr2 <- data.table(chr=chrOrder[2], newBreaks=chr2NewBreaks, oldBreaks=chr2OldBreaks)
                          temp <- rbind(temp, chr2)
                      }
                      
                      ## Get the start and stop for each chromosome
                      chr1End <- chr1Length
                      chr2End <- chr2Length
                      boundaries <- data.table(start=c(0, chr2Start), end=c(chr1End, chr2End))
                      
                      ##############################################################
                      ##### Plot the chromosome plot ###############################
                      ##############################################################
                      ## Get the cytoband data
                      coi <- dataset[Direction=="cytoband" & SV_Type != "centromere"]
                      coi$type <- "Chromosome"
                      suppressWarnings(coi$Height <- c(0.4, 0.6))
                      suppressWarnings(coi$color <- cytobandColor)
                      coi <- coi[!duplicated(coi$SV_Type)]
                      chrPlot <- ggplot() + geom_rect(data=coi, mapping=aes_string(xmin='Position',
                                                                                   xmax='Position2',
                                                                                   ymin=0,
                                                                                   ymax=1)) +
                          facet_grid(type ~ ., scales="fixed", space="fixed") +
                          scale_x_continuous(expand=c(0,0), limits=c(-5000000, max(coi$Position2) + 5000000)) + 
                          scale_y_continuous(expand=c(0,0)) + 
                          theme_bw() + 
                          geom_rect(data=coi, aes(xmin=Position, xmax=Position2, ymin=0, ymax=1, fill=Chromosome), fill=coi$color) +
                          geom_text(data=coi, aes(x=Midpoint, y=Height, label=SV_Type), angle=90, size=3) +
                          geom_vline(data=boundaries, aes(xintercept=start), linetype=2, color="Grey") + 
                          geom_vline(data=boundaries, aes(xintercept=end), linetype=2, color="Grey") + 
                          plotBLayers
                      
                      ## Get the centromeres
                      centromeres <- dataset[SV_Type=="centromere", c("Chromosome", "Position", "Position2")]
                      positions <- data.table::rbindlist(apply(centromeres, 1, function(x) {
                          midpoint <- round((as.numeric(as.character(x['Position'])) + 
                                                 as.numeric(as.character(x['Position2'])))/2, digits=0)
                          leftPositions <- data.table(x=c(x['Position'], midpoint, x['Position']),
                                                      y=c(0.10, 0.5, 0.90),
                                                      SV_Type="Chromosome",
                                                      id=paste(x['Chromosome'], "_left", sep=""))
                          rightPositions <- data.table(x=c(x['Position2'], midpoint, x['Position2']),
                                                       y=c(0.10, 0.5, 0.90),
                                                       SV_Type="Chromosome",
                                                       id=paste(x['Chromosome'], "_right", sep=""))
                          positions <- rbind(leftPositions, rightPositions)
                          positions$x <- as.numeric(positions$x)
                          return(positions)
                      }))
                      chrPlot <- chrPlot + geom_polygon(data=positions, mapping=aes(x=x, y=y, group=id), fill="red")
                      
                      ## Get the available SV events
                      availableSvTypes <- unique(dataset$SV_Type[-which(dataset$Direction=="cytoband")])
                      
                      ## Subset svWindow dataset to get DEL/DUP/INV/etc... and TRA/BND/etc...
                      sameChrSvWindow <- dataset[SV_Type=="DEL" | SV_Type=="DUP" | SV_Type =="INV" | SV_Type == "INS"]
                      sameChrSvWindow$SV_size <- sameChrSvWindow$Position2 - sameChrSvWindow$Position
                      diffChrSvWindow <- dataset[SV_Type=="BND" | SV_Type=="TRA"]
                      
                      ## Get the dataset for the gene text annotations
                      dataset <- dataset[Direction!="cytoband"]
                      gene_text <- dataset[,c("Midpoint", "Total_Read_Support", "gene", "SV_Type")]
                      ## If there is sv events for translocations, get the gene text
                      if (any(availableSvTypes %in% c("TRA", "BND"))) {
                          gene_text$Total_Read_Support[which(gene_text$SV_Type=="TRA")] <- 
                              as.numeric(as.character(gene_text$Total_Read_Support[which(gene_text$SV_Type=="TRA")])) + 
                              max(as.numeric(as.character(gene_text$Total_Read_Support[which(gene_text$SV_Type=="TRA")])))*.05   
                      }
                      ## If there is sv events for non-translocations, get the gene text
                      if (any(availableSvTypes %in% c("DEL", "DUP", "INV", "INS"))) {
                          gene_text$Total_Read_Support[which(gene_text$SV_Type!="TRA")] <- 
                              as.numeric(as.character(gene_text$Total_Read_Support[which(gene_text$SV_Type!="TRA")])) + 
                              max(as.numeric(as.character(gene_text$Total_Read_Support[which(gene_text$SV_Type!="TRA")])))*.05                          
                      }
                      gene_text$Total_Read_Support <- as.numeric(gene_text$Total_Read_Support)
                      gene_text <- gene_text[!duplicated(gene_text)]
                      if (!is.null(plotSpecificGene)) {
                          genes <- paste(plotSpecificGene, collapse="|")
                          gene_text <- gene_text[grep(genes, as.character(gene_text$gene))]
                          if (nrow(gene_text) == 0){
                              if (verbose) {
                                  message(paste0("The genes: ", plotSpecificGene, 
                                                 " could not be found. No genes will be shown on the plot."))
                              }
                              gene_text <- NULL
                              plotTraGenes <- FALSE
                              plotOtherGenes <- FALSE
                          }
                          if (!is.null(gene_text)) {
                              plotTraGenes <- any(c("TRA", "BND") %in% gene_text$SV_Type)
                              plotOtherGenes <- any(gene_text$SV_Type %in% c("DEL", "DUP", "INV", "INS"))
                          }
                      }
                      
                      ##########################################################
                      ##### Plot the translocation data ########################
                      ##########################################################
                      ## Get the start/end of chromosomes in the dataset if there is translocation data
                      ## TODO: Allow this to occur for intra-chromosomal translocations
                      if (any(availableSvTypes %in% c("TRA", "BND"))){
                          beziers <- data.frame(data.table::rbindlist(apply(diffChrSvWindow, 1, function(x) {
                              leftEnd <- data.table(position=as.numeric(x[2]), total_read_support=0, point="end", 
                                                    type="cubic", 
                                                    group=paste(as.character(x[2]), as.character(x[4]), as.character(x[5]),
                                                                 as.character(x[7]), as.character(x[8]), sep="_"), 
                                                    Sample=x[8], SV_Type=x[6],
                                                    Direction=x[5], sampleColor=x[13])
                              top <- data.table(position=as.numeric(x[10]), total_read_support=as.numeric(x[7])*2, point="control", 
                                                type="cubic", 
                                                group=paste(as.character(x[2]), as.character(x[4]), as.character(x[5]),
                                                             as.character(x[7]), as.character(x[8]), sep="_"), 
                                                Sample=x[8], SV_Type=x[6],
                                                Direction=x[5], sampleColor=x[13])
                              rightEnd <- data.table(position=as.numeric(x[4]), total_read_support=0, point="end", 
                                                     type="cubic", 
                                                     group=paste(as.character(x[2]), as.character(x[4]), as.character(x[5]),
                                                                  as.character(x[7]), as.character(x[8]), sep="_"),
                                                     Sample=x[8], SV_Type=x[6],
                                                     Direction=x[5], sampleColor=x[13])
                              final <- rbind(leftEnd, top, rightEnd)
                              return(final)
                          })))
                          beziers <- beziers[!duplicated(beziers),]
                          
                          beziers$Sample <- factor(beziers$Sample, levels=gtools::mixedsort(unique(beziers$Sample)))
                          
                          ## Plot the translocation data
                          traPlot <- ggplot() + geom_bezier(data=beziers, 
                                                            mapping=aes_string(x='position', y='total_read_support', group='group', 
                                                                               color='Sample', linetype='Direction')) +
                              facet_grid(SV_Type ~ ., scales="fixed", space="fixed") +
                              scale_x_continuous(expand=c(0,0), limits=c(-5000000, max(coi$Position2) + 5000000), 
                                                 breaks=temp$oldBreaks, labels=temp$newBreaks) + 
                              scale_y_continuous() +
                              geom_vline(data=boundaries, aes(xintercept=start), linetype=2, color="Grey") + 
                              geom_vline(data=boundaries, aes(xintercept=end), linetype=2, color="Grey") + 
                              theme_bw() + plotALayers + ylab("Total Read Support") +
                              geom_point(data=beziers[which(beziers$point=="control"),c("position","total_read_support", "Sample")], 
                                         aes(x=position, y=total_read_support/2, color=Sample)) 
                          if (plotTraGenes & !is.null(gene_text)) {
                              traPlot <- traPlot + geom_text(data=gene_text[SV_Type%in%c("TRA", "BND")], 
                                                             mapping=aes_string(x='Midpoint', y='Total_Read_Support', label='gene'))
                          }
                          ## Assign colors to sample
                          if (!is.null(sampleColor)) {
                              traPlot <- traPlot + scale_color_manual(name="Sample", values=sampleColor)
                          }
                      }
                      
                      ##############################################################
                      ##### Plot the non TRA sv events #############################
                      ##############################################################
                      if (any(availableSvTypes %in% c("DEL", "DUP", "INV", "INS"))) {
                          maxY <- max(as.numeric(as.character(sameChrSvWindow$Total_Read_Support))) + 30
                          sameChrSvWindow$Total_Read_Support <- as.numeric(sameChrSvWindow$Total_Read_Support)
                          nonTraPlot <- ggplot() + geom_point(data=sameChrSvWindow,
                                                              mapping=aes_string(x='Midpoint', y='Total_Read_Support', 
                                                                                 color="Sample"), size=2.5, alpha=0.75) + 
                              facet_grid(SV_Type ~ ., scales="fixed", space="fixed") +
                              scale_x_continuous(expand=c(0,0), limits=c(-5000000, max(coi$Position2) + 5000000), 
                                                 breaks=temp$oldBreaks, labels=temp$newBreaks) + 
                              scale_y_continuous(limits=c(0,maxY+maxY*0.05)) + 
                              geom_vline(data=boundaries, aes(xintercept=start), linetype=2, color="Grey") + 
                              geom_vline(data=boundaries, aes(xintercept=end), linetype=2, color="Grey") + 
                              theme_bw() + plotCLayers + ylab("Total Read Support") + xlab("Position")
                          if (plotOtherGenes & !is.null(gene_text)) {
                              nonTraPlot <- nonTraPlot + geom_text(data=gene_text[SV_Type%in%c("DEL", "DUP", "INV", "INS")], 
                                                                   mapping=aes_string(x='Midpoint', y='Total_Read_Support', label='gene')) 
                          }
                          ## Assign colors to sample
                          if (!is.null(sampleColor)) {
                              nonTraPlot <- nonTraPlot + scale_color_manual(name="Sample", values=sampleColor)
                          }
                      }
                      
                      ##############################################################
                      ##### Combine the 3 plots ####################################
                      ##############################################################
                      ## obtain the max width for relevant plots
                      if (nrow(diffChrSvWindow) > 0 & nrow(sameChrSvWindow) > 0) {
                          traPlot <- ggplotGrob(traPlot)
                          chrPlot <- ggplotGrob(chrPlot)
                          nonTraPlot <- ggplotGrob(nonTraPlot)
                          plotList <- list(traPlot, chrPlot, nonTraPlot)
                          sectionHeightsFinal <- sectionHeights
                      }
                      if (nrow(diffChrSvWindow) > 0 & nrow(sameChrSvWindow) == 0) {
                          traPlot <- ggplotGrob(traPlot)
                          chrPlot <- ggplotGrob(chrPlot)
                          plotList <- list(traPlot, chrPlot)
                          sectionHeightsFinal <- c(sectionHeights[1], sectionHeights[2])
                      }
                      if (nrow(diffChrSvWindow) == 0 & nrow(sameChrSvWindow) > 0) {
                          chrPlot <- ggplotGrob(chrPlot)
                          nonTraPlot <- ggplotGrob(nonTraPlot)
                          plotList <- list(chrPlot, nonTraPlot)
                          sectionHeightsFinal <- c(sectionHeights[2], sectionHeights[3])
                      }
                      plotList <- plotList[lapply(plotList, length) > 0]
                      plotWidths <- lapply(plotList, function(x) x$widths)
                      maxWidth <- do.call(grid::unit.pmax, plotWidths)
                      
                      ## Set the widths for all plots
                      for (i in 1:length(plotList)) {
                          plotList[[i]]$widths <- maxWidth
                      }
                      
                      ## Arrange the final plot
                      p1 <- do.call(gridExtra::arrangeGrob, 
                                    c(plotList, list(ncol=1, heights=sectionHeightsFinal)))
                      plot(p1)
                      
                      return(p1)
                  }, 
                  plotSpecificGene=plotSpecificGene, cytobandColor=cytobandColor, 
                  sectionHeights=sectionHeights, sampleColor=sampleColor, 
                  plotTraGenes=plotTraGenes, plotOtherGenes=plotOtherGenes, 
                  plotALayers=plotALayers, plotBLayers=plotBLayers, plotCLayers=plotCLayers))
                  
                  return(svPlots)
              }
              
              
          })