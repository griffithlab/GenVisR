################################################################################
##################### Public/Private Class Definitions #########################

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Public Class !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

#' Class VariantCallFormat
#' 
#' @name VariantCallFormat-class
#' @rdname VariantCallFormat-class
#' @slot path Character string specifying the path of the VCF file read in
#' @slot version Character string specifiying the version of the vcf file
#' @slot vcfObject vcf object which inherits from VCF_Virtual class
#' @exportClass VariantCallFormat
#' @include VCF_Virtual-class.R
#' @import methods
setClass("VariantCallFormat", 
         representation=representation(path="character",
                                       version="character",
                                       vcfObject="VCF_Virtual"),
         validity=function(object) {
             
         })

#' Constructor for the VCF container class
#' @name VariantCallFormat
#' @rdname VariantCallFormat-class
#' @param path String specifying the path to a VCF file. Can accept wildcards 
#' if multiple VCF files exist (see details).
#' @param data data.table object storing a VCF file. Overrides "path" if 
#' specified
#' @param version String specifying the version of the VCF file, if set to auto
#' the version will be obtained from the header in the VCF file
#' @param svCaller String specifying the structural variant caller used
#' @param paired Boolean specifiying if the svCaller was run with the paired option 
#' (i.e. tumor-normal)
#' @param tumorColumn Integer specifying the column number with the tumor read support information. 
#' Only used when paired=TRUE. 
#' @param verbose Bolean specifying if progress should be reported while reading
#' in the VCF file
#' @details When specifying a path to a VCF file, the option exists to either
#' specify the full path to a vcf file or to us wildcards to specify multiple 
#' files. When specifying a full path, the initializer will check if a column
#' named "sample" containing the relevant sample for each row exists. If such a 
#' column is not found, the initializer will assume this file correspnds to 
#' only one sample and populate a sample column accordingly. Alternatively, if 
#' multiple files are specified at once using a wildcard, the initializer will 
#' aggregate all the files and use the filenames minus any extension to 
#' populate the "sample" column.
#' @importFrom data.table fread
#' @importFrom data.table rbindlist
#' @importFrom data.table data.table
#' @export
VariantCallFormat <- function(path=NULL, data=NULL, version="auto", svCaller=NULL, paired=paired, 
                                     tumorColumn=tumorColumn, verbose=FALSE) {
    
    ## Check if both path and data are both null
    if (is.null(path) & is.null(data)) {
        memo <- paste("Both the path and data variable cannot be null.")
        stop(memo)
    }
    
    ## Get the data if the dataset is not provided
    if (is.null(data)) {
        ## Add wildcard to the path if it not present
        if (length(path) > 1 & strsplit(path, split="/")[[1]][length(strsplit(path, split="/")[[1]])] != "*") {
            memo <- paste("No wildcard found in the designated path. Please add wildcard to the path. For example:",  
                          "~/Desktop/StructuralVariants/* to use files in ~/Desktop/StructuralVariants/ directory.")
            stop(memo)
        }
        
        ## Grab all the files and assign to slot
        vcfFiles <- Sys.glob(path)
        
        ## Anonymous function to read in files
        a1 <- function(x, verbose) {
            ## Detect OS and remove slashes and extension
            if (.Platform$OS.type == "windows") {
                sampleName <- gsub("(.*/)||(.*\\\\)", "", x)
                sampleName <- gsub("\\.[^.]+$", "", x)
            } else {
                sampleName <- gsub("(.*/)", "", x)
                sampleName <- gsub("\\.[^.]+$", "", sampleName)
            }
            
            ## Read the header
            header <- readLines(con=x, n=400)[grep("^##", readLines(con=x, n=400))]
            
            ## Find where the headers stops and read the data
            skip <- length(header)
            vcfData <- suppressWarnings(data.table::fread(input=x,
                                                     stringsAsFactors=TRUE,
                                                     verbose=verbose, 
                                                     skip=skip))
            
            ## Set sample if it is not already in the data table
            if(any(colnames(vcfData) %in% "sample")){
                return(vcfData)
            } else {
                vcfData$sample <- sampleName
                return(list("data"=vcfData, "header"=header))
            }
        }
        
        ## Aggregate data into a single data table if necessary
        if(length(vcfFiles)==0) {
            memo <- paste("No files found using:", path)
            stop(memo)
        } else {
            ## Read in the information
            vcfInfo <- lapply(vcfFiles, a1, verbose)
            
            ## Extract header and data information
            vcfHeader <- lapply(vcfInfo, function(x) x[["header"]])
            vcfData <- data.table::rbindlist(lapply(vcfInfo, function(x) x[["data"]]))
        }
    } 
    
    ## If a dataset is provided: 
    if (!is.null(data)) {
        path <- "none"
        
        ## Check to see if the VCF version is present
        if (version == "auto") {
            memo <- paste("If a dataset is loaded into the function, the version of ", 
                          "the VCF that came from the sv caller must be specified as a CHARACTER VECTOR. ",
                          "The current VCF versions that are supported are: 4.1 and 4.2", sep="")
            stop(memo)
        }
        
        ## Convert dataset to data.table class
        if (is.data.table(data)) {
            vcfHeader <- data.table::data.table()
            vcfData <- data 
        }
        if (!is.data.table(data)) {
            memo <- paste("data is not of class data.table,",
                          "attempting to coerce")
            warning(memo)
            vcfHeader <- data.table::data.table()
            vcfData <- data.table::data.table(data)
        }
        
        ## Check to see if it has the necessary columns
        cnames <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL",
                    "FILTER", "INFO", "FORMAT", "sample")
        if (any(all(cnames %in% colnames(vcfData))==FALSE)){
            colnames(vcfData) %in% cnames
            memo <- paste("The columns in the input dataset do not match the required columns:",
                          paste(cnames, collapse=", "), ". This does not include the", 
                          "columns with sample read support information.")
            stop(memo)
        }
    }
    
    ## Check to see if there are any columns with read support information
    cnames <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL",
                "FILTER", "INFO", "FORMAT", "sample")
    readSupportColumns <- colnames(vcfData)[which(!colnames(vcfData) %in% cnames)]
    if (!(length(readSupportColumns) == 1 | length(readSupportColumns)==2)) {
        memo <- paste("There are", length(readSupportColumns), "columns that are not one of the",
                      "required columns. There must be 1 column with sample read support",
                      "information for unpaired data and 2 columns with sample read support",
                      "information for paired data.")
        stop(memo)
    }
    
    
    ## Grab the version and assign it
    a2 <- function(x) {
        ## Find the element which defines the VEP version
        x <- x[grepl("fileformat=", x)]
        
        ## Extract the version
        x <- regmatches(x, regexpr("[0-9]+\\.*[0-9]*",x))
        
        if (length(x) != 1) {
            memo <- paste("Expected 1 entry for VCF version found:",
                          length(x), "using", as.numeric(x[1]))
            warning(memo)
        }
        return(as.numeric(x[1]))
    }
    if (version == "auto") {
        version <- unique(unlist(lapply(vcfHeader, a2)))
        if (length(version) > 1) {
            version <- as.character(version[1])
            memo <- paste("Expected one version, the following versions were",
                          "found:", toString(version), "Using version",
                          version, "for parsing!")
            warning(memo)
        } else if (length(version) == 0) {
            memo <- paste("Cannot infer version from vcf headers",
                          "no versions found!")
            stop(memo)
        }
    }
    
    ## Perform quality check for the svCaller
    if (!svCaller %in% c("Manta")) {
        memo <- paste0("The specified svCaller: ", svCaller, " is not supported. ",
                       "Only the following callers are support: Manta. ",
                       "Make sure the svCaller is one of the ",
                       "supported callers listed with proper capitalization and spelling.")
        stop(memo)
    }
    
    ## Assign the vcfData to its slot
    if (version == "4.1" & svCaller == "Manta") {
        vcfObject <- VCF_Manta_v4.1(vcfData=vcfData, vcfHeader=vcfHeader, paired=paired, tumorColumn=tumorColumn)
    } else if (version == "4.2" & svcaller =="Manta") {
        vcfObject <- VCF_Manta_v4.2(vcfData=vcfData, vcfHeader=vcfHeader, paired=paired, tumorColumn=tumorColumn)
    } else {
        memo <- paste("Currently only VCF versions 4.1 and 4.2 for Manta are supported,",
                      "make a feature request on", 
                      "https://github.com/griffithlab/GenVisR!")
        stop(memo)
    }
    
    ## Initialize the object
    new("VariantCallFormat", path=path, vcfObject=vcfObject, version=as.character(version))
}

################################################################################
####################### Method function definitions ############################

#' @rdname getVcfData-methods
#' @name getVcfData
#' @aliases getVcfData
#' @noRd
#' @importFrom data.table data.table
setMethod(f="getVcfData",
          signature="VariantCallFormat",
          definition=function(object, filter, maxSvSize, svType, 
                              verbose, ...) {
              
              ## Print status message
              if (verbose) {
                  memo <- paste0("converting ", class(object), " to expected ",
                                 "StructuralVariant format")
                  message(memo)
              }
              
              object <- object@vcfObject@vcfData
              
              ## Filter out sv calls that are not "PASS"
              if (filter == TRUE) {
                  object <- object[FILTER=="PASS"]
              }
              
              ## Remove large SV
              if (is.null(maxSvSize) == FALSE) {
                  ## Get the difference in positions
                  temp <- suppressWarnings(data.table::rbindlist(apply(object, 1, function(x, maxSvSize){
                      if (x["svtype"] == "BND" | x["svtype"] == "TRA"){
                          x$diff <- maxSvSize - 1
                      } else {
                          x$diff <- as.numeric(x["position2"]) - as.numeric(x["position"])
                      }
                      return(x)
                  }, maxSvSize=maxSvSize)))
                  
                  ## Perform the subset
                  object <- temp[diff < maxSvSize, c(1:15)]
              }
              
              ## Remove sv types that are not necessary
              available_svTypes <- unlist(as.vector(object$svtype))
              if (length(svType) > 0) {
                  ## Check to see if the SV type is in the data.table
                  ## Perform the subset if svtype is available
                  if (all(svType %in% available_svTypes)) {
                      object <- object[svtype %in% svType]
                  }
                  if (!all(svType %in% available_svTypes)) {
                      memo <- paste0("Desired svtype is not found. Make sure ",
                                     "the specified svType is one of: ", 
                                     paste(available_svTypes, collapse=", "))
                      stop(memo)
                  }
              }
              
              ## Stop the
              return(object)
            })
