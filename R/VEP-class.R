#' Class VEP
#' 
#' An S4 class for Variant Effect Predictor input.
#' @name VEP-class
#' @rdname VEP-class
#' @slot path Character string specifying the paths of the VEP files read in.
#' @slot version Numeric value specifying the version of VEP used.
#' @slot vepObject vep object which inherits from VEP_Virtual class.
#' @exportClass VEP
#' @include VEP_Virtual-class.R
#' @import methods
setClass("VEP",
         representation=representation(path="character",
                                       version="numeric",
                                       vepObject="VEP_Virtual"))

#' Initalizer method for the VEP container class
#' 
#' @name VEP
#' @rdname VEP-class
#' @importFrom data.table fread
#' @importFrom data.table rbindlist
setMethod(
    f="initialize",
    signature="VEP",
    definition=function(.Object, path, version, verbose){
        cat("!!!!! VEP~Initalizer !!!!!\n")
        
        # Grab all files and assign to slot
        vepFiles <- Sys.glob(path)
        .Object@path <- vepFiles
        
        # anonymous function to read in files
        a <- function(x, verbose){
            # detect OS and remove slashes and extension
            if(.Platform$OS.type == "windows"){
                sampleName <- gsub("(.*/)||(.*\\\\)", "", x)
                sampleName <- gsub("\\.[^.]+$", "", x)
            } else {
                sampleName <- gsub("(.*/)", "", x)
                sampleName <- gsub("\\.[^.]+$", "", sampleName)
            }
            # read the header
            header <- readLines(con=x, n=400)
            header <- header[grepl("^##", header)]
            # find where headers stop and read the data
            skip <- length(header)
            vepData <- suppressWarnings(data.table::fread(input=x,
                                                          stringsAsFactors=TRUE,
                                                          verbose=verbose,
                                                          skip=skip))
            # set sample if it's not already in the data table
            if(any(colnames(vepData) %in% "sample")){
                return(vepData)
            } else {
                vepData$sample <- sampleName
                return(list("data"=vepData, "header"=header))
            }
        }
       
        # aggregate data into a single data table if necessary
        if(length(vepFiles) == 0){
            memo <- paste("No files found using:", path)
            stop(memo)
        } else {
            # Read in the information
            vepInfo <- lapply(vepFiles, a, verbose)
            
            # extract header and data information
            vepHeader <- lapply(vepInfo, function(x) x[["header"]])
            vepData <- lapply(vepInfo, function(x) x[["data"]])
            
            # aggregate the data
            vepData <- data.table::rbindlist(vepData, fill=TRUE)
        } 
        
        # grab the version and assign it
        a <- function(x){
            # find the element which defines the VEP version
            x <- x[grepl("VARIANT EFFECT PREDICTOR", x)]
            
            # extract the version
            x <- regmatches(x,regexpr("[0-9]+\\.*[0-9]*",x))
            
            if(length(x) != 1){
                memo <- paste("Expected 1 entry for VEP version, found:",
                              length(x), "using", as.numeric(x[1]))
                warning(memo)
            }
            return(as.numeric(x[1]))
        }
        if(version == "auto"){
            version <- lapply(vepHeader, a)
            version <- unique(unlist(version))
            if(length(version) > 1){
                version <- version[1]
                memo <- paste("Expect 1 version, the following versions were",
                              "found:", toString(version), "Using version",
                              version, "for parsing!")
                warning(memo)
            } else if(length(version) == 0){
                memo <- paste("Cannot infer version from vep headers",
                              "no versions found!")
                stop(memo)
            }
            .Object@version <- version
        } else {
            .Object@version <- version
        }
        
        # assign the vepData to it's slot
        if(version >= 88 || version < 89){
            .Object@vepObject <- VEP_v88(vepData=vepData, vepHeader=vepHeader)
        } else {
            memo <- paste("Currently only VEP version 88 is supported, make a",
                          "feature request on",
                          "https://github.com/griffithlab/GenVisR!")
            stop(memo)
        }
        
        return(.Object)
    })

#' Constructor for the VEP container class.
#' 
#' @name VEP
#' @rdname VEP-class
#' @param path String specifying the path to a VEP annotation file. Can accept
#' wildcards if multiple VEP annotation files exist (see details).
#' @param version String specifying the version of the VEP files, Defaults to
#' auto which will look for the version in the header.
#' @param verbose Boolean specifying if progress should be reported while
#' reading in the VEP files.
#' @details When specifying a path to a VEP annotation file the option exist to
#' either specify the full path to an annotation file or to use wildcards to
#' specify multiple files. When specifying a full path the initalizer will check
#' if a column named "sample" containg the relevant sample for each row exists.
#' If such a column is not found the initalizer will assume this file
#' corresponds to only one sample and populate a sample column accordingly.
#' Alternatively if multiple files are specified at once using a wildcard, the
#' initalizer will aggregate all the files and use the file names minus any
#' extension to populate sample names.
#' @export
VEP <- function(path, version="auto", verbose=FALSE){
    cat("!!!!! VEP~Constructor !!!!!\n")
    new("VEP", path=path, version=version, verbose=verbose)}

#' @rdname writeData-methods
#' @aliases writeData,VEP
#' @param object Object of class VEP.
#' @param file Character string specifying a file to send output to.
#' @param sep Delimiter used when writing output, defaults to tab.
setMethod(f="writeData",
          signature="VEP",
          definition=function(object, file, ...){
              writeData(object@vepObject, file, sep="\t")
          })