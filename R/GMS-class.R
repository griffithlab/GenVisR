#' Class GMS
#' 
#' An S4 class for Genome Modeling System annotation files.
#' @name GMS-class
#' @rdname GMS-class
#' @slot path Character string specifying the paths of the GMS files read in.
#' @slot version Numeric value specifying the version of the GMS annotation files.
#' @slot gmsObject gms object which inherits from gms_Virtual class.
#' @exportClass GMS
#' @include GMS_Virtual-class.R
#' @import methods
setClass("GMS",
         representation=representation(path="character",
                                       version="numeric",
                                       gmsObject="GMS_Virtual")
)

#' Initalizer method for the GMS container class
#' 
#' @name GMS
#' @rdname GMS-class
#' @importFrom data.table fread
#' @importFrom data.table rbindlist
setMethod(
    f="initialize",
    signature="GMS",
    definition=function(.Object, path, version, verbose){
        cat("!!!!! GMS~Initalizer !!!!!\n")
        
        # Grab all files
        gmsFiles <- Sys.glob(path)
        
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
            # read data
            gmsData <- suppressWarnings(data.table::fread(input=x,
                                                          stringsAsFactors=TRUE,
                                                          verbose=verbose))
            # set sample if it's not already in the data table
            if(any(colnames(gmsData) %in% "sample")){
                return(gmsData)
            } else {
                gmsData$sample <- sampleName
                return(gmsData)
            }
        }
        
        # aggregate data into a single data table if necessary
        if(length(gmsFiles) == 0){
            memo <- paste("No files found using:", path)
            stop(memo)
        } else {
            gmsData <- lapply(gmsFiles, a, verbose)
            gmsData <- data.table::rbindlist(gmsData)
        } 
        browser()
         
    }
)

#' Constructor for the GMS container class.
#' 
#' @name GMS
#' @rdname GMS-class
#' @param path String specifying the path to a GMS annotation file. Can accept
#' wildcards if multiple GMS annotation files exist (see details).
#' @param version String specifying the version of the GMS files, Defaults to
#' version 4.
#' @param verbose Boolean specifying if progress should be reported while
#' reading in the GMS files.
#' @details When specifying a path to a GMS annotation file the option exist to
#' either specify the full path to an annotation file or to use wildcards to
#' specify multiple files. When specifying a full path the initalizer will check
#' if a column named "sample" containg the relevant sample for each row exists.
#' If such a column is not found the initalizer will assume this file
#' corresponds to only one sample and populate a sample column accordingly.
#' Alternatively if multiple files are specified at once using a wildcard, the
#' initalizer will aggregate all the files and use the file names minus any
#' extension top populate sample names.
#' The version defaults to 4 which is the default value of the GMS annotator.
#' This value will need to be changed only if files were created using a
#' different GMS annotator version.
#' @export
GMS <- function(path, version=4, verbose=FALSE){
    cat("!!!!! GMS~Constructor !!!!!\n")
    new("GMS", path=path, version=version, verbose=verbose)
}