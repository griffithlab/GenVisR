#' Class MutationAnnotationFormat
#' 
#' An S4 class acting as a container for MutationAnnotationFormat version sub-classes.
#' @name MutationAnnotationFormat-class
#' @rdname MutationAnnotationFormat-class
#' @slot path Character string specifying the path of the MAF file read in.
#' @slot version Numeric value specifying the version of the MAF file.
#' @slot mafObject MutationAnnotationFormat object which inherits from
#' MutationAnnotationFormat_Virtual class.
#' @exportClass MutationAnnotationFormat
#' @include MutationAnnotationFormat_Virtual-class.R
#' @import methods

setClass("MutationAnnotationFormat",
         representation=representation(path="character",
                                       version="numeric",
                                       mafObject="MutationAnnotationFormat_Virtual")
)

#' Initalizer method for the MutationAnnotationFormat container class
#' 
#' @name MutationAnnotationFormat
#' @rdname MutationAnnotationFormat-class
#' @importFrom data.table fread

setMethod(
    f="initialize",
    signature="MutationAnnotationFormat",
    definition=function(.Object, path, version, verbose){
        
        cat("!!!!! MutationAnnotationFormat~Initalizer !!!!!\n")

        mafData <- suppressWarnings(data.table::fread(input=path,
                                                      stringsAsFactors=TRUE,
                                                      verbose=verbose))
        # grab the maf version
        if(toupper(version) == "AUTO"){
            # read the version
            mafVersion <- readLines(con=path, n=50)
            mafVersion <- mafVersion[which(grepl("^#version", mafVersion))]
            mafVersion <- as.numeric(as.character(gsub("#version\\W", "", mafVersion)))
            if(length(mafVersion) == 0) stop("Unable to infer the maf Version from header, please specify")
        } else {
            mafVersion <- version
        }

        ##### Obtain the appropriate MutationAnnotationFormat sub-class ######
        if(mafVersion == 1.0){
            mafObject <- MutationAnnotationFormat_v1.0(mafData=mafData)
        }else if(mafVersion == 2.0){
            mafObject <- MutationAnnotationFormat_v2.0(mafData=mafData)
        }else if(mafVersion == 2.1){
            mafObject <- MutationAnnotationFormat_v2.1(mafData=mafData)
        }else if(mafVersion == 2.3){
            mafObject <- MutationAnnotationFormat_v2.3(mafData=mafData)
        }else if(mafVersion == 2.4){
            mafObject <- MutationAnnotationFormat_v2.4(mafData=mafData)
        }else{
            memo <- paste("The maf version:", toString(mafVersion),
                          "is currently unsupported. Make a feature request on",
                          "https://github/griffithlab/GenVisR!")
            stop(memo)
        }
        browser()
        .Object@path <- path
        .Object@version <- mafVersion
        .Object@mafObject <- mafObject
        
        return(.Object)
    }
)

#' Constructor for the MutationAnnotationFormat container class.
#' 
#' @name MutationAnnotationFormat
#' @rdname MutationAnnotationFormat-class
#' @param path String specifying the path to a MAF file.
#' @param version String specifying the version of the MAF file, if set to auto
#' the version will be obtained from the header in the MAF file.
#' @param verbose Boolean specifying if progress should be reported while reading
#' in the MAF file.
#' @export

MutationAnnotationFormat <- function(path, version="auto", verbose=FALSE){
    cat("!!!!! MutationAnnotationFormat~Constructor !!!!!\n")
    new("MutationAnnotationFormat", path=path, version=version, verbose=verbose)
}
