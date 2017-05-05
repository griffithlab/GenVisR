#' Class VEP_v88
#' 
#' An S4 class to represent data in variant effect predictor version 88 format,
#' inherits from the VEP_Virtual class.
#' @name VEP_v88-class
#' @rdname VEP_v88-class
#' @slot header data.table object containing header information
#' @slot description data.table object containing column descriptions
#' @slot position data.table object containing column names "chromosome_name",
#' "start", "stop".
#' @slot mutation data.table object containing column names "reference",
#' "variant", "trv_type".
#' @slot sample data.table object containing columns names "sample".
#' @slot meta data.table object containing meta data.
#' @include GMS_Virtual-class.R
#' @import methods

setClass("VEP_v88",
         contains="VEP_Virtual",
         validity=function(object){
             cat("!!!!! VEP_v4~Inspector !!!!!\n")
             expecPositionNames <- c("chromosome_name", "start", "stop")
             expecMutationNames <- c("reference", "variant", "trv_type")
             expecSampleNames <- c("sample")
             
             if(!all(expecPositionNames %in% colnames(object@position))){
                 memo <- paste("Missing the following required columns in slot position:",
                               toString(expecPositionNames[!expecPositionNames %in% colnames(object@position)]))
                 stop(memo)
             }
             if(!all(expecMutationNames %in% colnames(object@mutation))){
                 memo <- paste("Missing the following required columns in slot mutation:",
                               toString(expecMutationNames[!expecMutationNames %in% colnames(object@mutation)]))
                 stop(memo)
             }
             if(!all(expecSampleNames %in% colnames(object@sample))){
                 memo <- paste("Missing the following required columns in slot sample:",
                               toString(expecSampleNames[!expecSampleNames %in% colnames(object@sample)]))
                 stop(memo)
             }
             return(TRUE)
         }
)

#' Initalizer method for the VEP_v88 sub-class
#' 
#' @name VEP_v88
#' @rdname VEP_v88-class
setMethod(
    f="initialize",
    signature="VEP_v88",
    definition=function(.Object, vepData, vepHeader){
        
        cat("!!!!! VEP_v88~Initalizer !!!!!\n")
        
        # set the columns descriptions for the object
        .Object@description <- parseDescription(.Object, vepHeader)
        
        # set the header for the object
        .Object@header <- parseHeader(.Object, vepHeader)
        
        # convert the "extra" field in vepData to separate columns
        vepData <- parseExtra(.Object, vepData)
        
        positionColNames <- c("chromosome_name", "start", "stop")
        .Object@position <- vepData[,positionColNames, with=FALSE]
        
        mutationColNames <- c("reference", "variant", "trv_type")
        .Object@mutation <- vepData[,mutationColNames, with=FALSE]
        
        sampleColNames <- c("sample")
        .Object@sample <- vepData[,sampleColNames, with=FALSE]
        
        metaColNames <- !colnames(vepData) %in% c(positionColNames, mutationColNames, sampleColNames)
        .Object@meta <- vepData[,metaColNames, with=FALSE]
        
        validObject(.Object)
        return(.Object)
    }
)

#' Constructor for the VEP_v4 sub-class
#' 
#' @name VEP_v88
#' @rdname VEP_v88-class
#' @param gmsData data.table object containing a VEP annotation file conforming
#' to the version 88 specifications.
#' @param vepHeader Object of class list containing character vectors for vep
#' header information.
VEP_v88 <- function(vepData, vepHeader){
    cat("!!!!! VEP_v88~Constructor !!!!!\n")
    new("VEP_v88", vepData=vepData, vepHeader=vepHeader)
}
