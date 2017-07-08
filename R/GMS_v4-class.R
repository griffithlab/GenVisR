#' Class GMS_v4
#' 
#' An S4 class to represent data in gms annotation version 4,
#' inherits from the GMS_Virtual class.
#' @name GMS_v4-class
#' @rdname GMS_v4-class
#' @slot position data.table object containing column names "chromosome_name",
#' "start", "stop".
#' @slot mutation data.table object containing column names "reference",
#' "variant", "trv_type".
#' @slot sample data.table object containing columns names "sample".
#' @slot meta data.table object containing meta data.
#' @include GMS_Virtual-class.R
#' @import methods

setClass("GMS_v4",
         contains="GMS_Virtual",
         validity=function(object){
             cat("!!!!! GMS_v4~Inspector !!!!!\n")
             expecPositionNames <- c("chromosome_name", "start", "stop")
             expecMutationNames <- c("reference", "variant", "trv_type")
             expecSampleNames <- c("sample")
             expecMetaNames <- c("gene_name")
             
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
             if(!all(expecMetaNames %in% colnames(object@meta))){
                 memo <- paste("Missing the following required columns in slot meta:",
                               toString(expecMetaNames[!expecMetaNames %in% colnames(object@meta)]))
                 stop(memo)
             }
             return(TRUE)
         }
)

#' Initalizer method for the GMS_v4 sub-class
#' 
#' @name GMS_v4
#' @param .Object object of class GMS_v4
#' @noRd
#' @rdname GMS_v4-class
setMethod(
    f="initialize",
    signature="GMS_v4",
    definition=function(.Object, gmsData){
        
        cat("!!!!! GMS_v4~Initalizer !!!!!\n")
        positionColNames <- c("chromosome_name", "start", "stop")
        .Object@position <- gmsData[,positionColNames, with=FALSE]
        
        mutationColNames <- c("reference", "variant", "trv_type")
        .Object@mutation <- gmsData[,mutationColNames, with=FALSE]
        
        sampleColNames <- c("sample")
        .Object@sample <- gmsData[,sampleColNames, with=FALSE]
        
        metaColNames <- !colnames(gmsData) %in% c(positionColNames, mutationColNames, sampleColNames)
        .Object@meta <- gmsData[,metaColNames, with=FALSE]
        
        validObject(.Object)
        return(.Object)
    }
)

#' Constructor for the GMS_v4 sub-class
#' 
#' @name GMS_v4
#' @rdname GMS_v4-class
#' @param gmsData data.table object containing a gms annotation file conforming
#' to the version 4 specifications.
GMS_v4 <- function(gmsData){
    cat("!!!!! GMS_v4~Constructor !!!!!\n")
    new("GMS_v4", gmsData=gmsData)
}