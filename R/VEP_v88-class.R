################################################################################
##################### Public/Private Class Definitions #########################

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Public Class !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

#' Class VEP_v88
#' 
#' An S4 class to represent data in variant effect predictor version 88 format,
#' inherits from the VEP_Virtual class, under development!!!
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
             expecPositionNames <- c("Location")
             expecMutationNames <- c("Allele", "Consequence")
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

#' Constructor for the VEP_v88 sub-class
#' 
#' @name VEP_v88
#' @rdname VEP_v88-class
#' @param vepData data.table object containing a VEP annotation file conforming
#' to the version 88 specifications.
#' @param vepHeader Object of class list containing character vectors for vep
#' header information.
VEP_v88 <- function(vepData, vepHeader){
    
    # set the columns descriptions for the object
    if(length(vepHeader) == 0){
        description <- data.table::data.table()
    } else {
        description <- parseDescription(vepHeader)
    }
    
    # set the header for the object
    if(length(vepHeader) == 0){
        header <- data.table::data.table()
    } else {
        header <- parseHeader(vepHeader)
    }
    
    # convert the "extra" field in vepData to separate columns
    vepData <- parseExtra(vepData)
    
    positionColNames <- c("Location")
    position <- vepData[,positionColNames, with=FALSE]
    
    mutationColNames <- c("Allele", "Consequence")
    mutation <- vepData[,mutationColNames, with=FALSE]
    
    sampleColNames <- c("sample")
    sample <- vepData[,sampleColNames, with=FALSE]
    
    metaColNames <- !colnames(vepData) %in% c(positionColNames, mutationColNames, sampleColNames)
    meta <- vepData[,metaColNames, with=FALSE]
    
    new("VEP_v88", header=header, description=description, position=position, mutation=mutation, sample=sample, meta=meta)
}
