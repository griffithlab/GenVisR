#' Class MutationAnnotationFormat_v2.1
#' 
#' An S4 class to represent data in mutation annotation format version 2.1,
#' inherits from the MutationAnnotationFormat_Virtual class.
#' @name MutationAnnotationFormat_v2.1-class
#' @rdname MutationAnnotationFormat_v2.1-class
#' @slot position data.table object containing column names "Chromosome",
#' "Start_Position", "End_Position, "Strand".
#' @slot mutation data.table object containing column names "Variant_Classification",
#' "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2".
#' @slot sample data.table object containing columns names "Tumor_Sample_Barcode".
#' @slot meta data.table object containing meta data.
#' @include MutationAnnotationFormat_Virtual-class.R
#' @import methods

setClass("MutationAnnotationFormat_v2.1",
         contains="MutationAnnotationFormat_Virtual",
         validity=function(object){
             cat("!!!!! MutationAnnotationFormat_v2.1~Inspector !!!!!\n")
             expecPositionNames <- c("Chromosome", "Start_Position", "End_Position", "Strand")
             expecMutationNames <- c("Variant_Classification", "Variant_Type", "Reference_Allele",
                                     "Tumor_Seq_Allele1", "Tumor_Seq_Allele2")
             expecSampleNames <- c("Tumor_Sample_Barcode")
             
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

#' Initalizer method for the MutationAnnotationFormat_v2.1 sub-class
#' 
#' @name MutationAnnotationFormat_v2.1
#' @rdname MutationAnnotationFormat_v2.1-class
#' @noRd
#' @param .Object object of class MutationAnnotationFormat_v2.1

setMethod(
    f="initialize",
    signature="MutationAnnotationFormat_v2.1",
    definition=function(.Object, mafData){
        
        cat("!!!!! MutationAnnotationVersion_v2.1~Initalizer !!!!!\n")
        positionColNames <- c("Chromosome", "Start_Position", "End_Position", "Strand")
        .Object@position <- mafData[,positionColNames, with=FALSE]
        
        mutationColNames <- c("Variant_Classification", "Variant_Type", "Reference_Allele",
                              "Tumor_Seq_Allele1", "Tumor_Seq_Allele2")
        .Object@mutation <- mafData[,mutationColNames, with=FALSE]
        
        sampleColNames <- c("Tumor_Sample_Barcode")
        .Object@sample <- mafData[,sampleColNames, with=FALSE]
        
        metaColNames <- !colnames(mafData) %in% c(positionColNames, mutationColNames, sampleColNames)
        .Object@meta <- mafData[,metaColNames, with=FALSE]
        
        validObject(.Object)
        return(.Object)
    }
)

#' Constructor for the MutationAnnotationFormat_v2.1 sub-class
#' 
#' @name MutationAnnotationFormat_v2.1
#' @rdname MutationAnnotationFormat_v2.1-class
#' @param mafData data.table object containing a maf file conforming to the
#' version 2.1 specification.

MutationAnnotationFormat_v2.1 <- function(mafData){
    cat("!!!!! MutationAnnotationFormat_v2.1~Constructor !!!!!\n")
    new("MutationAnnotationFormat_v2.1", mafData=mafData)
}

