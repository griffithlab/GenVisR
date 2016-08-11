#' Class MutationAnnotationFormat_v2.4
#' 
#' An S4 class to represent data in mutation annotation format version 2.4
#' @name MutationAnnotationFormat_v2.4-class
#' @rdname MutationAnnotationFormat_v2.4-class
#' @slot position test
setClass("MutationAnnotationFormat_v2.4",
         contains="MutationAnnotationFormat_Virtual",
         validity=function(object){
             cat("!!!!! MutationAnnotationFormat_v2.4~Inspector !!!!!\n")
             expecPositionNames <- c("Chromosome", "Start_Position", "End_Position", "Strand")
             expecMutationNames <- c("Variant_Classification", "Variant_Type", "Reference_Allele",
                                     "Tumor_Seq_Allele1", "Tumor_Seq_Allele2")
             expecSampleNames <- c("Tumor_Sample_Barcode")
             
             if(!all(expecPositionNames %in% colnames(object@position))){
                 memo <- paste("Missing the following required columns in slot position:",
                               toString(expected_col[!expecPositionNames %in% colnames(object@position)]))
                 stop(memo)
             }
             if(!all(expecMutationNames %in% colnames(object@mutation))){
                 memo <- paste("Missing the following required columns in slot mutation:",
                               toString(expected_col[!expecMutationNames %in% colnames(object@mutation)]))
                 stop(memo)
             }
             if(!all(expecSampleNames %in% colnames(object@sample))){
                 memo <- paste("Missing the following required columns in slot sample:",
                               toString(expected_col[!expecSampleNames %in% colnames(object@sample)]))
                 stop(memo)
             }
             return(TRUE)
         }
)

#' Initalizer method of MutationAnnotationFormat_v2.4 class
#' 
#' @name MutationAnnotationFormat_v2.4
#' @rdname MutationAnnotationFormat_v2.4-class
setMethod(
    f="initialize",
    signature="MutationAnnotationFormat_v2.4",
    definition=function(.Object, mafPath, mafVersion, mafData){

        cat("!!!!! MutationAnnotationVersion_v2.4~Initalizer !!!!!\n")
        positionColNames <- c("Chromosome", "Start_Position", "End_Position", "Strand")
        .Object@position <- mafData[,positionColNames, with=FALSE]
        
        mutationColNames <- c("Variant_Classification", "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2")
        .Object@mutation <- mafData[,mutationColNames, with=FALSE]
        
        sampleColNames <- c("Tumor_Sample_Barcode")
        .Object@sample <- mafData[,sampleColNames, with=FALSE]
        
        metaColNames <- !colnames(mafData) %in% c(positionColNames, mutationColNames, sampleColNames)
        .Object@meta <- mafData[,metaColNames, with=FALSE]

        validObject(.Object)
        return(.Object)
    }
)

#' Constructor for MutationAnnotationFormat_v2.4
#' 
#' @name MutationAnnotationFormat_v2.4
#' @rdname MutationAnnotationFormat_v2.4-class
#' @noRd
MutationAnnotationFormat_v2.4 <- function(mafPath, mafVersion, mafData){
    cat("!!!!! MutationAnnotationFormat_v2.4~Constructor !!!!!\n")
    new("MutationAnnotationFormat_v2.4", mafPath=mafPath, mafVersion=mafVersion,
        mafData=mafData)
}

