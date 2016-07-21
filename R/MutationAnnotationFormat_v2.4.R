#' Class MutationAnnotationFormat_v2.4
#' 
#' An S4 class to represent data in mutation annotation format version 2.4
#' @name MutationAnnotationFormat_v2.4-class
#' @rdname MutationAnnotationFormat_v2.4-class
#' @slot position test
setClass("MutationAnnotationFormat_v2.4",
         contains="MutationAnnotationFormat",
         validity=function(object){
             cat("!!!!! MutationAnnotationFormat_v2.4~Inspector !!!!!\n")
             # expected_col <- c("Hugo_Symbol", "Chromosome", "Start_Position", "End_Position",
             #                   "Variant_Classification", "Reference_Allele", "Tumor_Seq_Allele1",
             #                   "Tumor_Seq_Allele2", "Tumor_Sample_Barcode")
             # if(!all(expected_col %in% colnames(object@data_file))){
             #     memo <- paste("Missing the following required columns:",
             #                   toString(expected_col[!expected_col %in% colnames(object@data_file)]))
             #     stop(memo)
             #}else{}
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
        .Object@path <- mafPath
        .Object@version <- mafVersion
        
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

setMethod(
    f="getPath",
    signature="MutationAnnotationFormat_v2.4",
    definition=function(object){
        return(object@path)
    }
)

setMethod(
    f="getVersion",
    signature="MutationAnnotationFormat_v2.4",
    definition=function(object){
        return(object@version)
    }
)

setMethod(
    f="getPosition",
    signature="MutationAnnotationFormat_v2.4",
    definition=function(object){
        return(object@position)
    }
)

setMethod(
    f="getMutation",
    signature="MutationAnnotationFormat_v2.4",
    definition=function(object){
        return(object@mutation)
    }
)

setMethod(
    f="getSample",
    signature="MutationAnnotationFormat_v2.4",
    definition=function(object){
        return(object@Sample)
    }
)

setMethod(
    f="getMeta",
    signature="MutationAnnotationFormat_v2.4",
    definition=function(object){
        return(object@meta)
    }
)

