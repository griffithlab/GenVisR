#' Class MutationAnnotationFormat_Virtual
#' 
#' An S4 class to act as a virtual class for MutationAnnotationFormat version sub-classes.
#' @name MutationAnnotationFormat_Virtual-class
#' @rdname MutationAnnotationFormat_Virtual-class
#' @slot position data.table object holding genomic positions.
#' @slot mutation data.table object holding mutation status data.
#' @slot sample data.table object holding sample data.
#' @slot meta data.table object holding all other meta data.
#' @importFrom data.table data.table
#' @import methods
setClass(
    Class="MutationAnnotationFormat_Virtual",
    representation=representation(position="data.table",
                                  mutation="data.table",
                                  sample="data.table",
                                  meta="data.table", "VIRTUAL")
)

#' @rdname getPosition-methods
#' @aliases getPosition,MutationAnnotationFormat_Virtual
setMethod(f="getPosition",
          signature="MutationAnnotationFormat_Virtual",
          definition=function(object, ...){
              positions <- object@position
              colnames(positions) <- c("chromosome", "start", "end", "strand")
              return(positions)
          })

#' @rdname getMutation-methods
#' @aliases getMutation,MutationAnnotationFormat_Virtual
setMethod(
    f="getMutation",
    signature="MutationAnnotationFormat_Virtual",
    definition=function(object, ...){
        mutations <- object@mutation
        colnames(mutations) <- c("variant_classification", "variant_type", "reference_allele",
                                 "variant_allele_1", "variant_allele_2")
        return(mutations)
    }
)

#' @rdname getSample-methods
#' @aliases getSample,MutationAnnotationFormat_Virtual
setMethod(
    f="getSample",
    signature="MutationAnnotationFormat_Virtual",
    definition=function(object, ...){
        samples <- object@Sample
        colnames(samples) <- c("sample") 
        return(samples)
    }
)

#' @rdname getMeta-methods
#' @aliases getMeta,MutationAnnotationFormat_Virtual
setMethod(
    f="getMeta",
    signature="MutationAnnotationFormat_Virtual",
    definition=function(object, ...){
        return(object@meta)
    }
)
