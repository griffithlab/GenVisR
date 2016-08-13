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

setMethod(
    f="getPosition",
    signature="MutationAnnotationFormat_Virtual",
    definition=function(object){
        return(object@position)
    }
)

setMethod(
    f="getMutation",
    signature="MutationAnnotationFormat_Virtual",
    definition=function(object){
        return(object@mutation)
    }
)

setMethod(
    f="getSample",
    signature="MutationAnnotationFormat_Virtual",
    definition=function(object){
        return(object@Sample)
    }
)

setMethod(
    f="getMeta",
    signature="MutationAnnotationFormat_Virtual",
    definition=function(object){
        return(object@meta)
    }
)
