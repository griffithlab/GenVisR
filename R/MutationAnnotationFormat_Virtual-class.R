################################################################################
##################### Virutal Class Definition #################################

#' Class MutationAnnotationFormat_Virtual
#' 
#' An S4 class to act as a virtual class for MutationAnnotationFormat version sub-classes.
#' @name MutationAnnotationFormat_Virtual-class
#' @rdname MutationAnnotationFormat_Virtual-class
#' @slot position data.table object holding genomic positions.
#' @slot mutation data.table object holding mutation status data.
#' @slot sample data.table object holding sample data.
#' @slot meta data.table object holding all other meta data.
#' @importClassesFrom data.table data.table
#' @import methods
setClass(
    Class="MutationAnnotationFormat_Virtual",
    representation=representation(position="data.table",
                                  mutation="data.table",
                                  sample="data.table",
                                  meta="data.table", "VIRTUAL")
)

################################################################################
###################### Accessor function definitions ###########################

#' @rdname getPosition-methods
#' @aliases getPosition
setMethod(f="getPosition",
          signature="MutationAnnotationFormat_Virtual",
          definition=function(object, ...){
              positions <- object@position
              return(positions)
          })

#' @rdname getMutation-methods
#' @aliases getMutation
setMethod(f="getMutation",
          signature="MutationAnnotationFormat_Virtual",
          definition=function(object, ...){
              mutations <- object@mutation
              return(mutations)
          })

#' @rdname getSample-methods
#' @aliases getSample
setMethod(f="getSample",
          signature="MutationAnnotationFormat_Virtual",
          definition=function(object, ...){
              sample <- object@sample
              return(sample)
          })

#' @rdname getMeta-methods
#' @aliases getMeta
setMethod(f="getMeta",
          signature="MutationAnnotationFormat_Virtual",
          definition=function(object, ...){
              meta <- object@meta
              return(meta)
          })