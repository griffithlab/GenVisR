################################################################################
######################### Virtual Class Definitions ############################

#' Class VCF_Virtual
#' 
#' An S4 class to act as a virtual class for VCF version sub-classes
#' @name VCF_Virtual
#' @rdname VCF_Virtual-class
#' @slot vcf data.table object holding varscan data
#' @slot sample data.table object holding sample data
#' @importClassesFrom data.table data.table
#' @import methods
#' 
setClass(
    Class="VCF_Virtual",
    representation=representation(description="data.table",
                                  sample="data.table",
                                  vcfData="data.table",
                                  svType="data.table",
                                  "VIRTUAL")
)

################################################################################
###################### Accessor function definitions ###########################

#' @name getHeader
#' @rdname getHeader-methods
#' @aliases getHeader
setMethod(f="getHeader",
          signature="VCF_Virtual",
          definition=function(object, ...){
              header <- object@description
              return(header)
          })

#' @name getSample
#' @rdname getSample-methods
#' @aliases getSample
setMethod(f="getSample",
          signature="VCF_Virtual",
          definition=function(object, ...){
              sample <- object@sample
              return(sample)
          })

#' @rdname getMeta-methods
#' @aliases getMeta
setMethod(f="getMeta",
          signature="VCF_Virtual",
          definition=function(object, ...) {
              meta <- object@vcfData
              return(meta)
          })

#' @rdname getMutation-methods
#' @aliases getMutation
setMethod(f="getMutation",
          signature="VCF_Virtual",
          definition=function(object, ...) {
              mutation <- object@svType
              return(mutation)
          })