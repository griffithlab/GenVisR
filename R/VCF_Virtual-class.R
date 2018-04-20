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

#' @name getVcf
#' @rdname getVcf-methods
#' @aliases getVcf
setMethod(f="getVcf",
          signature="VCF_Virtual",
          definition=function(object, ...){
              vcf <- object@vcf
              return(vcf)
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