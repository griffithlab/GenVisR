################################################################################
##################### Virutal Class Definition #################################

#' Class VarScanFormat_Virtual
#' 
#' An S4 class to act as a virtual class for VarScanFormat version sub-classes.
#' @name VarScanFormat_Virtual-class
#' @rdname VarScanFormat_Virtual-class
#' @slot varscan data.table object holding varscan data.
#' @slot sample data.table object holding sample data.
#' @importClassesFrom data.table data.table
#' @import methods
#' 
setClass(
    Class="VarScanFormat_Virtual",
    representation=representation(varscan="data.table",
                                  sample="data.table",
                                  header="data.table",
                                  "VIRTUAL")
)

################################################################################
###################### Accessor function definitions ###########################

#' @rdname getVarScan-methods
#' @aliases getVarScan
setMethod(f="getVarScan",
          signature="VarScanFormat_Virtual",
          definition=function(object, ...){
              varscan <- object@varscan
              return(varscan)
          })

#' @rdname getSample-methods
#' @aliases getSample
setMethod(f="getSample",
          signature="VarScanFormat_Virtual",
          definition=function(object, ...){
              sample <- object@sample
              return(sample)
          })