#' Class GMS_Virtual
#' 
#' An S4 class to act as a virtual class for GMS version sub-classes.
#' @name GMS_Virtual-class
#' @rdname GMS_Virtual-class
#' @slot position data.table object holding genomic positions.
#' @slot mutation data.table object holding mutation status data.
#' @slot sample data.table object holding sample data.
#' @slot meta data.table object holding all other meta data.
#' @importClassesFrom data.table data.table
#' @import methods
setClass(
    Class="GMS_Virtual",
    representation=representation(position="data.table",
                                  mutation="data.table",
                                  sample="data.table",
                                  meta="data.table", "VIRTUAL")
)
