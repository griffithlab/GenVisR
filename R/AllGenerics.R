#' Method getPosition
#' 
#' @name getPosition
#' @rdname getPosition-methods
#' @return data.table object holding genomic positions and containing column names
#' "chromosome", "start", "end", "strand".
#' @exportMethod getPosition
setGeneric(
    name="getPosition",
    def=function(object, ...){standardGeneric("getPosition", ...)}
)

#' Method getMutation
#' 
#' @name getMutation
#' @rdname getMutation-methods
#' @return data.table object holding genomic mutations and containing column names
#' "variant_classification", "variant_type", "reference_allele", "variant_allele_1", "variant_allele_2".
#' @exportMethod getMutation
setGeneric(
    name="getMutation",
    def=function(object, ...){standardGeneric("getMutation", ...)}
)

#' Method getSample
#' 
#' @name getSample
#' @rdname getSample-methods
#' @return data.table object holding samples and containing column names
#' "sample".
#' @exportMethod getSample
setGeneric(
    name="getSample",
    def=function(object, ...){standardGeneric("getSample", ...)}
)

#' Method getMeta
#' 
#' @name getMeta
#' @rdname getMeta-methods
#' @return data.table object holding meta information.
#' @exportMethod getMeta
setGeneric(
    name="getMeta",
    def=function(object, ...){standardGeneric("getMeta", ...)}
)

#' Method toWaterfall
#'
#' @name toWaterfall
#' @rdname toWaterfall-methods
#' @noRd
#' @return data.table object formatted for input to waterfall
setGeneric(
    name="toWaterfall",
    def=function(object, ...){standardGeneric("toWaterfall", ...)}
)

